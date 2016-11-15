#include <iostream>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include "math.h"
#include "Point.h"
#include <chrono>
#include <stack>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <list>
#include <future>
#include <map>
#include "threadsafe_queue.h"
#include <sstream>
#include <cassert>
#include <random>

using namespace std::chrono;
using namespace std::chrono_literals;

struct waited_res {
	bool ready;
	double value;
};
enum State { childless, awaiting, used };
struct Range
{
	Point a, b;
	double eps;

	waited_res left_child;
	waited_res right_child;
	waited_res *res;
	State state;

	Range(Point const& a_, Point const& b_, double eps_, waited_res* par)
		: a{a_}, b{b_}, eps{eps_}, res{par}, left_child{false, -1.0}, right_child{false, 1.0},
		state{childless}
	{	}
	Range()
	{	}
};

Point middle(Range const& r)
{
	return middle(r.a, r.b);
}

double length(Range const& r)
{
	return dist(r.a, r.b);
}

std::ostream& operator<< (std::ostream& os, Range const& r)
{
	os << '(' << r.a << ' ' << r.b << ") " << r.eps;
	return os;
}

double integralStack(double left, double right, double eps);

using guard = std::lock_guard<std::mutex>;
struct chunk {
	double a, b;
	double eps;
};
std::mutex mio;
std::map<std::thread::id, int> usage;
std::map<std::thread::id, int> localPops;
std::map<std::thread::id, int> localPushes;
std::map<std::thread::id, int> globalPops;
std::map<std::thread::id, int> globalPushes;
struct pendingParams { int min; int max; std::size_t times; std::size_t sum; };

template <class T>
void rightRotate(std::list<T>& l)
{
	std::list<T> tmp;
	auto from = std::prev(l.end());
	tmp.splice(tmp.begin(), l, from, l.end());
	auto to = l.begin();
	l.splice(to, tmp);
}
class Counter {
	int cnt;
	int limit;
public:
	Counter (int lim) : limit(lim)
	{}
	void operator++() {
		cnt++;
		if (cnt == limit)
			cnt = 0;
	}
	operator int () {
		return cnt;
	}
};

struct Pool
{
	std::vector<Range> globalQ;
	std::mutex globalMx;
	std::atomic_bool done;
	std::vector<std::thread> threads;
	int max_local_size;

	using Pending = std::vector<Range>;
	using iter = Pending::iterator;
	static thread_local std::unique_ptr<Pending> pending;
	static thread_local iter top;
	static thread_local int idx;
	static thread_local int doGlobalPush;
	std::vector<bool> free;

	Pool(int nth, int m = 2000) : done{false}, max_local_size{m}
	{
		free.resize(nth, true);
		//globalQ.reserve(50);
		for (int i = 0; i < nth; ++i)
			threads.emplace_back(&Pool::worker, this, i);
	}
	void shutdown()
	{
		done = true;
		for (auto& t : threads)
			t.join();
	}
	~Pool()
	{
		if (std::none_of(begin(threads), end(threads), [](auto& t) {return t.joinable(); }))
			return;
		shutdown();
	}
	bool globalPop(Range& r)
	{
		guard l(globalMx);
		if (!globalQ.empty()) {
			r = globalQ.back();
			globalQ.pop_back();
			return true;
		}
		return false;
	}
	std::tuple<double, double, Point> estimation(Range& r)
	{
		const double s1 = trapezium_area(r.a, r.b);
		const Point c = middle(r.a, r.b);
		const double s2 = trapezium_area(r.a, c) + trapezium_area(c, r.b);
		const double change = abs(s2 - s1);
		return{change, s2, c};
	}
	bool both_ready(Range& r)
	{
		return r.left_child.ready && r.right_child.ready;
	}

	bool find_free()
	{
		bool some_free{false};
		for (int i = 0; i < free.size(); ++i)
		{
			if (i != idx && free[i]) {
				some_free = true;
				break;
			}
		}
		return some_free;
	}
	void push(Range & r, Point const& c)
	{
		if (find_free() && doGlobalPush)
		{
			doGlobalPush--;
			pending->emplace_back(r.a, c, r.eps / 2, &r.left_child);
			guard l(globalMx);
			globalQ.emplace_back(c, r.b, r.eps / 2, &r.right_child);
		}
		else {
		pending->emplace_back(r.a, c, r.eps / 2, &r.left_child);
		pending->emplace_back(c, r.b, r.eps / 2, &r.right_child);
		}
	}
	void worker(int idx_)
	{
		idx = idx_;
		auto myid = std::this_thread::get_id();
		pending.reset(new Pending);
		pending->reserve(10'000);
		top = pending->begin();
		doGlobalPush = 5000;
		int max_s{0};
		while (!done)
		{
			if (pending->empty()) free[idx] = true;
			if (free[idx]) {
				Range r;
				if (globalPop(r)) {
					pending->push_back(r);
					free[idx] = false;
					//{guard l(mio); std::cout << idx << " globalPop " << r.a.x() << ' ' << r.b.x() << '\n'; }
				}
				else {
					free[idx] = true;
					std::this_thread::yield();
					continue;
				}
				top = pending->end() - 1;
			}
			if (pending->size() > max_s)
				max_s = pending->size();
			Range& r = *top;
			if (r.state == childless)
			{
				double change, est; Point c;
				std::tie(change, est, c) = estimation(r);
				if (change < r.eps) {
					r.res->value = est;
					r.res->ready = true;
					r.state = used;
					continue;
				}
				r.state = awaiting;
				push(r, c);
				top = pending->end() - 1;
			}
			else if (both_ready(r) && r.state != used) {
				r.res->value = r.left_child.value + r.right_child.value;
				r.res->ready = true;
				r.state = used;
				continue;
			}
			else if (r.state == awaiting) {
				--top;
				if (top < pending->begin()) {
					Range r;
					if (globalPop(r)) pending->push_back(r);
					top = pending->end() - 1;
				}
				continue;
			}
			else if (r.state == used) {
				if (top == pending->end() - 1) { // 'used' on top -> delete
					pending->pop_back();
					top = pending->end() - 1;
					continue;
				}
				--top;
				if (top < pending->begin()) {
					Range r;
					if (globalPop(r)) pending->push_back(r);
					top = pending->end() - 1;
				}
			}
		}
		//{guard l(mio); std::cout << pending->capacity() <<' '<< max_s <<' '<< doGlobalPush << '\n'; }
	}
	void submit(double a, double b, double eps, waited_res* promise)
	{
		Range r(a, b, eps, promise);
		guard l(globalMx);
		globalQ.push_back(r);
	}
};

thread_local std::unique_ptr<Pool::Pending> Pool::pending;
thread_local Pool::iter Pool::top;
thread_local int Pool::idx;
thread_local int Pool::doGlobalPush;

double work(std::vector<chunk>& chunks, std::mutex& m_chunks)
{
	auto myid = std::this_thread::get_id();
	double th_local_sum = 0.0;
	std::vector<Range> tasks;
	tasks.reserve(140);

	while (true)
	{
		chunk c;
		{
			guard g(m_chunks);
			if (chunks.empty())
				break;
			c = chunks.back();
			chunks.pop_back();
		}
		waited_res s{0.0,false};
		tasks.emplace_back(c.a, c.b, c.eps, &s);
		while (!tasks.empty())
		{
			if (tasks.back().left_child.ready && tasks.back().right_child.ready) {
				tasks.back().res->ready = true;
				tasks.back().res->value = tasks.back().left_child.value + tasks.back().right_child.value;
				tasks.pop_back();
				continue;
			}

			const int top_id = tasks.size() - 1;
			Range& r = tasks[top_id];

			const double s1 = trapezium_area(r.a, r.b);
			const Point c = middle(r.a, r.b);
			const double s2 = trapezium_area(r.a, c) + trapezium_area(c, r.b);
			const double change = abs(s2 - s1);

			if (change < r.eps) {
				r.res->ready = true;
				r.res->value = s2;
				tasks.pop_back();
				continue;
			}
			tasks.emplace_back(r.a, c, r.eps / 2, &r.left_child);
			tasks.emplace_back(c, r.b, r.eps / 2, &r.right_child);
		}
		th_local_sum += s.value;
	}
	return th_local_sum;
}

double integralPar(double left, double right, double eps, int n_threads = 1, int chunks_per_thread = 1)
{
	int n_chunks = n_threads*chunks_per_thread;
	std::vector<chunk> chunks;
	chunks.reserve(n_chunks);

	double step = (right - left) / n_chunks;
	for (int i = 0; i < n_chunks; i++)
	{
		chunks.push_back({left + i*step, left + (i + 1)*step, eps / n_chunks});
	}

	double sum = 0.0;
	std::vector<std::future<double>> res(n_threads);
	std::mutex m_chunks;
	for (int i = 0; i < n_threads; ++i) {
		res[i] = std::async(std::launch::async, work, std::ref(chunks), std::ref(m_chunks));
	}
	for (int i = 0; i < n_threads; ++i) {
		sum += res[i].get();
	}

	return sum;
}

double integralRImpl(Point a, Point b, double myeps)
{
	double s1 = trapezium_area(a, b);
	Point c = middle(a, b);
	double s2 = trapezium_area(a, c) + trapezium_area(c, b);
	double change = abs(s2 - s1);

	if (change < myeps)
		return s2;
	else
		return integralRImpl(a, c, myeps / 2) + integralRImpl(c, b, myeps / 2);
}

double integralR(double left, double right, double myeps)
{
	Point a{left}, b{right};
	return integralRImpl(a, b, myeps);
}

double integralStack(double left, double right, double eps)
try {
	waited_res s{0.0,false};
	std::vector<Range> tasks;
	tasks.reserve(140);

	tasks.emplace_back(left, right, eps, &s);
	while (!tasks.empty())
	{
		if (tasks.back().left_child.ready && tasks.back().right_child.ready) {
			tasks.back().res->ready = true;
			tasks.back().res->value = tasks.back().left_child.value + tasks.back().right_child.value;
			tasks.pop_back();
			continue;
		}

		const int top_id = tasks.size() - 1;
		Range& r = tasks[top_id];

		const double s1 = trapezium_area(r.a, r.b);
		const Point c = middle(r.a, r.b);
		const double s2 = trapezium_area(r.a, c) + trapezium_area(c, r.b);
		const double change = abs(s2 - s1);

		if (change < r.eps) {
			r.res->ready = true;
			r.res->value = s2;
			tasks.pop_back();
			continue;
		}
		tasks.emplace_back(r.a, c, r.eps / 2, &r.left_child);
		tasks.emplace_back(c, r.b, r.eps / 2, &r.right_child);
	}
	return s.value;
}
catch (std::exception& e)
{
	auto* p = dynamic_cast<std::bad_alloc*>(&e);
	if (p)
		std::cerr << "memory allocation error\n";
	else
		std::cerr << e.what() << '\n';
}

template <class Fun>
std::chrono::duration<double, std::milli> time_stats(Fun& f, int times = 10)
{
	using namespace std::chrono;

	std::vector<duration<double, std::milli>> durs;
	for (int i = 0; i < times; ++i)
	{
		auto t1 = high_resolution_clock::now();
		f();
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> dur(t2 - t1);
		durs.push_back(dur);
	}
	std::sort(begin(durs), end(durs));

	for (auto d : durs)
		std::cout << d.count() << '\n';

	int div = durs.size() / 2;
	if (durs.size() % 2 == 0)
		return 0.5*(durs[div - 1] + durs[div]);
	else
		return durs[div];
}


int main()
{
	const double Eps = 1.0e-12;
	double a = 0.01, b = 1.0; double exact = 0.503981893175415; double S = 0.0;

	Pool p(4);
	waited_res r{false, 0.0};
	auto t1 = high_resolution_clock::now();

	p.submit(a, b, Eps, &r);
	while (!r.ready)
		std::this_thread::yield();
	//S = integralPar(a, b, Eps, 1, 1);
	//S = integralR(a, b, Eps);
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> dur = (t2 - t1);
	{ guard l(mio);
	std::cout << "ans = " << std::setprecision(-log10(Eps)) << std::fixed
		<< S << ' ' << abs(S - exact) << '\n';
 	std::cout << "elapsed " << dur.count  () << " ms\n\n";
	}

	//std::cout << "globalPops\n";
	//for (auto& p : globalPops)
	//	std::cout << p.first << ' ' << p.second << '\n';
	//
	//std::cout << "\nglobalPushes\n";
	//for (auto& p : globalPushes)
	//	std::cout << p.first << ' ' << p.second << '\n';
	//
	//std::cout << '\n';
	//p.shutdown();
	//
	//for (auto& p : p.pendingUsage) 
	//{
	//	std::cout << p.first << ": \n";
	//	std::cout << "\tmax " << p.second.max << '\n';
	//	std::cout << "\tave " << double(p.second.sum) / double(p.second.times) << '\n';
	//}
}


// benchmark
//const double Eps = 1.0e-12;
//double a = 0, b = M_PI; double exact = 2.0;
//
//duration<double, std::milli> min_dur(5h);
//int min_i, min_j;
//duration<double, std::milli> max_dur(0s);
//int max_i, max_j;
//
//const int max_nth = std::thread::hardware_concurrency();
//
//for (int nth = 1; nth <= max_nth; nth *= 2) {
//	for (int ch_per_th = 1; ch_per_th <= 32; ch_per_th *= 2) {
//		auto t1 = high_resolution_clock::now();
//		double s = integralPar(a, b, Eps, nth, ch_per_th);
//		auto t2 = high_resolution_clock::now();
//		duration<double, std::milli> dur = (t2 - t1);
//
//		if (dur < min_dur) {
//			min_dur = dur;
//			min_i = nth;
//			min_j = ch_per_th;
//		}
//		if (dur > max_dur) {
//			max_dur = dur;
//			max_i = nth;
//			max_j = ch_per_th;
//		}
//
//		std::cout << nth << ' ' << ch_per_th << ":\n";
//		std::cout << "ans = " << std::setprecision(-log10(Eps)) << std::fixed << abs(s - exact) << '\n';
//		std::cout << "elapsed " << dur.count() << " ms\n\n";
//	}
//}
//
//std::cout << min_dur.count() << ' ' << min_i << ' ' << min_j << '\n';
//std::cout << max_dur.count() << ' ' << max_i << ' ' << max_j << '\n';
//std::cout << "Smax = " << max_dur.count()/min_dur.count() << '\n';
//
//
//for (auto p : usage)
//	std::cout << p.first << ' ' << p.second << '\n';

//Pool p(1);
//int n = 1;
//double dx = (b - a) / n;
//std::vector<waited_res> results(n, {false,0.0});
//for (int i = 0; i < n; ++i)
//	p.submit(a + i*dx, a + (i + 1)*dx, Eps / n, &results[i]);
//
//for (int i = 0; i < n; ++i)
//	while (!results[i].ready)
//		std::this_thread::yield();
//
//for (auto& r : results)
//S += r.value;


//std::cout << "\nxs: " << xs.size() << '\n';
//std::sort(begin(xs), end(xs));
//auto p = std::unique(begin(xs), end(xs));
//if (p != end(xs))
//	std::cout << "!!!!!!!!!!!!!!\n";
//std::copy(begin(xs), end(xs), std::ostream_iterator<double>{std::cout, "\n"});