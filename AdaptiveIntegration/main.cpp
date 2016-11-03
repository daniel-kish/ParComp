#include <iostream>
#include "threadsafe_stack.h"
#include <future>
#include <chrono>
#include <map>
#include <stack>

//class Opened : public std::stack<range, std::vector<range>>
//{
//	using iter = container_type::iterator;
//
//public:
//	iter begin() {
//		return c.begin();
//	}
//	iter end() {
//		return c.end();
//	}
//};
//
//int n_iters = 0;
//std::map<std::thread::id, int> count;
//using lg = std::lock_guard<std::mutex>;
//std::mutex mGQ;
//
//using localQ = concurrent_list<range>;
//
//std::vector<localQ*> lqs;
//std::mutex m_lqs;
//
//bool done = false;
//std::mutex mio;
//std::condition_variable cv;
//int n_threads = 2;
//
//void est_and_branch1(range r, double& th_local_sum, concurrent_list<range>& lq)
//{
//	double c = (r.a + r.b)*0.5;
//	double fa = f(r.a), fb = f(r.b);
//	double fc = f(c);
//	double s1 = trapeze(fa, fb, r.b - r.a);
//	double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
//	double change = std::abs(s2 - s1);
//
//	if (change < eps)
//		th_local_sum += s2;
//	else {
//		double c = (r.a + r.b)*0.5;
//		double mid_ac = (r.a + c)*0.5;
//		double mid_cb = (c + r.b)*0.5;
//
//		lq.push_front({r.a,mid_ac});
//		lq.push_front({mid_ac,c});
//		lq.push_front({c,mid_cb});
//		lq.push_front({mid_cb,r.b});
//	}
//}
//
//std::map<std::thread::id, int> threads;
//std::map<std::thread::id, int> stealings;
//std::map<std::thread::id, int> poppings;
//std::map<std::thread::id, int> works;
//
//double work(concurrent_list<range>& GQ)
//{
//	localQ lq;
//	std::thread::id myid = std::this_thread::get_id();
//
//	// introduce ourselves
//	{lg lk(m_lqs);  lqs.push_back(&lq); }
//	cv.notify_all();
//
//	// wait for everybody
//	std::unique_lock<std::mutex> lk(m_lqs);
//	cv.wait(lk, [] { return lqs.size() == n_threads; });
//	lk.unlock();
//	{lg lk(mio);  std::cout << myid << ": start\n"; }
//	double th_local_sum = 0.0;
//	
//	while (true)
//	{
//		range r{};
//		if (GQ.try_pop_back(r)) {
//			lq.push_front(r);
//			{lg lk(mio);  std::cout << myid << ": popped\n"; }
//			{lg lk(mio); poppings[myid]++; }
//		}
//		else { // steal
//			lg lk(m_lqs);
//			for (auto p : lqs)
//			{
//				if (p == &lq)
//					continue;
//				if (p->try_pop_back(r,10)) {
//					lq.push_front(r);
//					{lg lk(mio); std::cout << myid << ": stealed\n"; }
//					{lg lk(mio); stealings[myid]++; }
//					break;
//				}
//			}
//		}
//		if (lq.empty())
//			break;
//		{lg lk(mio); std::cout << myid <<": "<< lq.size() << '\n'; }
//		while (true) {
//			range r;
//			if (!lq.try_pop_front(r))
//				break; // we're done locally
//			est_and_branch1(r, th_local_sum, lq);
//			{lg lk(mio); works[myid]++; }
//		}
//	}
//	{
//		lg lk(m_lqs);
//		auto my = std::find(begin(lqs), end(lqs), &lq);
//		lqs.erase(my);
//	}
//	return th_local_sum;
//}
//
//double integral(double left, double right)
//{
//	Opened opened;
//	opened.push({left,right});
//
//	double sum = 0.0;
//
//	while (!opened.empty())
//	{
//		range r = opened.top();
//		opened.pop();
//
//		double c = (r.a + r.b)*0.5;
//		double fa = f(r.a), fb = f(r.b);
//		double fc = f(c);
//		double s1 = trapeze(fa, fb, r.b - r.a);
//		double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
//		double change = std::abs(s2 - s1);
//
//		if (change < eps) 
//			sum += s2;
//		else {
//			double c = (r.a + r.b)*0.5;
//			double mid_ac = (r.a + c)*0.5;
//			double mid_cb = (c + r.b)*0.5;
//
//			opened.push({r.a,mid_ac});
//			opened.push({mid_ac,c});
//			opened.push({c,mid_cb});
//			opened.push({mid_cb,r.b});
//		}
//	}
//	return sum;
//}
//
//double integralPar(double left, double right)
//{
//	concurrent_list<range> opened;
//	int nsteps = 1;
//	double step = (right - left) / nsteps;
//
//	for (double x = left; x < right; x += step)
//		opened.push_back({x,x + step});
//	
//	double sum = 0.0;
//	const int n = n_threads;
//	std::vector<std::future<double>> res(n);
//
//	for (int i = 0; i < n; ++i) {
//		res[i] = std::async(std::launch::async, work, std::ref(opened));
//	}
//	for (int i = 0; i < n; ++i) {
//		sum += res[i].get();
//	}
//
//	return sum;
//}
//
//double integralR(double a, double b)
//{
//	double c = (a + b)*0.5;
//	double fa = f(a), fb = f(b);
//	double fc = f(c);
//	double s1 = trapeze(fa, fb, b - a);
//	double s2 = trapeze(fa, fc, c - a) + trapeze(fc, fb, b - c);
//	double change = std::abs(s2 - s1);
//
//	if (change < eps)      // tree leaf
//		return s2;
//	else {                 // branch
//		return integralR(a, c) + integralR(c, b);
//	}
//}

// stats gathering

std::map<std::thread::id, int> threadUsage;
std::map<std::thread::id, int> tasksDone;

class Pool
{
private:
	struct Range { float a; float b; };
	using prom = std::promise<float>;
	using res = std::future<float>;
	struct Task { Range range; prom promised_value; };
	struct pending_task { Task task; res f1; res f2; };
	
	threadsafe_stack<Task> l;
	std::vector<std::thread> threads;
	std::atomic_bool done;

public:
	Pool(int n_t = std::thread::hardware_concurrency())
		: done{false}
	{
		for (int i = 0; i < n_t; ++i)
			threads.emplace_back(&Pool::worker, this);
	}
	~Pool()
	{
		done = true;
		for (auto& t : threads)
			t.join();
	}
	std::tuple<float, bool> need_splitting(Pool::Range r) {
		auto len = r.b - r.a;
		if (len <= 2.5)	return{len,false};
		return{len,true};
	}
	std::pair<res, res> split(Task& t) 
	{
		float m = (t.range.a + t.range.b) / 2.0f;
		Task c1{{t.range.a, m}, prom{}};
		Task c2{{m, t.range.b}, prom{}};
		auto f1 = c1.promised_value.get_future();
		auto f2 = c2.promised_value.get_future();
		l.push_back(std::move(c1));
		l.push_back(std::move(c2));
		return{std::move(f1), std::move(f2)};
	}
	bool ready(std::future<float>& f)
	{
		return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
	}
	void worker()
	{
		using std::move; using std::tie;
		auto myid = std::this_thread::get_id();
		std::stack<pending_task> pend;
		while (!done)
		{
			if (!pend.empty() && ready(pend.top().f1) && ready(pend.top().f2)) {
				pend.top().task.promised_value.set_value(pend.top().f1.get() + pend.top().f2.get());
				pend.pop();
			}
			pending_task pt;
			if (!l.try_pop_back(pt.task)) continue;

			float len; bool tosplit;
			tie(len, tosplit) = need_splitting(pt.task.range);
			if (!tosplit) {
				pt.task.promised_value.set_value(len);
				continue;
			}

			tie(pt.f1, pt.f2) = split(pt.task);
			pend.push(move(pt));
		}
	}
	std::future<float> submit(Pool::Range r)
	{
		Task t{r, prom{}};
		auto f = t.promised_value.get_future();
		l.push_back(std::move(t));
		return f;
	}
};

template <class Fun>
std::chrono::duration<double, std::milli> time_stats(Fun f, int times = 10)
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

	int div = durs.size() / 2;
	if (durs.size() % 2 == 0)
		return 0.5*(durs[div - 1] + durs[div]);
	else
		return durs[div];
}

int main()
{
	using namespace std::chrono;

	Pool p(4);
	float len{};

	auto work = [&p, &len] {
		auto f = p.submit({1,137500});
		len = f.get();
	};
	std::cout << time_stats(work, 1).count() << '\n';
	std::cout << len << '\n';
}