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
#include <typeinfo>
#include <type_traits>

struct waited_res {
	bool ready;
	double res;
};

struct Range
{
	Point a, b;
	double eps;
	
	waited_res ans1;
	waited_res ans2;
	waited_res *res;

	Range(Point const& a_, Point const& b_, double eps_, waited_res* par)
		: a{a_}, b{b_}, eps{eps_}, res{par}, ans1{0.0, false}, ans2{0.0,false}
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
//std::mutex mq;
//std::mutex iom;
//using lg = std::lock_guard<std::mutex>;
//
//double work(Opened& GQ)
//{
//	Opened lq; // local queue
//	double th_local_sum = 0.0;
//	bool done = false;
//
//	while (!done)
//	{
//		{ // critical section
//			std::lock_guard<std::mutex> lk(mq);
//			if (!GQ.empty()) {
//				lq.push(GQ.top());
//				GQ.pop();
//			}
//			else {
//				break;
//			}
//		}
//
//		//{
//		//	lg lk(iom);
//		//	count[std::this_thread::get_id()]++;
//		//	//std::cout << std::this_thread::get_id() << ' ' << lq.top() << '\n';
//		//}
//		int iter = 0;
//		while (!lq.empty())
//		{
//			iter++;
//			range r = lq.top();
//			lq.pop();
//
//			double c = (r.a + r.b)*0.5;
//			double fa = f(r.a), fb = f(r.b);
//			double fc = f(c);
//			double s1 = trapeze(fa, fb, r.b - r.a);
//			double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
//			double change = std::abs(s2 - s1);
//
//			if (change < r.eps)
//				th_local_sum += s2;
//			else {
//				double c = (r.a + r.b)*0.5;
//				double mid_ac = (r.a + c)*0.5;
//				double mid_cb = (c + r.b)*0.5;
//
//				lq.push({r.a,mid_ac, r.eps*0.25});
//				lq.push({mid_ac,c, r.eps*0.25});
//				lq.push({c,mid_cb, r.eps*0.25});
//				lq.push({mid_cb,r.b, r.eps*0.25});
//			}
//		}
//	}
//	return th_local_sum;
//}
//
//double integral(double left, double right)
//{
//	Opened opened;
//	opened.push({left,right, eps});
//
//	double sum = 0.0;
//	int i = 0;
//	while (!opened.empty())
//	{
//		range r = opened.top();
//		opened.pop();
//
//		Point c = middle(r);
//		double s1 = trapeze(r.a.f(), r.b.f(), length(r));
//		double s2 = trapeze(r.a.f(), c.f(), dist(r.a,c)) + trapeze(c.f(), r.b.f(), dist(c,r.b));
//		double change = std::abs(s2 - s1);
//
//		if (change < r.eps)
//			sum += s2;
//		else {
//			opened.push({c, r.b, r.eps * 0.5});
//			opened.push({r.a, c, r.eps * 0.5});
//			//double mac = (r.a + c)*0.5;
//			//double mcb = (c + r.b)*0.5;
//			//opened.push({r.a, mac, r.eps * 0.25});
//			//opened.push({mac, c,   r.eps * 0.25});
//			//opened.push({c, mcb,   r.eps * 0.25});
//			//opened.push({mcb, r.b, r.eps * 0.25});
//		}
//	}
//	return sum;
//}
//
//double integralPar(double left, double right, int nt=1, int chunks_per_thread=1)
//{
//	Opened opened;
//	int chunks = nt*chunks_per_thread;
//	double step = (right - left) / chunks;
//	for (int i = 0; i < chunks; i++)
//		opened.push({left + i*step, left + (i + 1)*step, eps / chunks});
//
//	double sum = 0.0;
//	const int n_threads = nt;
//	std::vector<std::future<double>> res(n_threads);
//
//	for (int i = 0; i < n_threads; ++i) {
//		res[i] = std::async(std::launch::async, work, std::ref(opened));
//	}
//	for (int i = 0; i < n_threads; ++i) {
//		sum += res[i].get();
//	}
//
//	return sum;
//}
//
//double integralR(double left, double right, double myeps = eps)
//{
//	Point a(), b;
//	double fa = f(a), fb = f(b);
//	double fc = f(c);
//	double s1 = trapeze(fa, fb, b - a);
//	double s2 = trapeze(fa, fc, c - a) + trapeze(fc, fb, b - c);
//	double change = abs(s2 - s1);
//
//	if (change < myeps)      // tree leaf
//		return s2;
//	else                     // branch
//		return integralR(a, c, myeps*0.5) + integralR(c, b, myeps*0.5);
//}
//
//template <class Fun>
//std::chrono::duration<double, std::milli> time_stats(Fun f, int times = 10)
//{
//	using namespace std::chrono;
//
//	std::vector<duration<double, std::milli>> durs;
//	for (int i = 0; i < times; ++i)
//	{
//		auto t1 = high_resolution_clock::now();
//		f();
//		auto t2 = high_resolution_clock::now();
//		duration<double, std::milli> dur(t2 - t1);
//		durs.push_back(dur);
//	}
//	std::sort(begin(durs), end(durs));
//
//	for (auto& d : durs)
//		std::cout << d.count() << '\n';
//
//	int div = durs.size() / 2;
//	if (durs.size() % 2 == 0)
//		return 0.5*(durs[div - 1] + durs[div]);
//	else {
//		return durs[div];
//	}
//}


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
try{
	waited_res s{0.0,false};
	std::vector<Range> tasks;
	tasks.reserve(140);

	tasks.emplace_back(left,right,eps,&s);
	while (!tasks.empty()) 
	{
		if (tasks.back().ans1.ready && tasks.back().ans2.ready) {
			tasks.back().res->ready = true;
			tasks.back().res->res = tasks.back().ans1.res + tasks.back().ans2.res;
			tasks.pop_back();
			continue;
		}

		int top_id = tasks.size() - 1;
		Range& r = tasks[top_id];
		
		double s1 = trapezium_area(r.a, r.b);
		Point c = middle(r.a, r.b);
		double s2 = trapezium_area(r.a, c) + trapezium_area(c, r.b);
		double change = abs(s2 - s1);

		if (change < r.eps) {
			r.res->ready = true;
			r.res->res = s2;
			tasks.pop_back();
			continue;
		}
		tasks.emplace_back(r.a, c, r.eps / 2, &r.ans1);
		tasks.emplace_back(c, r.b, r.eps / 2, &r.ans2);
	}
	return s.res;
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

	int div = durs.size() / 2;
	if (durs.size() % 2 == 0)
		return 0.5*(durs[div - 1] + durs[div]);
	else
		return durs[div];
}

int main()
{
	using namespace std::chrono;
	const double Eps = 1.0e-15;
	double a = 0.0, b = M_PI; double exact = 2.0;

	std::cout << "Stack based\n";
	auto t1 = high_resolution_clock::now();
	double s = integralStack(a, b, Eps);
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> dur = (t2 - t1);
	
	std::cout << "ans = " << std::setprecision(-log10(Eps)) << std::fixed << abs(s-exact) << '\n';
	std::cout << "elapsed " << dur.count() << " ms\n";
	double d1 = dur.count();

	std::cout << "\nRecursive\n";
	t1 = high_resolution_clock::now();
	s = integralR(a, b, Eps);
	t2 = high_resolution_clock::now();
	dur = t2 - t1;
	
	std::cout << "ans = " << std::setprecision(-log10(Eps)) << std::fixed << abs(s - exact) << '\n';
	std::cout << "elapsed " << dur.count() << " ms\n";
	
	std::cout << "Stack to Recursive ratio: " << std::setprecision(2) << d1 / dur.count() << '\n';
	
	//std::cout << "\nxs: " << xs.size() << '\n';
	//std::sort(begin(xs), end(xs));
	//auto p = std::unique(begin(xs), end(xs));
	//if (p != end(xs))
	//	std::cout << "!!!!!!!!!!!!!!\n";
	//std::copy(begin(xs), end(xs), std::ostream_iterator<double>{std::cout, "\n"});
}