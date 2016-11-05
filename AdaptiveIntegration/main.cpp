#include <iostream>
#include <future>
#include <iomanip>
#include <map>
#include <queue>
#include <stack>
#include <list>
#include <chrono>
#include <random>
#include "Point.h"

const double eps = 1.0e-10;
const double Pi = 4.0*atan(1.0);





struct range
{ 
	Point a, b;
	double eps;
};


Point middle(range const& r)
{
	Point p{ (r.a.x() + r.b.x())*0.5 };
	return p;
}

Point middle(Point const& a, Point const& b)
{
	Point p{ (a.x() + b.x())*0.5 };
	return p;
}

double length(range const& r)
{
	return dist(r.a,r.b);
}

std::ostream& operator<< (std::ostream& os, range const& r)
{
	os << '(' << r.a << ' ' << r.b << ") " << r.eps;
	return os;
}

class Opened : public std::stack<range, std::vector<range>>
{
	using iter = container_type::iterator;

public:
	iter begin() {
		return c.begin();
	}
	iter end() {
		return c.end();
	}
};

int n_iters = 0;
std::map<std::thread::id, int> count;
std::mutex mq;
std::mutex iom;
using lg = std::lock_guard<std::mutex>;

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

double integral(double left, double right)
{
	Opened opened;
	opened.push({left,right, eps});

	double sum = 0.0;
	int i = 0;
	while (!opened.empty())
	{
		range r = opened.top();
		opened.pop();

		Point c = middle(r);
		double s1 = trapeze(r.a.f(), r.b.f(), length(r));
		double s2 = trapeze(r.a.f(), c.f(), dist(r.a,c)) + trapeze(c.f(), r.b.f(), dist(c,r.b));
		double change = std::abs(s2 - s1);

		if (change < r.eps)
			sum += s2;
		else {
			opened.push({c, r.b, r.eps * 0.5});
			opened.push({r.a, c, r.eps * 0.5});
			//double mac = (r.a + c)*0.5;
			//double mcb = (c + r.b)*0.5;
			//opened.push({r.a, mac, r.eps * 0.25});
			//opened.push({mac, c,   r.eps * 0.25});
			//opened.push({c, mcb,   r.eps * 0.25});
			//opened.push({mcb, r.b, r.eps * 0.25});
		}
	}
	return sum;
}

double integralPar(double left, double right, int nt=1, int chunks_per_thread=1)
{
	Opened opened;
	int chunks = nt*chunks_per_thread;
	double step = (right - left) / chunks;
	for (int i = 0; i < chunks; i++)
		opened.push({left + i*step, left + (i + 1)*step, eps / chunks});

	double sum = 0.0;
	const int n_threads = nt;
	std::vector<std::future<double>> res(n_threads);

	for (int i = 0; i < n_threads; ++i) {
		res[i] = std::async(std::launch::async, work, std::ref(opened));
	}
	for (int i = 0; i < n_threads; ++i) {
		sum += res[i].get();
	}

	return sum;
}

double integralR(double a, double b, double myeps = eps)
{
	double c = (a + b)*0.5;
	double fa = f(a), fb = f(b);
	double fc = f(c);
	double s1 = trapeze(fa, fb, b - a);
	double s2 = trapeze(fa, fc, c - a) + trapeze(fc, fb, b - c);
	double change = abs(s2 - s1);

	if (change < myeps)      // tree leaf
		return s2;
	else                     // branch
		return integralR(a, c, myeps*0.5) + integralR(c, b, myeps*0.5);
}

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

	for (auto& d : durs)
		std::cout << d.count() << '\n';

	int div = durs.size() / 2;
	if (durs.size() % 2 == 0)
		return 0.5*(durs[div - 1] + durs[div]);
	else {
		return durs[div];
	}
}

int main()
{
	using namespace std::chrono;

	//double s = integralPar(0.01, 1, 1, 10);
	//std::cout << std::setprecision(std::ceil(abs(log10(eps)))+1) << std::fixed << s << '\n';
	auto&& digits = std::setprecision(abs(log10(eps)));
	
	int times = 1;
	std::vector<double> res(times);
	int i = 0;
	auto work = [&res,&i] {
		res[i++] = integralR(0.01, 1);
	};

	std::cout << '\n' << time_stats(work, times).count() << '\n';
	for (auto d : res)
		std::cout << digits << std::fixed << abs(d - 0.503981893175415) << '\n';
}