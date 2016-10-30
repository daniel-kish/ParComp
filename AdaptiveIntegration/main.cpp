#include <iostream>
#include <future>
#include <iomanip>
#include <map>
#include <queue>
#include <stack>
#include <list>
#include <chrono>
#include "fun.h"




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
using lg = std::lock_guard<std::mutex>;
std::mutex mGQ;

using localQ = std::list<range>;
std::vector<localQ*> lqs;
std::mutex m_lqs;
bool done = false;
std::mutex mio;
std::condition_variable cv;

double work(Opened& GQ)
{
	{lg lk(mio); std::cout << std::this_thread::get_id() << '\n'; }
	localQ lq; // local queue
	{
		lg lk(m_lqs);
		{lg lk(mio); std::cout << std::this_thread::get_id() << " blocks 'lqs'\n"; }
		lqs.push_back(&lq);
	}
	{lg lk(mio); std::cout << std::this_thread::get_id() << " starts waiting\n"; }
	std::unique_lock<std::mutex> lk(m_lqs);
	cv.wait(lk, [] { return lqs.size() == 2; });
	cv.notify_all();

	{lg lk(mio); std::cout << std::this_thread::get_id() << " done waiting\n"; }

	double th_local_sum = 0.0;
	
	while (!done) 
	{
		{lg lk(mio); std::cout << std::this_thread::get_id() << " enters main loop\n"; }
		{
			lg lk(mGQ);
			{lg lk(mio); std::cout << std::this_thread::get_id() << " blocks GQ\n"; }
			if (!GQ.empty()) {
				lq.push_back(GQ.top());
				GQ.pop();
				{lg lk(mio); std::cout << std::this_thread::get_id() << " pops from GQ\n"; }
			}
		}
		if (lq.empty())
		{
			{lg lk(mio); std::cout << std::this_thread::get_id() << " tries stealing\n"; }
			lg lk(m_lqs);
			for (auto p : lqs)
			{
				localQ& q = *p;
				if (!q.empty()) {
					{lg lk(mio); std::cout << std::this_thread::get_id() << " steals\n"; }
					lq.push_back(q.back());
					q.pop_back();
					break;
				}
			}
		}
		if (lq.empty()) {
			done = true;
			continue;
		}

		int iter = 0;
		while (!lq.empty())
		{
			iter++;
			range r = lq.front();
			lq.pop_front();

			double c = (r.a + r.b)*0.5;
			double fa = f(r.a), fb = f(r.b);
			double fc = f(c);
			double s1 = trapeze(fa, fb, r.b - r.a);
			double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
			double change = std::abs(s2 - s1);

			if (change < eps)
				th_local_sum += s2;
			else {
				double c = (r.a + r.b)*0.5;
				double mid_ac = (r.a + c)*0.5;
				double mid_cb = (c + r.b)*0.5;

				lq.push_back({r.a,mid_ac});
				lq.push_back({mid_ac,c});
				lq.push_back({c,mid_cb});
				lq.push_back({mid_cb,r.b});
			}
		}
	}
	return th_local_sum;
}

double integral(double left, double right)
{
	Opened opened;
	opened.push({left,right});

	double sum = 0.0;

	while (!opened.empty())
	{
		range r = opened.top();
		opened.pop();

		double c = (r.a + r.b)*0.5;
		double fa = f(r.a), fb = f(r.b);
		double fc = f(c);
		double s1 = trapeze(fa, fb, r.b - r.a);
		double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
		double change = std::abs(s2 - s1);

		if (change < eps) 
			sum += s2;
		else {
			double c = (r.a + r.b)*0.5;
			double mid_ac = (r.a + c)*0.5;
			double mid_cb = (c + r.b)*0.5;

			opened.push({r.a,mid_ac});
			opened.push({mid_ac,c});
			opened.push({c,mid_cb});
			opened.push({mid_cb,r.b});
		}
	}
	return sum;
}

double integralPar(double left, double right)
{
	Opened opened;
	int nsteps=8;
	double step = (right - left) / nsteps;

	for (double x = left; x < right; x += step)
		opened.push({x,x + step});


	double sum = 0.0;
	const int n_threads = 2;
	std::vector<std::future<double>> res(n_threads);

	for (int i = 0; i < n_threads; ++i) {
		res[i] = std::async(std::launch::async, work, std::ref(opened));
	}
	for (int i = 0; i < n_threads; ++i) {
		sum += res[i].get();
	}

	return sum;
}

double integralR(double a, double b)
{
	double c = (a + b)*0.5;
	double fa = f(a), fb = f(b);
	double fc = f(c);
	double s1 = trapeze(fa, fb, b - a);
	double s2 = trapeze(fa, fc, c - a) + trapeze(fc, fb, b - c);
	double change = std::abs(s2 - s1);

	if (change < eps)      // tree leaf
		return s2;
	else {                 // branch
		return integralR(a, c) + integralR(c, b);
	}
}

int main()
{
	using namespace std::chrono;

	int times = 1;
	std::vector<double> r(times);
	auto t1 = high_resolution_clock::now();

	for (int i = 0; i < times; i++)
		r[i] = integralPar(0, 1);

	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> dur(t2 - t1);

	std::cout << std::setprecision(16) << std::fixed;
	for (auto v : r)
		std::cout << v << ' ';
	std::cout << '\n';

	std::cout << dur.count() / times << '\n';
	std::cout << n_iters << '\n';
	
	for (auto& p : count)
		std::cout << p.first << " / " << p.second << '\n';
	std::cout << count.size() << '\n';

	std::cout << eps << '\n';

	/*using namespace std::chrono;

	int times = 100;
	std::vector<double> results;
	double integ=0.0;
	auto t1 = high_resolution_clock::now();

	for(int i=0; i < times; i++)
		integ = integral(0.05, 1.0);

	auto t2 = high_resolution_clock::now();
	results.push_back(integ);
	duration<double, std::micro> d(t2 - t1);

	//std::cout << std::setprecision(16) << std::fixed << integ << '\n';
	std::cout << "time: " << d.count() / double(times) << '\n';

	//std::cout << level << '\n';
	//for (auto& p : count)
	//	std::cout << p.first << " / " << p.second << '\n';
	std::cout << count.size() << '\n';
	*/
}