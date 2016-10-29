#include <iostream>
#include <future>
#include <iomanip>
#include <map>
#include <queue>
#include <stack>
#include <list>
#include <chrono>
#include <random>

const double eps = 0.000000000001;
const double Pi = 4.0*atan(1.0);

inline double sq(double x) { return x*x; }

double f(double x)
{
	return sq(sin(Pi*x)) + sq(sin(2.0*Pi*x)) + cos(200.0*x*x*x*x);
}

double trapeze(double fa, double fb, double step)
{
	return (fa + fb)*0.5*step;
}

struct range { double a, b; };
std::ostream& operator<< (std::ostream& os, range const& r)
{
	os << '(' << r.a << ' ' << r.b << ')';
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

double work(Opened& GQ)
{
	Opened lq; // local queue
	double th_local_sum = 0.0;
	bool done = false;

	while (!done)
	{
		{ // critical section
			std::lock_guard<std::mutex> lk(mq);
			if (!GQ.empty()) {
				lq.push(GQ.top());
				GQ.pop();
			}
			else {
				break;
			}
		}

		//{
		//	lg lk(iom);
		//	count[std::this_thread::get_id()]++;
		//	//std::cout << std::this_thread::get_id() << ' ' << lq.top() << '\n';
		//}
		int iter = 0;
		while (!lq.empty())
		{
			iter++;
			range r = lq.top();
			lq.pop();

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

				lq.push({r.a,mid_ac});
				lq.push({mid_ac,c});
				lq.push({c,mid_cb});
				lq.push({mid_cb,r.b});
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

double integralPar(double left, double right, int nt=1)
{
	Opened opened;

	opened.push({left,right});
	int iters = 10;
	while (iters--)
	{
		range r = opened.top();
		opened.pop();

		double c = (r.a + r.b)*0.5;
		double mid_ac = (r.a + c)*0.5;
		double mid_cb = (c + r.b)*0.5;

		opened.push({r.a,mid_ac});
		opened.push({mid_ac,c});
		opened.push({c,mid_cb});
		opened.push({mid_cb,r.b});
	}
	//std::cout << "initial size: " << opened.size() << '\n';
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

	int times = 10;
	std::vector<double> r(times);
	auto t1 = high_resolution_clock::now();

	for (int i = 0; i < times; i++)
		r[i] = integralPar(0, 1, 8);

	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> dur(t2 - t1);

	std::cout << std::setprecision(16) << std::fixed << r.front() << '\n';
	std::cout << dur.count() / times << '\n';
	std::cout << n_iters << '\n';

	for (auto& p : count)
		std::cout << p.first << " / " << p.second << '\n';
	std::cout << count.size() << '\n';

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