#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <thread>

template <class Fun>
double t_integrate(double a, double b, Fun f, int n=0.6e7)
{
	double s = 0.0;
	const double h = (b - a) / n;

	for (int i = 1; i < n; ++i)
	{
		s += f(a + i*h);
	}

	s += (f(a) + f(b)) / 2.0;
	s *= h;
	return s;
}

double ptr_integrate(const double a, const double b, double(* const f)(double), const int n = 0.6e7)
{
	double s = 0.0;
	const double h = (b - a) / n;

	for (int i = 1; i < n; ++i)
	{
		s += f(a + i*h);
	}
	s += (f(a) + f(b)) / 2.0;
	s *= h;
	return s;
}

double fun(double x)
{
	return std::sin(x) + std::sin(2.0*x) + std::sin(4.0*x);
}

int main()
{	
	using namespace std;
	using namespace std::chrono;
	using namespace std::chrono_literals;

	double a = 0.0, b = 4.0*atan(1.0);
	int n = 1.0e7;

	auto f = [](double x) { return fun(x); };
	//volatile auto t1 = high_resolution_clock::now();
	volatile clock_t t1 = clock();

	double s = t_integrate(a, b, f, n);
	//double s = ptr_integrate(a, b, fun, n);
	
	volatile clock_t t2 = clock();
	//volatile auto t2 = high_resolution_clock::now();

	std::cout << s << " for " <<  1'000*(t2-t1)/CLOCKS_PER_SEC << "\n";
}