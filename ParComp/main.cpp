#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>

template <class Fun>
double integrate(const double a, const double b, Fun f, const int n=0.6e7)
{
	//std::cout << typeid(f).name() << '\n';
	double s = 0.0;
	const double h = (b - a) / n;
#pragma omp parallel for reduction(+:s) 
	for (int i = 1; i < n; ++i)
	{
		s += f(a + i*h);
	}
	s += (f(a) + f(b)) / 2.0;
	s *= h;
	return s;
}

int main()
{
	double a = 0.0, b = 4.0*atan(1.0);
	double omega = 3.0;
	auto f = [omega](double x) { return std::sin(omega*x) + std::sin(2.0*x) + std::sin(1.0*x); };

	using std::chrono::high_resolution_clock;
	using std::chrono::duration;

	omp_set_num_threads(4);
	
	double s = 0.0;
	auto t1 = high_resolution_clock::now();
	
	s = integrate(a, b, f, 1.0e8);

	auto t2 = high_resolution_clock::now();
	duration <double, std::milli> dur(t2 - t1);
	std::cout << std::setprecision(12) << std::fixed << s << '\n';
	std::cout << std::setprecision(2) << std::fixed;
	std::cout << dur.count() << '\n';
}
