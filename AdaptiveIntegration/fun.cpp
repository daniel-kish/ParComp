#include <cmath>
#include <iostream>
#include "fun.h"

const double eps = 0.000000000001;
const double Pi = 4.0*atan(1.0);

double f(double x)
{
	return sq(sin(Pi*x)) + sq(sin(2.0*Pi*x)) + cos(200.0*x*x*x*x);
}

double trapeze(double fa, double fb, double step)
{
	return (fa + fb)*0.5*step;
}

std::ostream& operator<< (std::ostream& os, range const& r)
{
	os << '(' << r.a << ' ' << r.b << ')';
	return os;
}