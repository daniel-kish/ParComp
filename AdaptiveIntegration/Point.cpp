#include "Point.h"
#define _USE_MATH_DEFINES
#include <math.h>

Point::Point()
{}

Point::Point(double x_) : mx{x_}, fv{fun(mx)}
{}

double Point::fun(double x)
{
	return sin(x);
}

double Point::x() const { return mx; }

double Point::f() const { return fv; }

void Point::set(double newx)
{
	mx = newx;
	fv = fun(mx);
}

std::ostream& operator<< (std::ostream& os, Point const& a)
{
	os << '{' << a.x() << '_' << a.f() << '}';
	return os;
}

double dist(Point const& p, Point const& q)
{
	return q.x() - p.x();
}


double trapeze(Point const& a, Point const& b)
{
	return (a.f() + b.f())*dist(a, b) / 2.0;
}