#include "Point.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

std::vector<double> xs;

Point::Point()
{}

Point::Point(double x_) : mx{x_}, fv{fun(mx)}
{}

double Point::fun(double x)
{
	//xs.push_back(x);
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
	os << '{' << a.x() << ' ' << a.f() << '}';
	return os;
}

double dist(Point const& p, Point const& q)
{
	return q.x() - p.x();
}

Point middle(Point const& a, Point const& b)
{
	Point p{(a.x() + b.x())*0.5};
	return p;
}

double trapezium_area(Point const& a, Point const& b)
{
	return (a.f() + b.f()) / 2.0 * dist(a, b);
}