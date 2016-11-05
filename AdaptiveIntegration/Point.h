#pragma once
#include <iostream>

class Point
{
	double mx;
	double fv;
	double fun(double x);
public:
	Point();
	Point(double x_);
	double x() const;
	double f() const;
	void set(double newx);
};

std::ostream& operator<< (std::ostream& os, Point const& a);
double dist(Point const& p, Point const& q);
double trapeze(Point const& a, Point const& b);