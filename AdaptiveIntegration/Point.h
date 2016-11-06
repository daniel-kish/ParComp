#pragma once
#include <iostream>
#include <vector>

class Point
{
	double mx;
	double fv;
	static double fun(double x);
public:
	Point();
	Point(double x_);
	double x() const;
	double f() const;
	void set(double newx);
};

std::ostream& operator<< (std::ostream& os, Point const& a);
double dist(Point const& p, Point const& q);
Point middle(Point const& a, Point const& b);
double trapezium_area(Point const& a, Point const& b);

extern std::vector<double> xs;