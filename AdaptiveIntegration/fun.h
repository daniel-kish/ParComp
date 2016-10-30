#pragma once
#include "concurrent_list.h"

inline double sq(double x) { return x*x; }

struct range { double a, b; };

std::ostream& operator<< (std::ostream& os, range const& r);

double trapeze(double fa, double fb, double step);
double f(double x);

extern const double eps;
extern const double Pi;

