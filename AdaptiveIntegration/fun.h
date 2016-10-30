#pragma once

inline double sq(double x) { return x*x; }

struct range { double a, b; };

std::ostream& operator<< (std::ostream& os, range const& r);