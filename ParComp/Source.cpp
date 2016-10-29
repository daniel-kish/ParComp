#include <cmath>

double linkSin(double x)
{
	return std::sin(x) + std::sin(2.0*x) + std::sin(4.0*x);
}