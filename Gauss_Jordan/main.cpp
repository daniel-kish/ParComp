#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <chrono>
#include <cassert>
#include <omp.h>

template <class T>
void MatMul(std::vector<T>& vC, const std::vector<T>& vA,
	const std::vector<T>& vB, int s)
{
	int cols = s + 1;
	for (int row = 0; row < s; row++)
	{
		T sum{};
		for (int col = 0; col < s; col++)
			sum += vA[row * cols + col] * vB[col];
		vC[row] = sum;
	}
}

void print(std::vector<double> const& mat, int h, int w=-1)
{
	if (w == -1) w = h + 1;

	std::cout << std::setprecision(3);
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
			std::cout << std::setw(6) << mat[i*w + j] << ' ';
		std::cout << '\n';
	}
	std::cout << '\n';
}

void gauss_fwd(/*std::vector<double>&*/double* m, int s)
{
	int cols = s + 1;
	for (int key_row = 0; key_row < s; key_row++)
	{
		double pivot = m[key_row*cols + key_row];
		for (int c = key_row; c < cols; c++)
			m[key_row*cols + c] /= pivot;
#pragma omp parallel for
		for (int row = key_row + 1; row < s; row++)
		{
			double k = -m[row*cols + key_row] / m[key_row*cols + key_row];
			for (int c = key_row; c < cols; c++)
				m[row*cols + c] += m[key_row*cols + c] * k;
		}
	}
}

void gauss_bwd(/*std::vector<double>&*/ double* m, int s)
{
	int cols = s + 1;
	for (int key_row = s-1; key_row >= 0; key_row--)
	{
#pragma omp parallel for
		for (int r = key_row - 1; r >= 0; r--)
		{
			double k = -m[r*cols + key_row];
			for (int c = key_row; c < s + 1; c++)
			{
				m[r*cols + c] += m[key_row*cols + c] * k;
			}
		}
	}
}

int main()
{
	using namespace std::chrono;

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 2.0);
	auto gen = [&] { return dist(mt); };

	int s = 2000;

	std::vector<double> mat(s*(s + 1));
	std::generate(mat.begin(), mat.end(), gen);
	std::vector<double> sol(s);
	std::generate(sol.begin(), sol.end(), gen);
	std::vector<double> rhs(sol.size());
	MatMul(rhs, mat, sol, s); // make rhs
	for (int row = 0; row < s; row++)
		mat[row*(s + 1) + s] = rhs[row]; // copy to mx
	
	std::cout << "gen done\n";
	omp_set_num_threads(4);
	auto t1 = high_resolution_clock::now();
	gauss_fwd(mat.data(), s);
	gauss_bwd(mat.data(), s);
	auto t2 = high_resolution_clock::now();
	std::cout << duration_cast<milliseconds>(t2 - t1).count() << '\n';


	double res = 0.0;
	for (int i = 0; i < s; i++)
		res += std::abs(sol[i] - mat[i*(s + 1) + s]);
	std::cout << res << '\n';
}