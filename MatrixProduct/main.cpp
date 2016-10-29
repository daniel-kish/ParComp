#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <mutex>

std::vector<int> access;
std::mutex mx;

template <class T>
void MatMul(std::vector<T>& vC, const std::vector<T>& vA,
	const std::vector<T>& vB, int height, int width, int length)
{
	// vA[height X width]
	// vB[width X length]
	// vC[height X length]

#pragma omp parallel for
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < length; col++)
		{
			T sum{};
			for (int i = 0; i < width; i++)
				sum += vA[row * width + i] * vB[i * length + col];
			vC[row * length + col] = sum;
#ifdef DEBUG
			{
				std::lock_guard<std::mutex> lk(mx);
				access[row*length + col] = omp_get_thread_num();
			}
#endif
		}
	}
}


//template <class T>
//void MatMulBlock(std::vector<T>& vC, const std::vector<T>& vA,
//	const std::vector<T>& vB, int height, int width, int length)
//{
//	// vA[height X width]
//	// vB[width X length]
//	// vC[height X length]
//
//	int n_blocks = 1;
//	int per_block = length / n_blocks;
//	for (int t = 0; t < nblocks - 1; ++t)
//	{
//		// do multiplication [t*per_block : (t+1)*per_block)
//	}
//	assert(t == n_blocks - 1);
//	// do multiplication [t*per_block : length)
//}

template <class T>
void print(std::vector<T>& v, int M, int N)
{
	for (int row = 0; row < M; row++)
	{
		for (int col = 0; col < N; col++)
		{
			std::cout << v[row*N + col] << ' ';
		}
		std::cout << '\n';
	}
}

int main()
{	
	using namespace std::chrono;

	int n = 1;
	int height = 1024/n;
	int width  = 1024;
	int length = 1024*n;

#ifdef DEBUG
	access.resize(height*length);
#endif
	std::vector<double> a(height*width,2.5), b(width*length,2), c(height*length);
	
	omp_set_num_threads(4);
	auto t1 = high_resolution_clock::now();
	
	MatMul(c, a, b, height, width, length);

	auto t2 = high_resolution_clock::now();

	std::cout << duration_cast<milliseconds>(t2 - t1).count() << '\n';
	
#ifdef DEBUG
	//print(c, height, length);
	print(access, height, length);
#endif
}