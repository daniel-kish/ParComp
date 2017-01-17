#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include <list>

using namespace std;
using namespace std::chrono;
int main()
{
	// 25 ms
	// 32 ms
	// 57 ms
	// 247 ms
	int times = 200;
	
	int sz = 50'000'000;
	vector<int> v(sz);
	list<int> l(sz);
	iota(begin(v), end(v), 100);
	iota(begin(l), end(l), 100);

	vector<int> ids(v.size());
	iota(begin(ids), end(ids), 0);

	auto t1 = high_resolution_clock::now();
	for (int i=0; i < times; i++) 
		for (auto& r : v)
			r *= 2;
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> d(t2 - t1);
	std::cout << d.count()/times << '\n';

	t1 = high_resolution_clock::now();
	for (int i = 0; i < times; i++)
		for (auto rp = v.rbegin(); rp != v.rend(); ++rp)
			*rp *= 2;
	t2 = high_resolution_clock::now();
	d = t2 - t1;
	std::cout << d.count()/times << '\n';

	t1 = high_resolution_clock::now();
	for (int i = 0; i < times; i++)
		for (auto id : ids)
			v[id] *= 2;
	t2 = high_resolution_clock::now();
	d = t2 - t1;
	std::cout << d.count()/times << '\n';

	t1 = high_resolution_clock::now();
	for (int i = 0; i < times; i++)
		for (auto &i : l)
			i *= 2;
	t2 = high_resolution_clock::now();
	d = t2 - t1;
	std::cout << d.count() / times << '\n';

}