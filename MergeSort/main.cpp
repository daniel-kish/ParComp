#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <list>
#include <random>
#include <chrono>
#include <thread>

using namespace std;
using namespace std::chrono;

template <class T>
ostream& operator<< (ostream& os, list<T> const& v)
{
	for (auto& e : v)
		os << e << ' ';
	return os << '\n';
}
template <class T>
ostream& operator<< (ostream& os, vector<T> const& v)
{
	for (auto& e : v)
		os << e << ' ';
	return os << '\n';
}

void my_merge(list<int>& p, list<int>& q, list<int>& l)
{
	while (true)
	{
		if (p.empty() && !q.empty()) {
			l.splice(l.end(), q);
			break;
		}
		else if (q.empty() && !p.empty()) {
			l.splice(l.end(), p);
			break;
		}
		else if (q.empty() && p.empty()) {
			break;
		}

		auto l_it = (l.size() > 0)? l.end() : l.begin();
		if (p.front() <= q.front())
			l.splice(l_it, p, p.begin());
		else
			l.splice(l_it, q, q.begin());
	}
}

void merge_sort(list<int>& l)
{
	auto mid = l.begin();
	std::advance(mid, l.size() / 2);

	list<int> p,q;
	p.splice(p.begin(), l, l.begin(), mid);
	q.splice(q.begin(), l, mid, l.end());

	p.sort();
	q.sort();
	//auto sorter = [] (list<int>& l){
	//	auto mid = l.begin();
	//	std::advance(mid, l.size() / 2);
	//	list<int> p, q;
	//	p.splice(p.begin(), l, l.begin(), mid);
	//	q.splice(q.begin(), l, mid, l.end());
	//	auto sorter = [](list<int>& l) {
	//		l.sort();
	//	};
	//	thread t1{ sorter, p }, t2{ sorter,q };
	//	t1.join(); t2.join();
	//	my_merge(p, q, l);
	//};
	//thread t1{ sorter, p }, t2{sorter,q};
	//t1.join(); t2.join();

	my_merge(p, q, l);
}

int main()
{
	list<int> l;
	
	random_device rd;
	mt19937 mt{rd()};
	uniform_int_distribution<int> d(-100, 100);
	for (int i = 0; i < 1'000'000; ++i)
		l.push_back(d(mt));

	auto t1 = high_resolution_clock::now();
	merge_sort(l);
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> dur(t2-t1);

	cout << boolalpha << is_sorted(begin(l), end(l)) << '\n';
	cout << "elapsed " << dur.count() << '\n';
}