#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <list>
#include <random>
#include <chrono>
#include <thread>
#include <future>

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

std::atomic<int> thread_count;
const int max_thread_count{4};

void merge_sort(list<int>& l)
{
	if (l.size() < 100) {
		l.sort();
		return;
	}
	auto mid = l.begin();
	std::advance(mid, l.size() / 2);

	list<int> p,q;
	p.splice(p.begin(), l, l.begin(), mid);
	q.splice(q.begin(), l, mid, l.end());

	std::thread t1, t2;
	if (thread_count < max_thread_count) {
		t1 = std::thread{merge_sort,std::ref(p)};
		thread_count++;
	}
	else
		merge_sort(p);

	if (thread_count < max_thread_count) {
		t2 = std::thread{merge_sort,std::ref(q)};
		thread_count++;
	}
	else
		merge_sort(q);

	if (t1.joinable()) {
		t1.join();
		thread_count--;
	}
	if (t2.joinable()) {
		t2.join();
		thread_count--;
	}

	my_merge(p, q, l);
}

// 443

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
	cout << "elapsed " << dur.count() << " ms" << '\n';

}