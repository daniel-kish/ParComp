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

int threadsNum = 0;

template <typename T> void swap(T *el1, T *el2)
{
	T temp = *el1;
	*el1 = *el2;
	*el2 = temp;
	return;
}
#define THREADS_NUM_MAX 0
mutex mio;
using guard = lock_guard<mutex>;

template <typename T>
void mergeSort(T *mass, int left, int right)
{
	{guard l(mio);  cout << this_thread::get_id() << ": " << left << ' ' << right << '\n'; }
	if (right - left <= 0)
	{
		return;
	}
	else if (right - left == 1)
	{
		if (mass[left] > mass[right])
		{
			swap(mass + left, mass + right);
		}
		return;
	}
	else
	{
		int i, j, k, newR = (right - left) / 2, newL = 1 + ((right - left) / 2);
		int n1 = newR - left + 1, n2 = right - newL + 1;
		
		if (threadsNum < THREADS_NUM_MAX)
		{
			std::future<void> merge = std::async(std::launch::async, mergeSort<T>, mass, newL, right);
			++threadsNum;
			mergeSort(mass, left, newR);
			merge.wait();
			--threadsNum;
		}
		else
		{
			mergeSort(mass, left, newR);
			mergeSort(mass, newL, right);
		}
		T *halfMass1 = new T[n1];
		T *halfMass2 = new T[n2];
		for (i = 0; i < n1; ++i)
		{
			halfMass1[i] = mass[left + i];
		}
		for (i = 0; i < n2; ++i)
		{
			halfMass2[i] = mass[newL + i];
		}
		i = 0; j = 0; k = left;
		do
		{
			if (i < n1 && j < n2)
			{
				if (halfMass1[i] < halfMass2[j])
				{
					mass[k] = halfMass1[i];
					++i;
				}
				else
				{
					mass[k] = halfMass2[j];
					++j;
				}
			}
			else
			{
				if (i == n1)
				{
					mass[k] = halfMass2[j];
					++j;
				}
				else
				{
					mass[k] = halfMass1[i];
					++i;
				}
			}
			++k;
		} while (k < n1 + n2);
		delete[](halfMass1);
		delete[](halfMass2);
		return;
	}
}

int main()
{
	//list<int> l;
	//
	//random_device rd;
	//mt19937 mt{rd()};
	//uniform_int_distribution<int> d(-100, 100);
	//for (int i = 0; i < 1'000'000; ++i)
	//	l.push_back(d(mt));
	//
	//auto t1 = high_resolution_clock::now();
	//merge_sort(l);
	//auto t2 = high_resolution_clock::now();
	//duration<double, std::milli> dur(t2-t1);
	//
	//cout << boolalpha << is_sorted(begin(l), end(l)) << '\n';
	//cout << "elapsed " << dur.count() << '\n';

	vector<int> v(20);
	iota(begin(v), end(v), 1);
	shuffle(begin(v), end(v), mt19937{});
	cout << v << '\n';
	mergeSort(v.data(), 0, v.size() - 1);
	cout << v << '\n';
}