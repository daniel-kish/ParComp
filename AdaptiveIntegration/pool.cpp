#include "pool.h"
#include <iostream>

void pool::worker()
{
	while (!done)
	{
		range r;
		if (tasks.try_pop_back(r)) {
			this->est_and_branch(r);
		}
		else
			std::this_thread::yield();
	}
}

void pool::est_and_branch(range r)
{
	double c = (r.a + r.b)*0.5;
	double fa = f(r.a), fb = f(r.b);
	double fc = f(c);
	double s1 = trapeze(fa, fb, r.b - r.a);
	double s2 = trapeze(fa, fc, c - r.a) + trapeze(fc, fb, r.b - c);
	double change = std::abs(s2 - s1);

	if (change < eps)
		this->sum += s2;
	else {
		double c = (r.a + r.b)*0.5;
		double mid_ac = (r.a + c)*0.5;
		double mid_cb = (c + r.b)*0.5;

		tasks.push_front({r.a,mid_ac});
		tasks.push_front({mid_ac,c});
		tasks.push_front({c,mid_cb});
		tasks.push_front({mid_cb,r.b});
	}
}

pool::pool(int n_threads) : done(false)
{
	if (n_threads > std::thread::hardware_concurrency())
		std::cerr << "warning: more threads have been acquired than there is cores\n";
	threads.reserve(n_threads);
	for (int i = 0; i < n_threads; i++)
		threads.emplace_back(&pool::worker,this);
}

pool::~pool()
{
	done = true;
	//for (auto& t : threads)
	//	t.join();
}

void pool::submit(range r)
{
	tasks.push_front(r);
}

double pool::get_result()
{
	for (auto& t : threads)
		t.join();
	return sum;
}