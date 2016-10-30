#pragma once
#include <atomic>
#include "concurrent_list.h"
#include "fun.h"
#include <vector>
#include <thread>
#include <cassert>

class pool
{
	std::atomic_bool done;
	concurrent_list<range> tasks;
	std::vector<std::thread> threads;
	double sum{0.0};
	
	void worker();
	void est_and_branch(range r);
public:
	pool(int n_threads);
	~pool();

	void submit(range r);
	double get_result();
};