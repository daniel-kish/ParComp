#pragma once
#include <deque>
#include <mutex>
#include <list>
#include <vector>

template <class T>
class thread_safe_stack
{
	using lg = std::lock_guard<std::mutex>;
	using iter = typename std::deque<T>::iterator;
	mutable std::mutex m;
	std::condition_variable cv;
	std::deque<T> data;
public:
	thread_safe_stack() {}
	thread_safe_stack(thread_safe_stack const& other)
	{
		lg lk(other.m);
		data = other.data;
	}
	thread_safe_stack& operator= (thread_safe_stack const&) = delete;

	void push_back(T const& v)
	{
		lg lk(m);
		data.push_back(std::cref(v));
	}
	void push_back(T&& v)
	{
		lg lk(m);
		data.push_back(std::move(v));
	}

	bool try_pop_back(T& value)
	{
		lg lk(m);

		if (data.empty())
			return false;
		value = std::move(data.back());
		data.pop_back();
		return true;
	}
	bool empty() const {
		lg lk(m);
		return data.empty();
	}
	bool size() const {
		lg lk(m);
		return data.size();
	}
	void wait_until_empty()
	{
		std::unique_lock<std::mutex> lk(m);
		cv.wait(lk, [this] {return data.empty(); });
	}
};