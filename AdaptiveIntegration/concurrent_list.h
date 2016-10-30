#pragma once
#include <deque>
#include <mutex>

template <class T, class Cont = std::deque<T>>
class concurrent_list
{
	using lg = std::lock_guard<std::mutex>;
	using iter = typename std::list<T>::iterator;
	mutable std::mutex m;
	Cont data;
public:
	concurrent_list() {}
	concurrent_list(concurrent_list const& other)
	{
		lg lk(other.m);
		data = other.data;
	}
	concurrent_list& operator= (concurrent_list const&) = delete;

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
	void push_front(T const& v)
	{
		lg lk(m);
		data.push_front(std::cref(v));
	}
	void push_front(T&& v)
	{
		lg lk(m);
		data.push_front(std::move(v));
	}

	bool try_pop_front(T& value)
	{
		lg lk(m);

		if (data.empty())
			return false;
		value = std::move(data.front());
		data.pop_front();
		return true;
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
	bool try_pop_back(T& value, int critical_size)
	{
		lg lk(m);

		if (data.size() < critical_size)
			return false;
		value = std::move(data.back());
		data.pop_back();
		return true;
	}
	iter begin() {
		return data.begin();
	}
	iter end() {
		return data.end();
	}
	bool empty() const {
		lg lk(m);
		return data.empty();
	}
	bool size() const {
		lg lk(m);
		return data.size();
	}
};