#pragma once
#include <memory>

template<typename T>
class threadsafe_queue
{
private:
	using guard = std::lock_guard<std::mutex>;
	struct node
	{
		std::shared_ptr<T> data;
		std::unique_ptr<node> next;
	};
	std::mutex head_mutex;
	std::unique_ptr<node> head;
	std::mutex tail_mutex;
	node* tail;

	node* read_tail()
	{
		 guard tail_lock{tail_mutex};
		 return tail;
	}

	std::unique_ptr<node> pop_head()
	{
		guard head_lock{head_mutex};
		if (head.get() == read_tail())
			return nullptr;
		std::unique_ptr<node> old_head = std::move(head);
		head = std::move(old_head->next);
		return old_head;
	}

public:
	threadsafe_queue() :
		head(new node), tail(head.get())
	{}
	threadsafe_queue(const threadsafe_queue& other) = delete;
	threadsafe_queue& operator=(const threadsafe_queue& other) = delete;
	
	bool empty()
	{
		return head.get() == read_tail();
	}

	std::shared_ptr<T> try_pop()
	{
		std::unique_ptr<node> old_head = pop_head();
		if (old_head)
			return old_head->data;
		return std::shared_ptr<T>{};
	}
	void push(T new_value)
	{
		std::shared_ptr<T> new_data = std::make_shared<T>(std::move(new_value));
		std::unique_ptr<node> p(new node);
		node* const new_tail = p.get();
		// critical from here to the end
		std::lock_guard<std::mutex> tail_lock{tail_mutex};
		tail->data = new_data;
		tail->next = std::move(p);
		tail = new_tail;
	}
};