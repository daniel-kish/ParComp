#include <iostream>
#include <memory>
#include <unordered_map>

std::unordered_map<void*, int> map;

template <class T>
struct sp
{
	T* p;
	explicit sp(T* ptr) : p{ptr} {
		map[p]++;
	}
	~sp() {
		map[p]--;
		if (map[p] == 0) {
			delete p;
			auto it = map.find(p);
			map.erase(it);
		}
	}
};

struct S {
	sp<S> i;
	S() : i(nullptr) {

	}
	~S()
	{
		std::cout << "dest " << this << '\n';
	}
};

int main()
{
	sp<S> s1(new S);
	sp<S> s2(new S);
	sp<S> s3(s2.p);
	std::cout << sizeof(s1) << '\n';
}