namespace YUI
{
template<typename T, typename L, typename R, typename OP>	
class VectorExpr;

template<typename T, int DIM>
class Vector{
private:

	T value[DIM];

public:

	Vector()
	{
	}
	
	Vector(const Vector<T, DIM> &rhs)
	{
		for (int i = 0; i < DIM; ++i) {
			this->value[i] = rhs[i];
		}
	}
	
	Vector(const T *value)
	{
		for (int i = 0; i < DIM; ++i) {
			this->value[i] = value[i];
		}
	}
	
	template<typename L, typename R, typename OP>
	Vector(const VectorExpr<T, L, R, OP> &rhs);
	
	
	virtual ~Vector()
	{
	}
	
	
	
	T &operator[](int i)
	{
		return value[i];
	}
	
	const T &operator[](int i) const
	{
		return value[i];
	}
	
	Vector<T, DIM> &operator=(const Vector<T, DIM> &rhs)
	{
		if (this != &rhs) {
			for (int i = 0; i < DIM; ++i) {
				(*this)[i] = rhs[i];
			}
		}
		
		return *this;
	}
	
	template<typename R>
	Vector &operator=(const R &rhs)
	{
		for (int i = 0; i < DIM; ++i) {
			value[i] = rhs[i];
		}
		return *this;
	}
};

template<typename T, typename L, typename R, typename OP>
class VectorExpr{
	
	const L &lhs_;
	const R &rhs_;
	
public:

	VectorExpr(const L &lhs, const R &rhs)
	:lhs_(lhs), rhs_(rhs)
	{
	}
	
	inline T operator[](int i) const
	{
		return OP::eval(lhs_[i], rhs_[i]);
	}
};

template<typename T>
class Plus{
public:

	static inline T eval(const T &lhs, const T &rhs)
	{
		return lhs + rhs;
	}
};

template<typename T, int DIM, typename R>
VectorExpr<T, Vector<T, DIM>, R, Plus<T> > operator+(const Vector<T, DIM> &lhs, const R &rhs)
{
	return VectorExpr<T, Vector<T, DIM>, R, Plus<T> >(lhs, rhs);
}

template<typename T, int DIM, typename L>
VectorExpr<T, L, Vector<T, DIM>, Plus<T> > operator+(const L &lhs, const Vector<T, DIM> &rhs)
{
	return VectorExpr<T, L, Vector<T, DIM>, Plus<T> >(lhs, rhs);
}

template<typename T, int DIM>
VectorExpr<T, Vector<T, DIM>, Vector<T, DIM>, Plus<T> > operator+(const Vector<T, DIM> &lhs, const Vector<T, DIM> &rhs)
{
	return VectorExpr<T, Vector<T, DIM>, Vector<T, DIM>, Plus<T> >(lhs, rhs);
}

/* inner product  */
template<typename T, int DIM, typename R>
T operator*(const Vector<T, DIM> &lhs, const R &rhs)
{
	T result = lhs[0]*rhs[0];
	for (int i = 1; i < DIM; ++i) {
		result += lhs[i]*rhs[i];
	}
	
	return result;
}
template<typename T, int DIM, typename L>
T operator*(const L &lhs, const Vector<T, DIM> &rhs)
{
	T result = lhs[0]*rhs[0];
	for (int i = 1; i < DIM; ++i) {
		result += lhs[i]*rhs[i];
	}
	
	return result;
}

template<typename T, typename T2, int DIM>
T operator*(const Vector<T, DIM> &lhs, const Vector<T2, DIM> &rhs)
{
	T result = lhs[0]*rhs[0];
	for (int i = 1; i < DIM; ++i) {
		result += lhs[i]*rhs[i];
	}
	
	return result;
}


template<typename T, int DIM>
template<typename L, typename R, typename OP>
Vector<T, DIM>::Vector<T, DIM>(const VectorExpr<T, L, R, OP> &rhs)
{
	for (int i = 0; i < DIM; ++i) {
		value[i] = rhs[i];
	}
}


} /* namespace YUI */


#include <iostream>

using namespace std;
using YUI;

int main()
{
	cout << "som started" << endl;
	
	Vector<int, 10> vec1((int[]){0,1,2,3,4,5,6,7,8,9});
	cout << vec1[3] << endl;
	Vector<int , 10> vec2;
	vec2 = vec1+vec1;
	cout << (vec1 + vec2)[3] << endl;
	Vector<int, 10> vec3;
	vec3 = vec2 + (vec1 + vec2);
	cout << vec3[3] << endl;
	Vector<int, 15> vec4;
	cout << vec3*vec2 << endl;
	cout << vec3*(vec2 + vec1) << endl;
	
	return 0;
}

