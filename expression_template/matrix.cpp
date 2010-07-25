#include <iostream>
using namespace std;

namespace YUI
{

template<typename T, int ROW, int COL, typename L, typename R, typename OP>
class MatrixExpr;

template<typename T, int ROW, int COL>
class Matrix{
public:

	typedef T ValueType;
	static const int ROW_ = ROW;
	static const int COL_ = COL;

private:

	T value[ROW][COL];

public:

	Matrix()
	{
	}
	
	Matrix(const T *value)
	{
		for (int i = 0; i < ROW; ++i) {
			for (int j = 0; j < COL; ++j) {
				this->value[i][j] = value[i*COL + j];
			}
		}
	}
	
	Matrix(const T **value)
	{
		for (int i = 0; i < ROW; ++i) {
			for (int j = 0; j < COL; ++j) {
				this->value[i][j] = value[i][j];
			}
		}
	}
	
	T &operator()(int i, int j)
	{
		return value[i][j];
	}
	
	const T &operator()(int i, int j) const
	{
		return value[i][j];
	}
	
	template<typename R>
	Matrix<T, ROW, COL> &operator=(const R &rhs)
	{
		if ((void*)this == (void*)&rhs ) {
			return *this;
		}
		
		for (int i = 0; i < ROW; ++i) {
			for (int j = 0; j < COL; ++j) {
				(*this)(i, j) = rhs(i, j);
			}
		}
	
		return *this;
	}
};

template<typename L, typename R, typename T>
class Plus{
public:

	static inline T eval(const L &lhs, const R &rhs, int i, int j)
	{
		return lhs(i, j) + rhs(i, j);
	}
};

template<typename L, typename R, typename T>
class Prod{
public:

	static inline T eval(const L &lhs, const R &rhs, int i, int j)
	{
		T result = lhs(i, 0) * rhs(0, j);
		
		for (int k = 1; k < L::COL_; ++k) {
			result += lhs(i, k) * rhs(k, j);
		} 
		
		return result;
	}
};


template<typename T, int ROW, int COL, typename L, typename R, typename OP>
class MatrixExpr{
	const L &lhs_;
	const R &rhs_;
	
public:

	MatrixExpr(const L &lhs, const R &rhs)
	: lhs_(lhs), rhs_(rhs)
	{
	}
	
	inline T operator()(int i, int j) const
	{
		return OP::eval(lhs_, rhs_, i, j);
	}
};


template<typename T, int ROW, int COL, typename T2, typename L, typename R, typename OP>
MatrixExpr<T, ROW, COL, Matrix<T, ROW, COL>, MatrixExpr<T2, ROW, COL, L, R, OP>, Plus<Matrix<T, ROW, COL>, MatrixExpr<T2, ROW, COL, L, R, OP>, T> >
operator+(const Matrix<T, ROW, COL> &lhs, const MatrixExpr<T2, ROW, COL, L, R, OP> &rhs)
{
	return MatrixExpr<T, ROW, COL, Matrix<T, ROW, COL>, MatrixExpr<T2, ROW, COL, L, R, OP>, Plus<Matrix<T, ROW, COL>, MatrixExpr<T2, ROW, COL, L, R, OP>, T> >(lhs, rhs);
}

template<typename T, int ROW, int COL, typename L, typename R, typename T2, typename OP>
MatrixExpr<T, ROW, COL, MatrixExpr<T, ROW, COL, L, R, OP>, Matrix<T2, ROW, COL>, Plus<MatrixExpr<T, ROW, COL, L, R, OP>, Matrix<T2, ROW, COL>, T> >
operator+(const MatrixExpr<T, ROW, COL, L, R, OP> &lhs, const Matrix<T2, ROW, COL> &rhs)
{
	return MatrixExpr<T, ROW, COL, MatrixExpr<T, ROW, COL, L, R, OP>, Matrix<T2, ROW, COL>, Plus<MatrixExpr<T, ROW, COL, L, R, OP>, Matrix<T2, ROW, COL>, T> >(lhs, rhs);
}

template<typename T, int ROW, int COL, typename T2, typename L1, typename R1, typename L2, typename  R2, typename OP1, typename OP2>
MatrixExpr<T, ROW, COL, MatrixExpr<T, ROW, COL, L1, R1, OP1>, MatrixExpr<T2, ROW, COL, L2, R2, OP2>, Plus<MatrixExpr<T, ROW, COL, L1, R1, OP1>, MatrixExpr<T2, ROW, COL, L2, R2, OP2>, T> >
operator+(const MatrixExpr<T, ROW, COL, L1, R1, OP1> &lhs, const MatrixExpr<T2, ROW, COL, L2, R2, OP2> &rhs)
{
	return MatrixExpr<T, ROW, COL, MatrixExpr<T, ROW, COL, L1, R1, OP1>, MatrixExpr<T2, ROW, COL, L2, R2, OP2>, Plus<MatrixExpr<T, ROW, COL, L1, R1, OP1>, MatrixExpr<T2, ROW, COL, L2, R2, OP2>, T> >(lhs, rhs);
}

template<typename T, int ROW, int COL, typename T2>
MatrixExpr<T, ROW, COL, Matrix<T, ROW, COL>, Matrix<T2, ROW, COL>, Plus<Matrix<T, ROW, COL>, Matrix<T2, ROW, COL>, T> >
operator+(const Matrix<T, ROW, COL> &lhs, const Matrix<T2, ROW, COL> &rhs)
{
	return MatrixExpr<T, ROW, COL, Matrix<T, ROW, COL>, Matrix<T2, ROW, COL>, Plus<Matrix<T, ROW, COL>, Matrix<T2, ROW, COL>, T> >(lhs, rhs);
}


template<typename T, int ROW, int COL, int COL2, typename T2>
MatrixExpr<T, ROW, COL2, Matrix<T, ROW, COL>, Matrix<T2, COL, COL2>, Prod<Matrix<T, ROW, COL>, Matrix<T2, COL, COL2>, T> >
operator*(const Matrix<T, ROW, COL> &lhs, const Matrix<T2, COL, COL2> &rhs)
{
	return MatrixExpr<T, ROW, COL2, Matrix<T, ROW, COL>, Matrix<T2, COL, COL2>, Prod<Matrix<T, ROW, COL>, Matrix<T2, COL, COL2>, T> >(lhs, rhs);
}

}

#include <iostream>
using namespace std;
using namespace YUI;

int main()
{
	Matrix<int, 10, 10> mat;
	
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			mat(i,j) = i*j;
		}
	}
	
	Matrix<int, 10, 10> mat2;
	mat2 = mat + mat;
	
	cout << mat2(1,2) << endl;
	
	cout << (mat + mat2 + (mat + mat))(1,2) << endl;
	cout << (mat*mat2)(1,2) << endl;

	
	return 0;
}

