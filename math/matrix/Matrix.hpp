#ifndef __MATRIX_HPP
#define __MATRIX_HPP

#include <iostream>

namespace math {
	
	
template<class T = double, size_t rows=3, size_t cols=3>
class Matrix 
	:: public Tensor<T>
{ 
// Public data members
public: 
	

// Protected data members
protected: 
	// the physical data
	std::vector<T> data[rows*cols];

// Public member functions
public: 

	size_t numel() {
		return rows*cols;
	}
	
	size_t rows() {
		return rows;
	}
	size_t columns() { 
		return cols; 
	}
	
	

	// Operator overloads
	
	template<class T2, class T3>
	Matrix<T3, rows, cols2>
	operator* (Matrix<T2,cols,rows2> &other) {
	}
	
	Matrix
	operator+ (Matrix &other) {
		
	}
	
	Matrix
	operator- (Matrix &other) {
	}

	Matrix
	operator+ () { 
		
	}
	
	Matrix
	operator- () {
		
	}

	

// Protected member functions
protected:

};


// Some matrix generators
template<class T, size_t rows, size_t cols>
Matrix<T,rows,cols> ones(size_t rows, size_t cols) {
	
}

template<class T, size_t rows, size_t cols>
Matrix<T,rows,cols> zeros(size_t rows, size_t cols) {
	
}

template<class T, size_t rows, size_t cols>
Matrix<T,rows,cols> NaN(size_t rows, size_t cols) {
	
}


template<class T, size_t rows, size_t cols>
Matrix<T,rows,cols> Inf(size_t rows, size_t cols) {
	
}


// Some handy aliases
typedef math::Matrix<T, 3ul,1ul>    Vector<T>;
typedef math::Matrix<T, 3ul,1ul>    RowVector<T>;
typedef math::Matrix<T, 1ul,3ul>    ColumnVector<T>;



// Vector-specific operations

template <class T1, class T2>
Vector<T2> 
cross(Vector<T1> &other) { 
	
}

template <class T1, class T2>
Vector <T2>
dot (Vector<T1> &other) {
	
}

norm()


// MatrixMask  (logical indexing)
template<class T>
class MatrixMask 
{
	public: 
	protected:
	
	
};

} // namespace math





#endif
