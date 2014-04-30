#ifndef _MATRIX_HPP
#define _MATRIX_HPP


#include <stdexcept>
#include <valarray>
#include <cmath>
#include <limits.h>
#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <utility>

#include "PhysUnits.hh"


// Forward declarations
// --------------------------------------------------------------------------
class Slice;

template<typename T> class Vector;
template<typename T> class Matrix;

template<> class Vector<bool>;
template<> class Matrix<bool>;


// Handy defines, typedefs
// --------------------------------------------------------------------------
#define   START     0
#define   END       INT_MAX
#define   ALL       Matrix<>::Slice()


// class definitions
#if 1


// Matrix template class
// --------------------------------------------------------------------------
#if 1

template<typename T = double>
class Matrix
{
protected:
    class Cref;

public:

    // We have a friends!
    friend class Matrix<bool>;
    friend class Cref;


    // Slice class
    // --------------------------------------------------------------------------
    class Slice
    {
    public:

        // data member binds
        const unsigned int
            &start, &end;
        const bool
            &all;
        const int
            &step;

        // constructors
        Slice();
        Slice(unsigned int start, unsigned int end);
        Slice(unsigned int start, int step, unsigned int end);

        // return number of elements implied by the slice
        unsigned int
            length();

    private:

        // start, step, end
        unsigned int
             _start, _end;  // start and end must be positive
        int  _step;         // step may also be negative
        bool _all;          // useful to slice "all" columns or rows

    };

    // getters
    const unsigned int
        &rows, &cols, &numel;
    const bool
        &isempty, &isvector, &isscalar;
    const char*
        &datatype;

    // Matrix size
    Matrix<unsigned int>
        size() const;

    // constructor/destructor
    Matrix<T>();                                // empty matrix

    Matrix<T>(unsigned int r, unsigned int c);  // r-by-c matrix of zeros
    Matrix<T>(unsigned int i);                  // i-by-i square matrix of zeros
    Matrix<T>(unsigned int r, unsigned int c,
              T contents[]);                    // r-by-c matrix, with data in contents[]
    Matrix<T>(unsigned int r, unsigned int c,
              std::vector<T> contents);         // r-by-c matrix, with data in std::vector<T>

    Matrix<T>(unsigned int i, T contents[]);    // i-by-i square matrix, with data in contents[]
    Matrix<T>(unsigned int i,
              std::vector<T> contents);         // i-by-i square matrix, with data in std::vector<T>


    Matrix<T>(std::vector<T> contents);         // contents.size()-by-1 matrix, with data in std::vector<T>

    template<typename D>
    Matrix<T>(const Matrix<D> &M2);             // copy matrix of other type


    // type casting
    operator Matrix<float>    () { return Matrix<float>   (*this); }
    operator Matrix<double>   () { return Matrix<double>  (*this); }
    operator Matrix<char>     () { return Matrix<char>    (*this); }
    operator Matrix<int>      () { return Matrix<int>     (*this); }
    operator Matrix<uint8_t>  () { return Matrix<uint8_t> (*this); }
    operator Matrix<uint16_t> () { return Matrix<uint16_t>(*this); }
    operator Matrix<uint32_t> () { return Matrix<uint32_t>(*this); }
    operator Matrix<uint64_t> () { return Matrix<uint64_t>(*this); }
    operator Vector<T>        () { return Vector<T>(*this);        }

    operator T () {
        if (_isscalar)
            return M[0];
        throw std::runtime_error("Conversion to <T> only allowed for scalars.");
    }

    // compound assignments
    Matrix<T>
        & operator++(),    // prefix
        & operator++(int), // suffix
        & operator--(),    // prefix
        & operator--(int), // suffix
        & operator+=(const Matrix<T> &M2),  & operator+=(T a),
        & operator-=(const Matrix<T> &M2),  & operator-=(T a),
        & operator/=(const Matrix<T> &M2),  & operator/=(T a),
        & operator*=(const Matrix<T> &M2),  & operator*=(T a),
        & operator^=(const Matrix<T> &M2),  & operator^=(T a);

    // matrix multiplication, addition, etc.
    Matrix<T>
        operator*(T a) const,    operator*(const Matrix<T> &M2) const,
        operator/(T a) const,    operator/(const Matrix<T> &M2) const,
        operator^(T a) const,    operator^(const Matrix<T> &M2) const,
        operator+() const,       operator+(T a) const,
        operator-() const,       operator-(T a) const,
        operator+(const Matrix<T> &M2) const,
        operator-(const Matrix<T> &M2) const;

    // Matrix comparisons
    Matrix<bool>
        operator==(const Matrix<T> &M2) const,   operator==(T a) const,
        operator!=(const Matrix<T> &M2) const,   operator!=(T a) const,
        operator< (const Matrix<T> &M2) const,   operator< (T a) const,
        operator<=(const Matrix<T> &M2) const,   operator<=(T a) const,
        operator> (const Matrix<T> &M2) const,   operator> (T a) const,
        operator>=(const Matrix<T> &M2) const,   operator>=(T a) const;

    // assignment
    Matrix<T>&
        operator=(const Matrix<T> &M2);

    // transpose
    Matrix<T>
        transpose(),
        operator~();

    // indexing
    const T
        & operator()(unsigned int r, unsigned int c), // subscripts
        & operator()(unsigned int ind),               // linear index
        & operator[](unsigned int r, unsigned int c), // subscripts, NO bounds checks
        & operator[](unsigned int ind);               // linear index, NO bounds checks

         Cref operator()(Slice i, unsigned int j);       // slicing rows, assign data
    Matrix<T> operator()(Slice i, unsigned int j) const; // slicing rows, extract data

    Matrix<T>
        operator()(unsigned int i, Slice j) const, // slicing cols
        operator()(Slice i, Slice j) const,        // slicing both
        operator()(Slice i) const,                 // vectorizing
        operator()(const Matrix<bool> &M2) const,  // logical indexing
        operator()(const Matrix<unsigned int> &R, const Matrix<unsigned int> &C) const, // indexing with integer matrix
        operator()(const Matrix<unsigned int> &I) const; // logical indexing with integer matrix

    // first, last elements
    T& start() const { return M[0]; }
        T& first() const { return M[0]; }
        T& begin() const { return M[0]; }
    T& end()   const { return M[numel-1]; }
        T& final() const { return M[numel-1]; }
        T& last()  const { return M[numel-1]; }

    // reshape matrix
    Matrix<T>
        reshape(unsigned int i, unsigned int j);

    // print matrix to stdout
    void show();
        void print()  { show(); }
        void display(){ show(); }

    // Advanced operations
    T
        trace() const,
        determinant() const;

    Matrix<T>
        rot90() const,
        flipud() const,
        fliprl() const,
        circshift(Matrix<T> &dir),

        lowerTriangular() const,
        upperTriangular() const,

        diag() const,

        inv() const,
        invTimes() const,
        timesInv() const,

        LU() const,
        Cholesky() const,

        eig() const,
        eigs(unsigned int N) const,

        eigenVectors() const;


protected:

    // the Matrix' data is itself just a C-array
    T *M;
    // the Matrix' data is itself just a valarray
    std::valarray<T> *M;

    // matrix' size
    unsigned int
        _rows,
        _cols,
        _numel;

    // Matrix charachteristics
    bool
        _isempty,
        _isvector,
        _isscalar;

    // class of Matrix
    const char*
        _class;
    void
        _get_class(T&);


    // Cref class
    // Helper class to allow sliced assignments
    // --------------------------------------------------------------------------
    class Cref
    {
    public:

        operator       Matrix<T>  (){ return M; }
        operator const Matrix<T>& (){ return M; }
        operator const Matrix<T>  (){ return M; }

        const Matrix<T> & operator=(const Matrix<T> &M2);

    private:

        Matrix<T>& M;
        uint8_t mode;

        union rowind_t;  rowind_t rowind;
        union colind_t;  colind_t colind;

        // Mode 0: sliced row, constant column
        Cref(Matrix<T> &M, const Slice &R, unsigned int c);

        // Mode 1: constant row, sliced column
        Cref(Matrix<T> &M, unsigned int r, const Slice &C);

        // Mode 2: both indices sliced
        Cref(Matrix<T> &M, const Slice &R, const Slice &C);

        // Mode 3: linear slice
        Cref(Matrix<T> &M, const Slice &I);

        // Mode 4: boolmatrix
        Cref(Matrix<T> &M, const Matrix<bool> &I);

        // Mode 5: row with int matrix, constant col
        Cref(Matrix<T> &M, const Matrix<unsigned int> &iR, unsigned int c);

        // Mode 6: constant row, col with int matrix
        Cref(Matrix<T> &M, unsigned int r, const Matrix<unsigned int> &iC);

        // Mode 7: linear intmatrix
        Cref(Matrix<T> &M, const Matrix<unsigned int> &I);

    };

};

#endif



// Matrix class:: bool specialization
// --------------------------------------------------------------------------
#if 1

template<>
class Matrix<bool>
{
public:

    // Matrix<> is our friend!
    template<typename T> friend class Matrix;

    // getters
    unsigned int
        &rows, &cols, &numel;
    bool
        &isempty, &isvector;

    // constructor/destructor
    Matrix<bool>();                                // empty matrix
    Matrix<bool>(unsigned int r, unsigned int c);  // r-by-c matrix of false's
    Matrix<bool>(unsigned int rc);                 // r-by-r square matrix of false's
    Matrix<bool>(const Matrix<bool> &M2);          // copy matrix
   ~Matrix<bool>();                                // destroy matrix

    // Indexing, assignment
    Matrix<bool>&
        operator=(const Matrix<bool> &M2);
    bool& operator()(unsigned int r, unsigned int c) const;
    bool& operator()(unsigned int ind) const;
    Matrix<bool>
        operator()(const Matrix<bool> &M2) const;

    // Compound assignments (including logicals)
    Matrix<bool>
        operator&=(const Matrix<bool> &M2),
        operator&=(bool s),
        operator|=(const Matrix<bool> &M2),
        operator|=(bool s);

    // logical operators
    Matrix<bool>
        operator! ()                       const,
        operator&&(const Matrix<bool> &M2) const,
        operator&&(bool s)                 const,
        operator||(const Matrix<bool> &M2) const,
        operator||(bool s)                 const;


    // print matrix to stdout
    void show();
    void print()  {show();}
    void display(){show();}


private:

    // the Matrix' data itself is an ordinary C-array
    bool
        *M;
    // Matrix' size
    unsigned int
        _rows,_cols,_numel;
    bool
        _isempty, _isvector;

    // class of Matrix
    const char*
        _class;
    void
        _get_class(bool&);

    // try to create Matrix from Vector (for casting)
    Matrix<bool>(const Vector<bool> &V2);
};

#endif



// Vector template class
// --------------------------------------------------------------------------
#if 1

template<typename T = double>
class Vector
    : public Matrix<T>
{
public:

    // constructor/destructor
    Vector<T>()                        // empty Vector
        : Matrix<T>()
    {}

    Vector<T>(unsigned int N)          // N-element vector of zeros
        : Matrix<T>(N,1)
    {}

    Vector<T>(unsigned int N,          // N-element vector, with data in contents[]
              T contents[])
        : Matrix<T>(N,contents)
    {}
    Vector<T>(std::vector<T> contents) // data in std::vector<T> contents
        : Matrix<T>(contents)
    {}

    Vector<T>(const Matrix<T> &M)      //
        : Matrix<T>(M)
    {
        if (this->rows != 1 || this->cols != 1)
            throw std::runtime_error("Vector<T>::Vector<T>: Cannot cast matrix to vector.");
    }

    // Vector operations
    Vector<T>
        cross(const Vector<T> &V2) const;

    T
        dot(const Vector<T> &V1) const,
        mag() const,
            norm() const { return mag(); }

    Vector<T>&
        normalize();

    // casts
    operator std::vector<T> () {
        std::vector<T> V;
        for (unsigned int i=0; i< Matrix<T>::numel; ++i)
            V.push_back(this->M[i]);
        return V;
    }

};

#endif // vector class


#endif // class definitions





// Handy typedefs
// --------------------------------------------------------------------------
#if 1

typedef Matrix<double>    dMatrix;
typedef Matrix<float>     fMatrix;
typedef Matrix<int>       iMatrix;
typedef Matrix<uint8_t>   uiMatrix;
typedef Matrix<uint8_t>   ui8Matrix;
typedef Matrix<uint16_t>  ui16Matrix;
typedef Matrix<uint32_t>  ui32Matrix;
typedef Matrix<uint32_t>  ui64Matrix;
typedef Matrix<char>      cMatrix;
typedef Matrix<bool>      bMatrix;

typedef Vector<double>    dVector;
typedef Vector<float>     fVector;
typedef Vector<int>       iVector;
typedef Vector<uint8_t>   uiVector;
typedef Vector<uint8_t>   ui8Vector;
typedef Vector<uint16_t>  ui16Vector;
typedef Vector<uint32_t>  ui32Vector;
typedef Vector<uint32_t>  ui64Vector;
typedef Vector<char>      cVector;
typedef Vector<bool>      bVector;


typedef std::vector<double>        Vec;
typedef std::vector<int>           VecI;
typedef std::vector<unsigned int>  VecUI;


#endif // typedefs





// Matrix base template: constructors
// --------------------------------------------------------------------------
#if 1

// class getter
template<typename T>
void Matrix<T>::_get_class(T&) {
    _class = "unknown";
}

// empty matrix
// TODO: there must be a better way to do this...
#define MATRIX_INIT \
    , _numel     (_rows*_cols) \
    , _isempty   (_numel == 0) \
    , _isvector  ((_rows==1)||(_cols==1))\
    , _isscalar  (_numel == 1) \
                               \
    , rows       (_rows)       \
    , cols       (_cols)       \
    , numel      (_numel)      \
    , datatype   (_class)      \
    , isempty    (_isempty)    \
    , isvector   (_isvector)   \
    , isscalar   (_isscalar)   \
{                              \
    T a; _get_class(a);

template<typename T>
Matrix<T>::Matrix()
    : _rows   (0)
    , _cols   (0)
    MATRIX_INIT
}

// r-by-c matrix of zeros
template<typename T>
Matrix<T>::Matrix(unsigned int r, unsigned int c)
    : _rows      (r)
    , _cols      (c)
    MATRIX_INIT

    // Allocate new M safely
    M = new T[r*c];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    // initailize to zero
    for (unsigned int m=0; m<numel; ++m)
        M[m] = (T)0.0;
}

// r-by-c matrix, with data in contents[]
template<typename T>
Matrix<T>::Matrix(unsigned int r, unsigned int c, T contents[])
    : _rows      (r)
    , _cols      (c)
    MATRIX_INIT

    M = new T[r*c];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    for (unsigned int i=0; i<numel; ++i)
        M[i] = contents[i];
}

// r-by-c matrix, with data in std::vector<T> contents
template<typename T>
Matrix<T>::Matrix(unsigned int r, unsigned int c, std::vector<T> contents)
    : _rows      (r)
    , _cols      (c)
    MATRIX_INIT

    M = new T[r*c];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    for(unsigned int i=0; i<numel; ++i)
        M[i] = contents.at(i); // DO use bounds checking
}

// i-by-i square matrix of zeros
template<typename T>
Matrix<T>::Matrix(unsigned int i)
    : _rows      (i)
    , _cols      (i)
    MATRIX_INIT

    M = new T[i*i];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    // initalize to zero
    for (i=0; i<numel; ++i)
        M[i] = (T)0.0;
}

// i-by-i square matrix, with data in contents[]
template<typename T>
Matrix<T>::Matrix(unsigned int i, T contents[])
    : _rows      (i)
    , _cols      (i)
    MATRIX_INIT

    M = new T[i*i];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    for (i=0; i<numel; ++i)
        M[i] = contents[i];
}

// i-by-i square matrix, with data in std::vector<T> contents
template<typename T>
Matrix<T>::Matrix(unsigned int i, std::vector<T> contents)
    : _rows  (i)
    , _cols  (i)
    MATRIX_INIT

    M = new T[i*i];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    for(i=0; i<numel; ++i)
        M[i] = contents.at(i); // DO check bounds
}

// N-by-1 matrix, from data in std::vector<T> contents
template<typename T>
Matrix<T>::Matrix(std::vector<T> contents)
    : _rows  (contents.size())
    , _cols  (1)
    MATRIX_INIT

    M = new T[_rows];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    for (unsigned int i=0; i<numel; ++i)
        M[i] = contents[i]; // do NOT check bounds
}



// copy matrix of other type
template<typename T> template<typename D>
Matrix<T>::Matrix(const Matrix<D> &M2)
    : _rows (M2.rows)
    , _cols (M2.cols)
    MATRIX_INIT

    M = new T[numel];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");

    // Copy data unsafely
    for (unsigned int m=0; m<numel; ++m)
        M++ = (T)M2.M++;
}


#undef MATRIX_INIT

#endif // Matrix constructors




// Matrix base template: operator overloading
// --------------------------------------------------------------------------
#if 1

// indexing
#if 1

// with subscripts
template<typename T>
T&
Matrix<T>::operator()(unsigned int r, unsigned int c) const
{
    if (r>rows-1 || c>cols-1)
        throw std::runtime_error("Matrix<T>::operator()() - Index out of bounds.");
    return M[r+rows*c];
}

// with subscripts (no bounds check)
template<typename T>
T&
Matrix<T>::operator()(unsigned int r, unsigned int c, bool dummy) const
{
    return M[r+rows*c];
}

// with linear index
template<typename T>
T&
Matrix<T>::operator()(unsigned int ind) const
{
    if (ind>numel-1)
        throw std::runtime_error("Matrix<T>::operator()() - Index out of bounds.");

    return M[ind];
}

// with linear index (no bounds check)
template<typename T>
T&
Matrix<T>::operator()(unsigned int ind, bool dummy) const
{
    return M[ind];
}

// slicing rows, assignment
template<typename T>
auto
Matrix<T>::operator()(Slice r, unsigned int c) -> Cref
{
    return Cref(this->operator()(r,c));
}

// slicing rows, data access
template<typename T>
Matrix<T>
Matrix<T>::operator()(Slice r, unsigned int c) const
{
    if (c > cols-1)
        throw std::runtime_error("Matrix<T>::operator()() - Index out of bounds.");

    unsigned int
        i, k=0,
        L = r.length();
    if (L > rows-1){
        L = 0; for (i=r.start; i<=rows; i+=r.step) ++L;}

    Matrix<T> A(L,1);
    for (i=r.start; i<(r.end>rows?rows:r.end); i+=r.step)
        A.M[k++] = M[i+rows*c];

    return A;
}

// slicing cols
template<typename T>
Matrix<T>
Matrix<T>::operator()(unsigned int r, Slice c) const
{
    if (r > rows-1)
        throw std::runtime_error("Matrix<T>::operator()() - Index out of bounds.");

    unsigned int
        j,k=0,
        L = c.length();
    if (L > cols){
        L = 0; for (j=c.start; j<=cols; j+=c.step) ++L;}

    Matrix<T> A(1,L);
    for (j=c.start; j<(c.end>cols?cols:c.end); j+=c.step)
        A.M[k++] = M[r+rows*j];

    return A;
}

// slicing both
template<typename T>
Matrix<T>
Matrix<T>::operator()(Slice r, Slice c) const
{
    unsigned int
        i,j,k=0,m=0,
        L1 = r.length(),
        L2 = c.length();
    if (L1 > rows){
        L1 = 0; for (i=c.start; i<=rows; i+=r.step) ++L1;}
    if (L2 > cols){
        L2 = 0; for (j=c.start; j<=cols; j+=c.step) ++L2;}

    Matrix<T> A(L1,L2);
    for (j=c.start; j<(c.end>cols?cols:c.end); j+=c.step){
        for (i=r.start; i<(r.end>rows?rows:r.end); i+=r.step)
            A(k++,m,false) = M[i+rows*j];
        ++m; k=0;
    }
    return A;
}

// vectorize
template<typename T>
Matrix<T>
Matrix<T>::operator()(Slice rc) const
{
    unsigned int
        ij,k=0,
        L = rc.length();
    if (L > numel){
        L = 0; for (ij=rc.start; ij<=numel; ij+=rc.step) ++L;}

    Matrix<T> A(L,1);
    for (ij=rc.start; ij<(rc.end>numel?numel:rc.end); ++ij)
        A.M[k++] = M[ij];

    return A;
}

// indexing with boolmatrix
template<typename T>
Matrix<T>
Matrix<T>::operator()(const Matrix<bool> &M2) const
{
    if (M2.rows != rows || M2.cols != cols)
        throw std::runtime_error("Matrix<T>::operator()() - Matrix dimensions must agree.");

    Matrix<T> A(numel,1);
    unsigned int i,j,counter = 0;
    for(i=0; i<M2.numel; ++i){
        if (M2.M[i]){
            A.M[counter] = M[i];
            ++counter;
        }
    }
    if (counter == numel)
        return A;

    else {
        Matrix<T> B(counter,1);
        for (j=0; j<counter; ++j) B.M[j] = A.M[j];
        return B;
    }
}
#endif

// assignment
#if 1

template<typename T>
Matrix<T>&
Matrix<T>::operator=(const Matrix<T> &M2)
{
    if (this != &M2){
        if (numel != 0)
            if ((rows != M2.rows) || (cols != M2.cols))
                throw std::runtime_error("Matrix:(assign): matrix dimensions must be the same.");
        else {
            M = new T[M2.numel];
            _rows = M2.rows; _cols = M2.cols; _numel = M2.numel;
            _isempty = (_numel == 0); _isvector = M2.isvector;
        }
        for (unsigned int i=0; i < M2.numel; ++i)
            M[i] = M2.M[i];
    }
    return *this;
}

#endif


// compound assignments
#if 1

// prefix increment
template<typename T>
Matrix<T>&
Matrix<T>::operator++()
{
    *this += (T)1.0;
    return *this;
}
// suffix increment
template<typename T>
Matrix<T>&
Matrix<T>::operator++(int dummy)
{
    Matrix<T> M3(*this);
    *this += (T)1.0;
    return M3;
}

// prefix decrement
template<typename T>
Matrix<T>&
Matrix<T>::operator--()
{
    *this -= (T)1.0;
    return *this;
}
// suffix decrement
template<typename T>
Matrix<T>&
Matrix<T>::operator--(int dummy)
{
    Matrix<T> M3(*this);
    *this -= (T)1.0;
    return M3;
}

// add-assign
template<typename T>
Matrix<T>&
Matrix<T>::operator+=(const Matrix &M2)
{
    if ((_rows != M2.rows ) || (_cols != M2.cols))
        throw std::runtime_error("Matrix::(plus assignment): matrix dimensions must be the same.");

    for(unsigned int i=0; i<numel; ++i)
        M[i] += M2.M[i];

    return *this;
}
template<typename T>
Matrix<T>&
Matrix<T>::operator+=(T a) {
    for(unsigned int i=0; i<numel; ++i)
        M[i] += a;
    return *this;
}

// subtract-assign
template<typename T>
Matrix<T>&
Matrix<T>::operator-=(const Matrix &M2)
{
    if ((_rows != M2.rows ) || (_cols != M2.cols))
        throw std::runtime_error("Matrix::(plus assignment): matrix dimensions must be the same.");

    for(unsigned int i=0; i<numel; ++i)
        M[i]-= M2.M[i];
    return *this;
}
template<typename T>
Matrix<T>&
Matrix<T>::operator-=(T a) {
    for(unsigned int i=0; i<numel; ++i)
        M[i] -= a;
    return *this;
}

// divide-assign
template<typename T>
Matrix<T>&
Matrix<T>::operator/=(const Matrix &M2)
{
    if ((_rows != M2.rows ) || (_cols != M2.cols))
        throw std::runtime_error("Matrix::(divide assignment): matrix dimensions must be the same.");

    for(unsigned int i=0; i<numel; ++i)
        M[i] /= M2.M[i];
    return *this;
}
template<typename T>
Matrix<T>&
Matrix<T>::operator/=(T a) {
    for(unsigned int i=0; i<numel; ++i)
        M[i] /= a;
    return *this;
}

// multiply-assign
template<typename T>
Matrix<T>&
Matrix<T>::operator*=(T a) {
    for(unsigned int i=0; i<numel; ++i)
        M[i] *= a;
    return *this;
}
template<typename T>
Matrix<T>&
Matrix<T>::operator*=(const Matrix<T> &M2) {
    this = (this->operator*(*this,M2));
    return *this;
}


// exponentiate-assign
template<typename T>
Matrix<T>&
Matrix<T>::operator^=(T a){
    for(unsigned int i=0; i<numel; ++i)
        M[i] = std::pow(M[i],a);
    return *this;
}
template<typename T>
Matrix<T>&
Matrix<T>::operator^=(const Matrix<T> &M2)
{
    if ((rows != M2.rows ) || (cols != M2.cols))
        throw std::runtime_error("Matrix::(exponentiate-assignment): matrix dimensions must be the same.");

    for(unsigned int i=0; i<numel; ++i)
        M[i] = std::pow(M[i], M2.M[i]);
    return *this;
}



#endif


// Subtract, add, divide, multiply
#if 1

// subtract
template<typename T>
Matrix<T>
Matrix<T>::operator-(const Matrix<T> &M2) const {
    return Matrix<T>(*this) -= M2;
}
template<typename T>
Matrix<T>
Matrix<T>::operator-(T a) const {
    return Matrix<T>(*this) -= a;
}
template<typename T>
Matrix<T>
operator-(T a, const Matrix<T> &M2) {
    return -Matrix<T>(M2) += a;
}

// negate
template<typename T>
Matrix<T>
Matrix<T>::operator-() const
{
    Matrix<T> A(*this);
    for(unsigned int i=0; i<numel; i++)
        A.M[i] = -M[i];
    return A;
}

// add
template<typename T>
Matrix<T>
Matrix<T>::operator+(const Matrix<T> &M2) const {
    return Matrix<T>(*this) += M2;
}
template<typename T>
Matrix<T>
Matrix<T>::operator+(T a) const {
    return Matrix<T>(*this) += a;
}
template<typename T>
Matrix<T>
operator+(T a, const Matrix<T> &M2) {
    return M2+a;
}

// "posate"
template<typename T>
Matrix<T> Matrix<T>::operator+() const {
    return *this;
}

// divide
template<typename T>
Matrix<T>
Matrix<T>::operator/(T a) const{
    return Matrix<T>(*this) /= a;
}
template<typename T>
Matrix<T>
operator/(T a, const Matrix<T> &M2) {
    Matrix<T> M1(M2);
    for(unsigned int i=0; i < M2.numel; ++i) M1.M[i] = a/M1.M[i];
    return M1;
}
template<typename T>
Matrix<T>
Matrix<T>::operator/(const Matrix<T> &M2) const {
    return Matrix<T>(*this) /= M2;
}

// multiply
template<typename T>
Matrix<T>
Matrix<T>::operator*(T a) const {
    return Matrix<T>(*this) *= a;
}
template<typename T>
Matrix<T>
operator*(T a, const Matrix<T> &M2) {
    return Matrix<T>(M2) *= a;
}

// Matrix multiplication
template<typename T>
Matrix<T>
Matrix<T>::operator*(const Matrix<T> &M2) const
{
    if ( (isempty) || (M2.isempty) )
        throw std::runtime_error("Matrix::times() - Product of empty matrix is undefined.");
    if (cols != M2.rows)
        throw std::runtime_error("Matrix::times() - Matrix dimensions do not agree.");

    Matrix<T> M3(rows,M2.cols);
    unsigned int i,j,k;
    for(i = 0; i < rows; ++i)
        for(j = 0; j < M2.cols; ++j)
            for(k = 0; k < cols; ++k)
                M3(i,j) += this->operator()(i,k,false) * M2(k,j);
    return M3;
}

// exponentiate
template<typename T>
Matrix<T>
Matrix<T>::operator^(T a) const{
    return Matrix<T>(*this) ^= a;
}
template<typename T>
Matrix<T>
operator^(T a, const Matrix<T> &M2) {
    Matrix<T> M1(M2);
    for(unsigned int i=0; i<M2.numel; ++i)
        M1.M[i] = std::pow(a, M1.M[i]);
    return M1;
}
template<typename T>
Matrix<T>
Matrix<T>::operator^(const Matrix<T> &M2) const {
    return Matrix<T>(*this) ^= M2;
}



#endif

// Comparisons
#if 1

template<typename T>
Matrix<bool>
Matrix<T>::operator==(const Matrix<T> &M2) const
{
    if ((M2.cols != cols) || (M2.rows != rows))
        throw std::runtime_error("Matrix::eq() - Matrix dimensions must agree.");

    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] == M2.M[i]);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator==(T a) const
{
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] == a);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator<(const Matrix<T> &M2) const
{
    if ((M2.cols != cols) || (M2.rows != rows))
        throw std::runtime_error("Matrix::lt() - Matrix dimensions must agree.");
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] < M2.M[i]);
    return A;
}
template<typename T>
Matrix<bool>
Matrix<T>::operator<(T a) const
{
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] < a);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator<=(const Matrix<T> &M2) const
{
    if ((M2.cols != cols) || (M2.rows != rows))
        throw std::runtime_error("Matrix::le() - Matrix dimensions must agree.");

    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] <= M2.M[i]);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator<=(T a) const
{
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] <= a);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator>(const Matrix<T> &M2) const
{
    if ((M2.cols != cols) || (M2.rows != rows))
        throw std::runtime_error("Matrix::gt() - Matrix dimensions must agree.");

    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] > M2.M[i]);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator>(T a) const
{
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] > a);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator>=(const Matrix<T> &M2) const
{
    if ((M2.cols != cols) || (M2.rows != rows))
        throw std::runtime_error("Matrix::ge() - Matrix dimensions must agree.");

    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] >= M2.M[i]);
    return A;
}

template<typename T>
Matrix<bool>
Matrix<T>::operator>=(T a) const
{
    Matrix<bool> A(rows,cols);
    for(unsigned int i=0; i<numel; ++i)
        A.M[i] = (M[i] >= a);
    return A;
}

#endif

// casting
#if 1

// FIXME: DOES NOT WORK

//// cast scalars to matrix type
//template <>
//operator Matrix<> (double p) {
    //return Matrix<>(1,1, &p);
//}

//template <>
//operator Matrix<int>  (int p) {
    //return Matrix<int>(1,1, &p);
//}


#endif // casting

#endif // Matrix operator overloading




// Matrix base template: basic ops
// --------------------------------------------------------------------------
#if 1

// Transpose
template<typename T>
Matrix<T>
Matrix<T>::operator~()
{
    // TODO: there must be a smarter way to do this...
    Matrix A(cols,rows);
    for (unsigned int i=0; i<rows; ++i)
        for (unsigned int j=0; j<cols; ++j)
            A(j,i,false) = this->operator()(i,j,false);
    return A;
}
template<typename T>
Matrix<T>
Matrix<T>::transpose() {
    return this->operator~();
}

// Matrix trace
template<typename T>
T
Matrix<T>::trace() const
{
    if (rows != cols)
        throw std::runtime_error("Matrix<T>::trace() - Matrix must be square.");

    T t = 0.0;
    for (unsigned int i=0; i<rows; ++i)
        t += this->operator()(i,i,false);
    return t;
}

// Matrix diagonal
template<typename T>
Matrix<T>
Matrix<T>::diag() const
{
    if (rows != cols)
        throw std::runtime_error("Matrix<T>::diag() - Matrix must be square.");

    Matrix<T> A(rows,1);
    for (unsigned int i=0; i<rows; ++i)
        A.M[i] = this->operator()(i,i,false);

    return A;
}

// Matrix size
template<typename T>
Matrix<unsigned int>
Matrix<T>::size() const {
    Matrix<unsigned int> M(1,2);
    M(1) = rows;
    M(2) = cols;
    return M;
}

#endif // Matrix basic ops





// Vector base template: vector-specific operations
// --------------------------------------------------------------------------
#if 1

// cross product
template<typename T>
Vector<T>
Vector<T>::cross(const Vector &V2) const
{
    if ( (V2.numel != 3) || (this->numel != 3) )
        throw std::runtime_error("Vector::cross() - Cross product is only defined for two 3-element Matrices.");

    Vector<T> A(V2);
    A.M[0] = V2.M[1]*this->M[2] - V2.M[2]*this->M[1];
    A.M[1] = V2.M[2]*this->M[0] - V2.M[0]*this->M[2];
    A.M[2] = V2.M[0]*this->M[1] - V2.M[1]*this->M[0];
    return A;
}
template<typename T>
T
cross(const Vector<T> &V1, const Vector<T> &V2) {
    return V1.cross(V2);
}


// dot product
template<typename T>
T
Vector<T>::dot(const Vector &V2) const
{
    if (V2.numel != this->numel)
        throw std::runtime_error("Vector<T>::dot() - Vectors must have an equal number of elements.");

    T dp = 0.0;
    for (unsigned int i=0; i<this->numel; ++i)
        dp += V2.M[i]*this->M[i];

    return dp;
}
template<typename T>
T
dot(const Vector<T> &V1, const Vector<T> &V2) {
    return V1.dot(V2);
}


// magnitude
template<typename T>
T
Vector<T>::mag() const {
    return (T)sqrt(dot(*this));
}
template<typename T>
T
mag(const Vector<T> &V1) {
    return V1.mag();
}
template<typename T>
T
norm(const Vector<T> &V) {
    return V.mag();
}


// normalize vector
template<typename T>
Vector<T>&
Vector<T>::normalize()
{
    // error check is done in mag()
    T VM = mag();
    if (VM == (T)0.0)
        throw std::runtime_error("Vector::normalize() - magnitude is zero.");

    for (unsigned int i=0; i<this->numel; ++i)
        this->M[i] /= VM;

    return *this;
}
template<typename T>
Vector<T>
normalize(const Vector<T> &V) {
    return Vector<> (V).normalize();
}


#endif // Vector-specific operations




// Matrix base template: Special ops
// --------------------------------------------------------------------------
#if 1

// all
template<typename T>
bool
all(const Matrix<T> &M1, unsigned int direction = 1)
{
    for (unsigned int i=0; i<M1.numel; ++i)
        if (!M1(i,false))
            return false;
    return true;
}

// any
template<typename T>
bool
any(const Matrix<T> &M1, unsigned int direction = 1)
{
    for (unsigned int i=0; i<M1.numel; ++i)
        if (M1(i,false))
            return true;
    return false;
}

// isnan, isfinite, isinf
template<typename T>
Matrix<bool>
isfinite(const Matrix<T> &M){
    Matrix<bool> B(M.rows,M.cols);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::isfinite(M(i,false));
    return B;
}

template<typename T>
Matrix<bool>
isnan(const Matrix<T> &M){
    Matrix<bool> B(M.rows,M.cols);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::isnan(M(i,false));
    return B;
}

template<typename T>
Matrix<bool>
isinf(const Matrix<T> &M){
    Matrix<bool> B(M.rows,M.cols);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::isinf(M(i,false));
    return B;
}



// print matrix to stdout
template<typename T>
void
Matrix<T>::show()
{
    if (isempty)
        std::cout << "[]" << std::endl;

    else {
        unsigned int i,j;
        for (i=0; i<rows; ++i){
            std::cout << std::endl;
            for (j=0; j<cols; ++j)
                std::cout << M[i+rows*j];
        }
    }
    std::cout << std::endl << std::endl;
}

#endif // special ops




 Basic Matrix generators
// --------------------------------------------------------------------------
#if 1

// zeros, ones
template<typename T = double>
Matrix<T>
ones(unsigned int r, unsigned int c) {
    Matrix<T> A(r,c);
    for (unsigned int i=0; i<A.numel; ++i)
        A(i,false) = (T)1.0;
    return A;
}
template<typename T = double>
Matrix<T>
ones(unsigned int i) {
    Matrix<T> A(i,i);
    for (i=0; i<A.numel; ++i)
        A(i,false) = (T)1.0;
    return A;
}

template<typename T = double>
Matrix<T>
zeros(unsigned int r, unsigned int c) {
    return Matrix<T>(r,c);
}
template<typename T = double>
Matrix<T>
zeros(unsigned int i) {
    return Matrix<T>(i,i);
}


// identity matrix
template<typename T = double>
Matrix<T>
eye(unsigned int i) {
    Matrix<T> A(i,i);
    for (i=0; i<A.rows; ++i)
        A(i,i,false) = (T)1.0;
    return A;
}

// Matrix trace
template<typename T>
T
trace(const Matrix<T> &A) {
    return A.trace();
}

// Matrix diagonal
// TODO: there's more functionality to make here...
template<typename T>
Matrix<T>
diag(const Matrix<T> &A) {
    return A.diag();
}

// horizontal concatenation of matrices
template<typename T>
Matrix<T>
hcat(const Matrix<T> &A, const Matrix<T> &B)
{
    if (A.rows != B.rows)
        throw std::runtime_error("Matrix<T>::hcat() - Matrix dimensions must agree.");

    Matrix<T> C(A.rows, B.cols+A.cols);
    for (unsigned int i=0; i<A.numel; ++i)
        C.M[i] = A.M[i];
    for (unsigned int j=0; j<B.numel; ++j)
        C.M[j+A.numel] = B.M[j];

    return C;
}

// vertical concatenation of matrices
template<typename T>
Matrix<T>
vcat(const Matrix<T> &A, const Matrix<T> &B)
{
    if (A.cols != B.cols)
        throw std::runtime_error("Matrix<T>::vcat() - Matrix dimensions must agree.");

// FIXME: THIS IS INCORRECT
    Matrix<T> C(A.rows+B.rows, A.cols);
    for (unsigned int i=0; i<A.numel; ++i)
        C.M[i] = A.M[i];
    for (unsigned int j=0; j<B.numel; ++j)
        C.M[j+A.numel] = B.M[j];

    return C;
}

// min, max, abs, etc.
template<typename T>
T
min(T a, T b){
    return (a<b ? a : b);
}
template<typename T>
Matrix<T>
min(const Matrix<T> &A, const Matrix<T> &B = {}, unsigned int dimension = 1)
{
    // TODO: row/column wise

    T
        mm = (T)INFINITY;

// FIXME: NOT DONE YET

    Matrix<> sz = A.size();
    sz(dimension) = 1;
    Matrix<> C(sz(1), sz(2));

    for (unsigned int r=0; r<A.rows; ++r){
        if (dimension != 1) {
            mm = (T)INFINITY;
        }
        for (unsigned int c=0; c<A.cols; ++c){
            if (dimension != 2) {
                mm = (T)INFINITY;
            }
            T tmp = A(r,c,false);
            if (tmp < mm)
                mm = tmp;
        }
    }
    return mm;
}

template<typename T>
T
max(T a, T b){
    return (a<b ? b : a);
}
template<typename T>
Matrix<T>
max(const Matrix<T> &A)
{
// FIXME: NOT DONE YET

    T
        MM = (T)(-INFINITY),
        tmp;
    for (unsigned int r=0; r<A.rows; ++r){
        for (unsigned int c=0; c<A.cols; ++c){
            tmp = A(r,c,false);
            if (tmp > MM) MM = tmp;
        }
    }
    return MM;
}

template<typename T>
T
abs(T val){
    return (val<(T)0 ? -val : +val);
}
template<typename T>
Matrix<T>
abs(const Matrix<T> &A) {
    Matrix<> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = abs(A(i,false));
    return B;
}


template<typename T>
int8_t
sign(T val){
    return (int)(val>(T)0)-(val<(T)0);
}
template<typename T>
Matrix<int8_t>
sign(const Matrix<T> &A)
{
    Matrix<int8_t> B(A.rows,A.cols);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = sign(B(i,false));
    return B;
}


// trig
template<typename T>
Matrix<T>
cos(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::cos(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
sin(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::sin(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
tan(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::tan(B(i,false));
    return B;
}

// TODO: sec, csc, cot

template<typename T>
Matrix<T>
acos(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::acos(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
asin(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::asin(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
atan(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::atan(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
atan2(const Matrix<T> &A1, const Matrix<T> &A2)
{
    if ( (A1.numel != A2.numel && A1.numel != 1 && A2.numel != 1) ||
         (A1.rows != A2.rows) )
        throw std::runtime_error("atan2: Matrix dimensions must agree.");

    Matrix<T> B(A1);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::atan2(A2(i,false), B(i,false));
    return B;
}

// TODO: asec, acsc, acot

// TODO: cosh, sinh, tanh, sech, csch, coth
// TODO: acosh, asinh, atanh, asech, acsch, acoth


// exponentiation
template<typename T>
Matrix<T>
log(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::log(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
exp(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::exp(B(i,false));
    return B;
}

template<typename T>
Matrix<T>
pow(const Matrix<T> &A, const double pwr) {
    return Matrix<T>(A).operator^(pwr);
}

template<typename T>
Matrix<T>
sqrt(const Matrix<T> &A) {
    Matrix<T> B(A);
    for (unsigned int i=0; i<B.numel; ++i)
        B(i,false) = std::sqrt(B(i,false));
    return B;
}


#endif // matrix generators



#endif // _MATRIX_HPP






