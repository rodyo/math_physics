#include "Matrix.hpp"


/*****************************************************************************
 Slice class (for array slicing into the Matrix)
*****************************************************************************/
#if 1

#define SLICE_INIT \
    , start  (_start) \
    , end    (_end)   \
    , all    (_all)   \
    , step   (_step)

// empty constructor: slice over all rows/columns
template <typename T>
Matrix<T>::Slice::Slice()
    : _start (0)
    , _step  (1)
    , _end   (0)
    , _all   (true)
    SLICE_INIT
{}

// use only start, end (step = 1)
template <typename T>
Matrix<T>::Slice::Slice(unsigned int start, unsigned int end)
    : _start (start)
    , _step  (1)
    , _end   (end)
    , _all   (false)
    SLICE_INIT
{}

// full call: start,step,end
template <typename T>
Matrix<T>::Slice::Slice(unsigned int start, int step, unsigned int end)
    : _start (start)
    , _step  (step)
    , _end   (end)
    , _all   (false)
    SLICE_INIT
{ }

// return number of elements inferred by the slice
template <typename T>
unsigned int
Matrix<T>::Slice::length() {
    return ( end==INT_MAX ? INT_MAX : (int)floor((end-start)/step)+1 );
}


/*****************************************************************************
 Cref class
 * for assignmetns done with slices or int/boolmatrices
*****************************************************************************/

union Matrix<T>::Cref::rowind_t
{
    const Slice R;
    const unsigned int r;
    const Matrix<unsigned int> iR;
    const Matrix<bool> I;
};

union Matrix<T>::Cref::colind_t
{
    const Slice C;
    const unsigned int c;
    const Matrix<unsigned int> iC;
};

const Matrix<T> &
Matrix<T>::Cref::operator=(const Matrix<T> &M2)
{
    if ( (M2.numel != M.numel) || (M2.rows != M.rows) )
        throw std::runtime_error("Subscripted assignment dimension mismatch.");

// TODO

    switch (mode)
    {
        // Mode 0: sliced row, constant column
        case 0:
            break;

        // Mode 1: constant row, sliced column
        case 1:
            break;

        // Mode 2: both indices sliced
        case 2:
            break;

        // Mode 3: linear slice
        case 3:
            break;

        // Mode 4: boolmatrix
        case 4:
            break;

        // Mode 5: row with int matrix, constant col
        case 5:
            break;

        // Mode 6: constant row, col with int matrix
        case 6:
            break;

        // Mode 7: linear intmatrix
        case 7:
            break;
    }

    return M;
}



// Mode 0: sliced row, constant column
Matrix<T>::Cref::Cref(Matrix<T> &M, const Slice &R, unsigned int c)
    : M   (M)
    , mode(0) {
    rowind.R = R;
    colind.c = c;
}

// Mode 1: constant row, sliced column
Matrix<T>::Cref::Cref(Matrix<T> &M, unsigned int r, const Slice &C)
    : M   (M)
    , mode(1) {
    rowind.r = r;
    colind.C = C;
}

// Mode 2: both indices sliced
Matrix<T>::Cref::Cref(Matrix<T> &M, const Slice &R, const Slice &C)
    : M   (M)
    , mode(2) {
    rowind.R = R;
    colind.C = C;
}

// Mode 3: linear slice
Matrix<T>::Cref::Cref(Matrix<T> &M, const Slice &I)
    : M   (M)
    , mode(3) {
    rowind.R = I;
}

// Mode 4: boolmatrix
Matrix<T>::Cref::Cref(Matrix<T> &M, const Matrix<bool> &I)
    : M   (M)
    , mode(3) {
    rowind.I = I;
}


#endif



/*****************************************************************************
 Matrix, class getters
*****************************************************************************/
#if 1

template<> void Matrix<double>::_get_class(double&) { _class = "double";}
template<> void Matrix<float>::_get_class (float&)  { _class = "float"; }
template<> void Matrix<int>::_get_class   (int&)    { _class = "int";   }
template<> void Matrix<char>::_get_class  (char&)   { _class = "char";  }
           void Matrix<bool>::_get_class  (bool&)   { _class = "bool";  }

#endif // class getters



/*****************************************************************************
 Matrix bool specialization:
*****************************************************************************/
#if 1

#endif

/*****************************************************************************
 Matrix<bool> specialization: constructors
*****************************************************************************/
#if 1

// empty bool matrix
Matrix<bool>::Matrix() :
    rows(_rows), cols(_cols), numel(_numel), isempty(_isempty), isvector(_isvector)
{
     _rows = 0; _cols = 0; _numel = 0;
    _isempty = true; _isvector = false;
}

// r-by-c bool matrix of false's
Matrix<bool>::Matrix(unsigned int r, unsigned int c) :
    rows(_rows), cols(_cols), numel(_numel), isempty(_isempty), isvector(_isvector)
{
    M = new bool[r*c];
    if (!M)
        throw std::runtime_error("Matrix<bool>::Matrix() - Matrix could not be created.");
    _rows = r; _cols = c; _numel = r*c;
    _isempty = (_numel == 0); _isvector = ((_rows==1)||(_cols==1));
    for(unsigned int m = 0; m < numel; m++) M[m] = false;
}

Matrix<bool>::Matrix(unsigned int i) :
    rows(_rows), cols(_cols), numel(_numel), isempty(_isempty), isvector(_isvector)
{
    M = new bool[i*i];
    if (!M)
        throw std::runtime_error("Matrix<bool>::Matrix() - Matrix could not be created.");
    _rows = i; _cols = i; _numel = i*i;
    _isempty = (i==0); _isvector = (i==1);
    for(unsigned int m = 0; m < numel; m++) M[m] = false;
}

Matrix<bool>::Matrix(const Matrix<bool> &M2) :
    rows(_rows), cols(_cols), numel(_numel), isempty(_isempty), isvector(_isvector)
{
    M = new bool[M2.numel];
    if (!M)
        throw std::runtime_error("Matrix<T>::Matrix() - Matrix could not be created.");
    _rows = M2.rows; _cols = M2.cols; _numel = M2.numel;
    _isempty = (_numel == 0); _isvector = ((_rows==1)||(_cols==1));
    for(unsigned int m = 0; m < numel; m++) M[m] = M2.M[m];
}

// destructor
Matrix<bool>::~Matrix()
    {delete[] M;}

#endif



/*****************************************************************************
 Matrix<bool> specialization: operator overloads
*****************************************************************************/
#if 1

// assignment
Matrix<bool>& Matrix<bool>::operator=(const Matrix<bool> &M2)
{
    if (this != &M2){
        if ((rows != M2.rows ) || (cols != M2.cols))
            throw std::runtime_error("Matrix:(assign): matrix dimensions must be the same.");
        for(unsigned int i=0; i<M2.numel; i++) M[i] = M2.M[i];
    }
    return *this;
}

// indexing with subscripts
bool& Matrix<bool>::operator()(unsigned int i, unsigned int j) const
{
    if (i>rows-1 || j>cols-1)
        throw std::range_error("Matrix<bool>::operator()() - Index out of bounds.");
    return M[i+rows*j];
}
// indexing with linear index
bool& Matrix<bool>::operator()(unsigned int i) const
{
    if (i>numel-1)
        throw std::range_error("Matrix<bool>::operator()() - Index out of bounds.");
    return M[i];
}
// indexing with boolmatrix
Matrix<bool> Matrix<bool>::operator()(const Matrix<bool> &M2) const
{
    if (M2.rows != rows || M2.cols != cols)
        throw std::range_error("Matrix<bool>::operator()() - Matrix dimensions must agree.");
    Matrix<bool> A(numel,1);
    unsigned int i,j,counter=0;
    for(i=0; i<M2.numel; i++){
        if (M2.M[i]){
            counter++;
            A.M[i] = M[i];
        }
    }
    if (counter == numel)
        return A;
    else {
        Matrix<bool> B(counter,1);
        for(j=0; j<counter; j++) B.M[j] = A.M[j];
        return B;
    }
}

#endif



/*****************************************************************************
Matrix bool specialization: Basic ops
*****************************************************************************/
#if 1

// print matrix
void Matrix<bool>::show()
{
    for (unsigned int i=0; i<rows; i++){
        printf("\n");
        for (unsigned int j=0; j<cols; j++)
            printf("%d ", (int)this->operator()(i,j));
    }
    printf("\n\n");
}


//bool all(const Matrix<bool> &M1)
//{
    //for(unsigned int r=0; r<M1.rows; r++)
        //for(unsigned int c=0; c<M1.cols; c++)
            //if (M1(r,c,false) == 0.0) return false;
    //return true;
//}

#endif




