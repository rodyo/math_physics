#include "Matrix.hpp"







//
///*****************************************************************************
// Matrix, class getters
//*****************************************************************************/
//
//
//template<> void Matrix<double>::_get_class(double&) { _class = "double";}
//template<> void Matrix<float>::_get_class (float&)  { _class = "float"; }
//template<> void Matrix<int>::_get_class   (int&)    { _class = "int";   }
//template<> void Matrix<char>::_get_class  (char&)   { _class = "char";  }
//           void Matrix<bool>::_get_class  (bool&)   { _class = "bool";  }
//



/*****************************************************************************
 Matrix bool specialization:
*****************************************************************************/


/*****************************************************************************
 Matrix<bool> specialization: constructors
*****************************************************************************/


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





/*****************************************************************************
 Matrix<bool> specialization: operator overloads
*****************************************************************************/


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





/*****************************************************************************
Matrix bool specialization: Basic ops
*****************************************************************************/

//// print matrix
//void Matrix<bool>::show()
//{
//    for (unsigned int i=0; i<rows; i++){
//        printf("\n");
//        for (unsigned int j=0; j<cols; j++)
//            printf("%d ", (int)this->operator()(i,j));
//    }
//    printf("\n\n");
//}


//bool all(const Matrix<bool> &M1)
//{
    //for(unsigned int r=0; r<M1.rows; r++)
        //for(unsigned int c=0; c<M1.cols; c++)
            //if (M1(r,c,false) == 0.0) return false;
    //return true;
//}






