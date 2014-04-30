#include "math.h"

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Matrix class definition
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

 class Matrix{
 	
	private:
	
		// data members
		unsigned rows_, cols_;    // Matrix size
		double*  data_;            // matrix contents
		bool     transposed_;      // Is the matrix transposed?
		int      indexBase_;       // Base of indexing
		char     type_;            // type of the data (c = complex double, {d = double}, i = int, b = bool)
		
		// private members
		void checkindex(unsigned row, unsigned col) const;
		void checkindex(unsigned element) const;

	// all member functions
	public:
	
		// Matrix operators
		double& operator() (unsigned row, unsigned col);       // subscript assignment   
		double  operator() (unsigned row, unsigned col) const; // subscript reference
		double& operator() (unsigned element);                 // linear indexing assignment
		double  operator() (unsigned element) const;           // linear indexing reference
		
		Matrix& operator=  (const Matrix& m);       void operator~  (void);     // Assignment operator / transpose
		Matrix& operator*  (const Matrix& B);    Matrix& operator*  (double b); // multiply
		  void  operator*= (const Matrix& B);      void  operator*= (double b); // multiply-assign
		Matrix& operator+  (const Matrix& B);    Matrix& operator+  (double b); // add
		  void  operator+= (const Matrix& B);      void  operator+= (double b); // add-assign
		Matrix& operator-  (const Matrix& B);    Matrix& operator-  (double b); //subtract
		  void  operator-= (const Matrix& B);      void  operator-= (double b); // subtract-assign
		Matrix& operator^  (const double p);                                    // power
		
		// Comparison
		Matrix& operator<  (const Matrix& B);    Matrix& operator<  (double b); // less than
		Matrix& operator<= (const Matrix& B);    Matrix& operator<= (double b); // less than or equals
		Matrix& operator>  (const Matrix& B);    Matrix& operator>  (double b); // greater than
		Matrix& operator>= (const Matrix& B);    Matrix& operator>= (double b); // greater than or equals
		Matrix& operator== (const Matrix& B);    Matrix& operator== (double b); // equals exactly
		Matrix& operator!= (const Matrix& B);    Matrix& operator!= (double b); // does not equal exactly
		Matrix& operator&& (const Matrix& B);    Matrix& operator&& (double b); // AND (for logical matrices)
		Matrix& operator|| (const Matrix& B);    Matrix& operator|| (double b); // OR  (for logical matrices)
			   
		// sets/gets
		double rows(){return rows_;}
		double cols(){return cols_;}
		void   setIndexbase(unsigned int b){indexBase_ = b;}
	   
		// more member functions
		Matrix(unsigned rows, unsigned cols);       // Constructor
		~Matrix();                                  // Destructor
		Matrix(const Matrix& m);                    // Copy matrix
		
		
 };
 
 // constructor
 //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 inline Matrix::Matrix(unsigned rows, unsigned cols) : rows_ (rows), cols_ (cols)   
 {
   // error trap
   if (rows == 0 || cols == 0)
     throw BadIndex("Matrix has 0 size.");
	// initialize data
    data_ = new double[rows * cols];
	// initialize everything else
	transposed = false;
	indexBase = 1;
 }
 
 // destructor
 //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 inline Matrix::~Matrix(){delete[] data_;}
 
 // Check subscripts
 void Matrix::checkindex(unsigned row, unsigned col) const
 {
	if (row >= rows_ || col >= cols_)
		throw BadIndex("Subscripts out of range."); 
 }
 // Check linear index
 void Matrix::checkindex(unsigned element) const;
 {
	if ( element >= (rows_ * cols_))
		throw BadIndex("Linear index out of range."); 
 }
 
 // subscript assignment 
 //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 inline double& Matrix::operator() (unsigned row, unsigned col)
 {
	// Check indices
	checkindex(row,col);
	// return proper index
	if transposed
		return data_[cols_*(col-indexBase) + (row-indexBase)];
	else
		return data_[cols_*(row-indexBase) + (col-indexBase)];		
 }
 
 // subscript reference 
 //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 inline double Matrix::operator() (unsigned row, unsigned col) const
 {
   // Check indices
	checkindex(row,col);
   // return proper index
   if transposed
		return data_[cols_*(col-indexBase) + (row-indexBase)];
   else
		return data_[cols_*(row-indexBase) + (col-indexBase)];
 }
 
 // linear indexing assignment
 inline double& Matrix::operator() (unsigned element) 
 {
	// Check indices
	checkindex(element);
	// return proper index
	if transposed
		return data_[element-indexBase];
	else
		return data_[element-indexBase];
 }
 
 // linear indexing reference
 inline double Matrix::operator() (unsigned element) const
 {
	// Check indices
	checkindex(element);
	// return proper index
	if transposed
		return data_[element-indexBase];
	else
		return data_[element-indexBase];
 }
 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Matrix properties
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// number of elements in a matrix
inline unsigned int numel(const Matrix& A){return (A.rows*A.cols);}

// Is there any element nonzero/true?
inline bool any(const Matrix& A)
{
   for (int i=0; i<numel(A); i++)
      {if (A(i) != 0.0) return true;}
}

// Is the Matrix empty?
inline bool isempty(Matrix A)
   {return (A.getCols==0 || A.getRows==0)};

// Are all elements nonzero/true?
bool all(const Matrix& A)
{
   for (int i=0; i<numel(A); i++)
      {if (A(i) == 0.0) return false;}
   return true;
}

// maximum element
double max(const Matrix& A)
{	
	double M = -inf;	
	for (int i=0; i<numel(A); i++)
      {if (A(i) > M) M = A(i);}
	return M;
}

// minimum element
double min(const Matrix& A)
{
	double m = +inf;	
	for (int i=0; i<numel(A); i++)
      {if (A(i) < m) m = A(i);}
	return m;
}

// number of nonzero elements
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
inline unsigned int nnz(const Matrix& A)
{
    unsigned int nz = 0;
    for (int i=0; i<numel(A); i++)
       {if (A(i)!=0.0) nz+=1;}
    return nz;
}

// reshape Matrix
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& reshape(const Matrix& A, int i, int j);
Matrix& reshape(const Matrix& A, int i, char j);
Matrix& reshape(const Matrix& A, char i, int j);

// Matrix of nonzero-elements
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& nonzeros(const Matrix& A);


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Standard matrices
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Matrix of ones
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& ones(unsigned rows, unsigned cols)
{
	Matrix A = Matrix(rows,cols);
	for (int i=0; i<numel(A&); i++)
       {A(i)=1.0;}
	return A;
}
Matrix& ones(unsigned squaresize)
{
	return ones(squaresize,squaresize);
}

// Matrix of zeros
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix zeros(unsigned rows, unsigned cols)
{
	Matrix A = Matrix(rows,cols);
	for (int i=0; i<numel(A&); i++)
       {A(i)=0.0;}
	return A;
}
Matrix zeros(unsigned squaresize);
{
	return zeros(squaresize,squaresize);
}

// Matrix of inf's
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix inf(unsigned rows, unsigned cols)
{
	Matrix A = Matrix(rows,cols);
	for (int i=0; i<numel(A&); i++)
       {A(i)=inf;}
	return A;
}
Matrix inf(unsigned squaresize)
{
	return inf(squaresize,squaresize);
}

// Matrix of NaN's
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix NaN(unsigned rows, unsigned cols)
{
	Matrix A = Matrix(rows,cols);
	for (int i=0; i<numel(A&); i++)
       {A(i)=NaN;}
	return A;
}
Matrix NaN(unsigned squaresize)
{
	return NaN(squaresize,squaresize);
}

// Matrix of random doubles
Matrix rand(unsigned rows, unsigned cols);
Matrix rand(unsigned squaresize);

// Matrix of random integers
Matrix randi(unsigned rows, unsigned cols);
Matrix randi(unsigned squaresize);

// Matrix of normally distributed random doubles
Matrix randn(unsigned rows, unsigned cols);
Matrix randn(unsigned squaresize);

// Matrix of t-student distributed random doubles
Matrix randt(unsigned rows, unsignedcols);
Matrix randt(unsigned squaresize);

// Matrix of Chi-square distributed random doubles
Matrix randchi(unsigned rows, unsigned cols);
Matrix randchi(unsigned squaresize);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Elementary matrix operations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// Matrix of eigenvalues
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& eig(const Matrix& A);

// Matrix of eigenvectors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& eigv(const Matrix& A);

// Matrix inverse
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& inv(const Matrix& A);

// Singular value decomposition
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& svd(const Matrix& A);

// QR decomposition
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& QR(const Matrix& A);

// LU decomposition
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& QR(const Matrix& A);

// Kronecker product
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& kron(const Matrix& A, const Matrix& B);


// Find
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix find(const Matrix& B){
       for (register int i=0; i < B.numel(); i++){
           }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Overloaded operators
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// transpose matrix 
void Matrix::operator~()
{
	// swap rows and cols
	int tmp=rows_; rows_=cols_; cols_=tmp;
	// set transposed to inverse value
	transposed_ = !transposed_;
}




// matrix-matrix (or matrix-vector) multiply
Matrix Matrix::operator*(Matrix B){

       // create new Matrix
       Matrix C(rows, B.getCols);

       for (int i=0; i<*rows; i++){
           for (int j=0; j<*cols; j++){
               contents++;
               
}
// matrix-scalar multiply
Matrix Matrix::operator*(double b){
       for (register int i=0; i<this.numel(); i++)
           contents[i] *= b;
       return *this;
}
// Same things, but then for "*="
//Matrix Matrix::operator*=(Matrix B){
//       for (int i=0; i<(*rows * *cols); i++)
//           contents[i] *= b;
//}
Matrix Matrix::operator*=(double b){
       for (int i=0; i<(*rows * *cols); i++)
           contents[i] *= b;
}

// matrix-matrix (or matrix-vector) addition
//Matrix Matrix::operator+(Matrix B){
//}
// matrix-scalar addition
Matrix Matrix::operator+(double b){
       for (int i=0; i<(*rows * *cols); i++)
           contents[i] += b;
}

// matrix-matrix (or matrix-vector) subtraction
/*
Matrix Matrix::operator-(Matrix B){
       // first check sizes
       // if ( rows != B.getRows) && (cols != B.getCols) )

       // if all is well, execute the subtraction
       for (int i=0; i<rows*cols; i++)
           contents[i] += b;
}
*/
// matrix-scalar subtraction
Matrix Matrix::operator-(double b){
   for (int i=0; i<(*rows * *cols); i++)
      contents[i] -= b;
   return *this;
}

// matrix power
Matrix Matrix::operator^(double b){
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Elementary functions for Matrix input
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

// square root of a matrix: element wise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& sqrt(const Matrix& A){
       // copy matrix
       Matrix B = Matrix(A);
       // loop through A's contents
       for (int element = 0; element < numel(B); element++)
			{B(element) = std::sqrt(A(element));}
       // done
       return B;
}

// exponential of a matrix: element wise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& exp(const Matrix A);

// trigonometrics: all element wise
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Matrix& sin(Matrix A);
Matrix& cos(Matrix A);
Matrix& tan(Matrix A);
Matrix& sin(Matrix A);
Matrix& cos(Matrix A);
Matrix& tan(Matrix A);
Matrix& csc(Matrix A);
Matrix& sec(Matrix A);
Matrix& cot(Matrix A);
Matrix& sinh(Matrix A);
Matrix& cosh(Matrix A);
Matrix& tanh(Matrix A);
Matrix& csch(Matrix A);
Matrix& sech(Matrix A);
Matrix& coth(Matrix A);