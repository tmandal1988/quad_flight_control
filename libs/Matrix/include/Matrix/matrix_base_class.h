#ifndef MATRIXBASE_H
#define MATRIXBASE_H
#endif

#include<iostream>
#include<string>

using namespace std;

template <typename T>
class MatrixBase{
	
	public:
		// constructors
		MatrixBase();
		MatrixBase(size_t num_rows, size_t num_cols);
		MatrixBase(size_t num_rows, size_t num_cols, string type);
		~MatrixBase();

		// copy constructor
		MatrixBase(const MatrixBase<T>& matrix_to_copy);
		// member functions
		void PrintMatrix();
		// operator overloading
		MatrixBase operator= (const MatrixBase<T>& matrix_to_copy);
		// assign elements at row idx and col idx
		T& operator() (size_t r_idx, size_t c_idx);
		// read elements at row idx and col idx
		T operator() (size_t r_idx, size_t c_idx) const;
		// multiply two matrices
		MatrixBase operator* (const MatrixBase& matrix_to_multipy);
		// add two matrices
		MatrixBase operator+ (const MatrixBase& matrix_to_add);
		// subtract two matrices
		MatrixBase operator- (const MatrixBase& matrix_to_subtract);
		// add a scalar
		MatrixBase operator+ (const T x);
		// subtract a scalar
		MatrixBase operator- (const T x);
		// multiply a scalar
		MatrixBase operator* (const T x);
		// divide by a scalar
		MatrixBase operator/ (const T x);

	private:
		size_t nrows_, ncols_;
		T **matrix_;


};