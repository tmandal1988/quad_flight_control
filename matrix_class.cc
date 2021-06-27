#include<iostream>
#include<string>

using namespace	std;


template <typename T>
class Matrix{
private:
	size_t nrows_, ncols_;
	T **matrix_;
public:
	Matrix();
	Matrix(size_t num_rows, size_t num_cols);
	Matrix(size_t num_rows, size_t num_cols, string type);
	~Matrix();
	// copy constructor
	Matrix(const Matrix<T>& matrix_to_copy);
	// member functions
	void PrintMatrix();
	Matrix Transpose();
	Matrix Inverse();
	// operator overloading
	Matrix operator= (const Matrix<T>& matrix_to_copy);
	// assign elements at row idx and col idx
	T& operator() (size_t r_idx, size_t c_idx);
	// read elements at row idx and col idx
	T operator() (size_t r_idx, size_t c_idx) const;
	// multiply two matrices
	Matrix operator* (const Matrix& matrix_to_multipy);
	// add two matrices
	Matrix operator+ (const Matrix& matrix_to_add);
	// subtract two matrices
	Matrix operator- (const Matrix& matrix_to_subtract);
	// add a scalar
	Matrix operator+ (const T x);
	// subtract a scalar
	Matrix operator- (const T x);
	// multiply a scalar
	Matrix operator* (const T x);
	// divide by a scalar
	Matrix operator/ (const T x);

};

template <typename T>
Matrix<T>::Matrix(){
	nrows_ = 0;
	ncols_ = 0;
	matrix_ = nullptr;
}

template <typename T>
Matrix<T>::Matrix(size_t num_rows, size_t num_cols){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_ = new T*[nrows_];
	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize matrix with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = 0;
		}
	}

}

template <typename T>
Matrix<T>::Matrix(size_t num_rows, size_t num_cols, string type){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_ = new T*[nrows_];
	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}

	if (type == "eye"){
		if (num_rows != num_cols){
			throw invalid_argument("To create identity matrix number of rows and columns should be equal\n");
		}
			// initialize matrix with zeros and then assign 1 to diagonal
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
				if (idx_r == idx_c)
					matrix_[idx_r][idx_c] = 1;
				else
					matrix_[idx_r][idx_c] = 0;
			}
		}
	}

}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix_to_copy){
	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_ = new T*[nrows_];
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize matrix with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
}

template <typename T>
Matrix<T>::~Matrix(){
	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			delete[] matrix_[idx_r];
		}
	}

	if(nrows_ > 0){
		delete[] matrix_;
	}
}

// member functions
template <typename T>
void Matrix<T>::PrintMatrix(){
	printf("*****************************\n");
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			printf("%g\t", matrix_[idx_r][idx_c]);
		}
		printf("\n");
	}
	printf("*****************************\n");
}

template <typename T>
Matrix<T> Matrix<T>::Transpose(){
	Matrix<T> transpose_matrix(ncols_, nrows_);
	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			transpose_matrix(idx_c, idx_r) = matrix_[idx_r][idx_c];
		}
	}
	return transpose_matrix;
}

template <typename T>
Matrix<T> Matrix<T>::Inverse(){
	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");

	T orig_matrix[ncols_][nrows_];
	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			orig_matrix[idx_r][idx_c] = matrix_[idx_r][idx_c];
	}

	Matrix<T> inverse_matrix(nrows_, ncols_, "eye");

	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		T diag_val = orig_matrix[idx_r][idx_r];
		if ( diag_val == 0){
			for( size_t idx_non_zero = idx_r + 1; idx_non_zero < nrows_; idx_r ++){
				if ( orig_matrix[idx_non_zero][idx_r] != 0){
					// swap rows
					for (size_t idx_c = 0; idx_c < ncols_; ++idx_c)
					{
						T temp_val_orig = orig_matrix[idx_r][idx_c];
						orig_matrix[idx_r][idx_c] = orig_matrix[idx_non_zero][idx_c];
						orig_matrix[idx_non_zero][idx_c] = temp_val_orig;

						temp_val_orig = inverse_matrix(idx_r, idx_c);
						inverse_matrix(idx_r, idx_c) = inverse_matrix(idx_non_zero, idx_c);
						inverse_matrix(idx_non_zero, idx_c) = temp_val_orig;
					}
					diag_val = orig_matrix[idx_r][idx_r];
					break;
				}
			}
		}

		// divide all the columns for idx_r row by diag_val
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			orig_matrix[idx_r][idx_c] = orig_matrix[idx_r][idx_c]/diag_val;
			inverse_matrix(idx_r, idx_c) = inverse_matrix(idx_r, idx_c)/diag_val;
		}
		for(size_t idx2_r = 0; idx2_r < nrows_; idx2_r++){
			if (idx2_r != idx_r){
				T value_to_multiply = orig_matrix[idx2_r][idx_r];
				for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
					orig_matrix[idx2_r][idx_c] = orig_matrix[idx2_r][idx_c] - orig_matrix[idx_r][idx_c]*value_to_multiply;
					inverse_matrix(idx2_r, idx_c) = inverse_matrix(idx2_r, idx_c) - inverse_matrix(idx_r, idx_c)*value_to_multiply;

				}
			}
		}
	}

	return inverse_matrix;
}


// operator overloading
template <typename T>
Matrix<T> Matrix<T>::operator= (const Matrix<T>& matrix_to_copy){
	if( &matrix_to_copy == this){
		return *this;
	}

	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			delete[] matrix_[idx_r];
		}
	}

	if(nrows_ > 0){
		delete[] matrix_;
	}

	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_ = new T*[nrows_];
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize matrix with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
	return *this;
}

template <typename T>
T& Matrix<T>::operator()(size_t r_idx, size_t c_idx){
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
T Matrix<T>::operator() (size_t r_idx, size_t c_idx) const {
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
Matrix<T> Matrix<T>::operator* (const Matrix& matrix_to_multipy){
	if (ncols_ != matrix_to_multipy.nrows_)
		throw invalid_argument("Matrix dimensions invalid for multiplication");
	Matrix<T> return_matrix(nrows_, matrix_to_multipy.ncols_);
	for (size_t idx_c = 0; idx_c < ncols_; idx_c++) {
    	for (size_t idx_r = 0; idx_r < nrows_; idx_r++) {
        	for (size_t idx2_c = 0; idx2_c < matrix_to_multipy.ncols_; idx2_c++) {
            	return_matrix.matrix_[idx_r][idx2_c] += matrix_[idx_r][idx_c] * matrix_to_multipy.matrix_[idx_c][idx2_c];
        	}
    	}
	}

	return return_matrix;
}

template <typename T>
Matrix<T> Matrix<T>::operator+ (const Matrix& matrix_to_add){
	// check dimensions
	if (nrows_ != matrix_to_add.nrows_ || ncols_ != matrix_to_add.ncols_)
		throw invalid_argument("Matrix dimensions invalid for addition");

	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] + matrix_to_add.matrix_[idx_r][idx_c];
	}

	return return_matrix;

}

template <typename T>
Matrix<T> Matrix<T>::operator- (const Matrix& matrix_to_subtract){
	// check dimensions
	if (nrows_ != matrix_to_subtract.nrows_ || ncols_ != matrix_to_subtract.ncols_)
		throw invalid_argument("Matrix dimensions invalid for subtraction");

	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] - matrix_to_subtract.matrix_[idx_r][idx_c];
	}

	return return_matrix;

}

template <typename T>
Matrix<T> Matrix<T>::operator+ (const T x){
	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] + x;
	}

	return return_matrix;
}

template <typename T>
Matrix<T> Matrix<T>::operator- (const T x){
	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] - x;
	}

	return return_matrix;
}

template <typename T>
Matrix<T> Matrix<T>::operator* (const T x){
	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] * x;
	}

	return return_matrix;
}

template <typename T>
Matrix<T> Matrix<T>::operator/ (const T x){
	Matrix<T> return_matrix(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_matrix(idx_r, idx_c) = matrix_[idx_r][idx_c] / x;
	}

	return return_matrix;
}

int main()
{
	Matrix<float> m(3, 3);
	m(0, 0) = 1;
	m(0, 1) = 5;
	m(0, 2) = 3;
	m(1, 0) = 10;
	m(1, 1) =  7;
	m(1, 2) = 1;
	m(2, 0) = 2;
	m(2, 1) = 0;
	m(2, 2) = 13;
	m.PrintMatrix();
	Matrix<float> m2 = m;
	m2.PrintMatrix();
	Matrix<float> m3;
	m3 = m2;
	m3.PrintMatrix();
	m3(1, 1) = 10;
	m3.PrintMatrix();
	printf("%g\n", m3(1,1));
	Matrix<float> m4;
	m4 = m3.Transpose();
	m4.PrintMatrix();
	Matrix<float> m5(3, 3, "eye");
	m5(0, 0) = 1;
	m5(0, 1) = 0;
	m5(0, 2) = 5;

	m5(1, 0) = 2;
	m5(1, 1) = 10;
	m5(1, 2) = 1;

	m5(2, 0) = -1;
	m5(2, 1) = 3;
	m5(2, 2) = -2;
	m5.PrintMatrix();
	Matrix<float> m6 = m5.Inverse();
	m6.PrintMatrix();
	m5(0, 0) = 0;
	m5(1, 1) = 0;
	m5(2, 2) = 0;
	m5.Inverse().PrintMatrix();
	Matrix<float> m7(3, 3);
	m7 = m5.Inverse()*m5;
	m7.PrintMatrix();
	m7 = m7 + m7;
	m7.PrintMatrix();
	m7 = m7 - m7;
	m7.PrintMatrix();
	m7 = m7 + 1;
	m7.PrintMatrix();
	m7 = m7 - 2;
	m7.PrintMatrix();
	m7 = m7/2;
	m7.PrintMatrix();
	m7 = m7*5;
	m7.PrintMatrix();
}