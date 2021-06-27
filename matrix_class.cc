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
	Matrix(const Matrix& matrix_to_copy);
	// member functions
	void PrintMatrix();
	Matrix Transpose();
	Matrix Inverse();
	// operator overloading
	Matrix operator= (const Matrix& matrix_to_copy);
	T& operator() (size_t r_idx, size_t c_idx);
	T operator() (size_t r_idx, size_t c_idx) const;

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
			printf("here\n");
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
Matrix<T>::Matrix(const Matrix& matrix_to_copy){
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
					printf("%lu\n", idx_non_zero);
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
Matrix<T> Matrix<T>::operator= (const Matrix& matrix_to_copy){
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
}