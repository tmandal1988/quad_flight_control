#include "matrix_inv_class.h"

template <typename T>
MatrixInv<T> MatrixInv<T>::Inverse(){
	size_t nrows_ = this->nrows_;
	size_t ncols_ = this->ncols_;
	vector< vector<T> > matrix_ = this->matrix_;
	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");

	T orig_matrix[ncols_][nrows_];
	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			orig_matrix[idx_r][idx_c] = matrix_[idx_r][idx_c];
	}

	MatrixInv<T> inverse_matrix(nrows_, ncols_, "eye");

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

// Explicit template instantiation
template class MatrixInv<float>;
template class MatrixInv<double>;