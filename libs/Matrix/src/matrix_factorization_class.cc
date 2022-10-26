#include "matrix_factorization_class.h"

template<typename T>
array<MatrixFact<T>, 2> MatrixFact<T>::Lu(){
	size_t nrows_ = this->nrows_;
	size_t ncols_ = this->ncols_;
	vector< vector<T> > matrix_   = this->matrix_;
	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");
	array<MatrixFact<T>, 2> lu_decomp = {MatrixFact<T>(nrows_, ncols_, "eye"), MatrixFact<T>(nrows_, ncols_)};
	T temp_sum;
	for (size_t idx_c = 0; idx_c < ncols_; idx_c ++){
		for (size_t idx_t = idx_c; idx_t < ncols_; idx_t++){
			temp_sum = 0;
			for (size_t idx_k = 0; idx_k < idx_c; idx_k++)
				temp_sum += lu_decomp[0](idx_c, idx_k)*lu_decomp[1](idx_k, idx_t);

			lu_decomp[1](idx_c, idx_t) = matrix_[idx_c][idx_t] - temp_sum;
		}

		for (size_t idx_t = (idx_c + 1); idx_t < ncols_; idx_t++){
			temp_sum = 0;
			for (size_t idx_k = 0; idx_k < idx_c; idx_k++)
				temp_sum += lu_decomp[0](idx_t, idx_k) * lu_decomp[1](idx_k, idx_c);
			lu_decomp[0](idx_t, idx_c) = (matrix_[idx_t][idx_c] - temp_sum)/lu_decomp[1](idx_c, idx_c);
		}
	}
	return lu_decomp;
}

// Explicit template instantiation
template class MatrixFact<float>;
template class MatrixFact<double>;