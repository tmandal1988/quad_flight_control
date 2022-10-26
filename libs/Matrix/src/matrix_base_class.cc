#include "matrix_base_class.h"

template <typename T>
MatrixBase<T>::MatrixBase(){
	nrows_ = 0;
	ncols_ = 0;
}

template <typename T>
MatrixBase<T>::MatrixBase(initializer_list< initializer_list<T> > initial_data){
	// get the number of rows
	nrows_ = initial_data.size();
	matrix_.resize(nrows_);

	// get the number of columns
	const initializer_list<T> *initial_data_first_row = initial_data.begin();
	ncols_ = (*initial_data_first_row).size();
	size_t ncols_temp = ncols_;

	// iterate over the initialization list and fill out the matrix
	size_t idx_r = 0;
	for (initializer_list<T> row : initial_data) {
        ncols_temp = row.size();
        if (ncols_temp != ncols_)
        	throw invalid_argument("Initializer list dimensions are not consistent accross rows or columns");
        	matrix_[idx_r].resize(ncols_);
        	size_t idx_c = 0;
        	for(T data : row){
        		matrix_[idx_r][idx_c] =  data;
        		idx_c++;
        	}
        	idx_r++;
    }
}

template <typename T>
MatrixBase<T>::MatrixBase(size_t num_rows, size_t num_cols){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_.resize(nrows_);

	// initialize MatrixBase with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		// assign memory for the columns
		matrix_[idx_r].resize(ncols_, 0);
	}

}

template <typename T>
MatrixBase<T>::MatrixBase(size_t num_rows, size_t num_cols, string type){
	if (type == "eye"){
		if (num_rows != num_cols){
			throw invalid_argument("To create identity Matrix number of rows and columns should be equal\n");
		}
	}

	nrows_ = num_rows;
	ncols_ = num_cols;

	matrix_.resize(nrows_);
	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r].resize(ncols_, 0);
		matrix_[idx_r][idx_r] = 1;
	}

}

template <typename T>
MatrixBase<T>::MatrixBase(const MatrixBase<T>& matrix_to_copy){
	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_.resize(nrows_);
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		// assign memory for the columns
		matrix_[idx_r].resize(ncols_, 0);
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
}

template <typename T>
MatrixBase<T>::~MatrixBase(){
	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			matrix_[idx_r].clear();
		}
	}

	if(nrows_ > 0){
		matrix_.clear();
	}

}

// member functions
template <typename T>
void MatrixBase<T>::PrintMatrix() const{
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
void MatrixBase<T>::PrintRow(size_t r_idx) const{
	printf("*****************************\n");
	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
		printf("%g\t", matrix_[r_idx][idx_c]);
	}
	printf("*****************************\n");
}

template <typename T>
void MatrixBase<T>::PrintCol(size_t c_idx) const{
	printf("*****************************\n");
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		printf("%g\n", matrix_[idx_r][c_idx]);
	}
	printf("*****************************\n");
}

template <typename T>
MatrixBase<T> MatrixBase<T>::Transpose(){
	MatrixBase<T> transpose_matrix(ncols_, nrows_);
	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			transpose_matrix(idx_c, idx_r) = matrix_[idx_r][idx_c];
		}
	}
	return transpose_matrix;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::GetRow(size_t r_idx){
	MatrixBase<T> row_matrix(1, ncols_);
	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			row_matrix(0, idx_c) = matrix_[r_idx][idx_c];
	}
	return row_matrix;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::GetCol(size_t c_idx){
	MatrixBase<T> col_matrix(nrows_, 1);
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			col_matrix(idx_r, 0) = matrix_[idx_r][c_idx];
	}
	return col_matrix;
}

// operator overloading
template <typename T>
MatrixBase<T> MatrixBase<T>::operator= (const MatrixBase<T>& matrix_to_copy){
	if( &matrix_to_copy == this){
		return *this;
	}

	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			matrix_[idx_r].clear();
		}
	}

	if(nrows_ > 0){
		matrix_.clear();
	}

	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_.resize(nrows_);
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r].resize(ncols_);
	}
	// initialize Matrix with the matrix_to_copy values
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
	return *this;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator= (const initializer_list< initializer_list<T> > list_data){
// get the number of columns
	size_t nrows_in = list_data.size();
	const initializer_list<T> *list_data_first_row = list_data.begin();
	size_t ncols_in = (*list_data_first_row).size();
	size_t ncols_temp = ncols_;

	if( (nrows_ == nrows_in) && (ncols_ == ncols_in) ){
		// iterate over the initialization list and fill out the matrix
		size_t idx_r = 0;
		for (initializer_list<T> row : list_data) {
        	ncols_temp = row.size();
        	if (ncols_temp != ncols_)
        		throw invalid_argument("Initializer list dimensions are not consistent accross rows or columns1");
        		size_t idx_c = 0;
        		for(T data : row){
        			matrix_[idx_r][idx_c] = data;
        			idx_c++;
        		}
        	idx_r++;
    	}
	}else{
		if (ncols_ > 0){
			for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
				matrix_[idx_r].clear();
			}
		}

		if(nrows_ > 0){
			matrix_.clear();
		}

		// get the number of rows
		nrows_ = nrows_in;
		matrix_.resize(nrows_);
		ncols_ = ncols_in;
	

		// iterate over the initialization list and fill out the matrix
		size_t idx_r = 0;
		for (initializer_list<T> row : list_data) {
        	ncols_temp = row.size();
        	if (ncols_temp != ncols_)
        		throw invalid_argument("Initializer list dimensions are not consistent accross rows or columns2");
        		matrix_[idx_r].resize(ncols_);
        		size_t idx_c = 0;
        		for(T data : row){
        			matrix_[idx_r][idx_c] = data;
        			idx_c++;
        		}
        	idx_r++;
    	}
	}

    return *this;

}

template <typename T>
T& MatrixBase<T>::operator()(size_t r_idx, size_t c_idx){
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
T& MatrixBase<T>::operator()(size_t r_idx){
	if (r_idx >= nrows_ )
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][0];
}

template <typename T>
T MatrixBase<T>::operator() (size_t r_idx, size_t c_idx) const {
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
T MatrixBase<T>::operator() (size_t r_idx) const {
	if (r_idx >= nrows_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][0];
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator* (const MatrixBase& matrix_to_multipy){
	if (ncols_ != matrix_to_multipy.nrows_)
		throw invalid_argument("MatrixBase dimensions invalid for multiplication");
	MatrixBase<T> return_MatrixBase(nrows_, matrix_to_multipy.ncols_);
	for (size_t idx_c = 0; idx_c < ncols_; idx_c++) {
    	for (size_t idx_r = 0; idx_r < nrows_; idx_r++) {
        	for (size_t idx2_c = 0; idx2_c < matrix_to_multipy.ncols_; idx2_c++) {
            	return_MatrixBase.matrix_[idx_r][idx2_c] += matrix_[idx_r][idx_c] * matrix_to_multipy.matrix_[idx_c][idx2_c];
        	}
    	}
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator+ (const MatrixBase& matrix_to_add){
	// check dimensions
	if (nrows_ != matrix_to_add.nrows_ || ncols_ != matrix_to_add.ncols_)
		throw invalid_argument("MatrixBase dimensions invalid for addition");

	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] + matrix_to_add.matrix_[idx_r][idx_c];
	}

	return return_MatrixBase;

}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator- (const MatrixBase& matrix_to_subtract){
	// check dimensions
	if (nrows_ != matrix_to_subtract.nrows_ || ncols_ != matrix_to_subtract.ncols_)
		throw invalid_argument("MatrixBase dimensions invalid for subtraction");

	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] - matrix_to_subtract.matrix_[idx_r][idx_c];
	}

	return return_MatrixBase;

}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator+ (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] + x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator- (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] - x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator* (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] * x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator/ (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] / x;
	}

	return return_MatrixBase;
}

// getters and setters
template<typename T>
size_t MatrixBase<T>::get_ncols(){
	return ncols_;
}

template<typename T>
size_t MatrixBase<T>::get_nrows(){
	return nrows_;
}

template<typename T>
T MatrixBase<T>::get_element(size_t idx_r, size_t idx_c){
	return matrix_[idx_r][idx_c];
}

template<typename T>
void MatrixBase<T>::set_element(size_t idx_r, size_t idx_c, T val){
	matrix_[idx_r][idx_c] = val;
}

template<typename T>
bool MatrixBase<T>::is_empty(){
	if (matrix_.empty()){
		return true;
	}else{
		return false;
	}
}

// Explicit template instantiation
template class MatrixBase<float>;
template class MatrixBase<double>;