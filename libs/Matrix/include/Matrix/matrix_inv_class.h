#ifndef MATRIXINV_H
#define MATRIXINV_H

#include "matrix_base_class.h"

template <typename T>
class MatrixInv : public MatrixBase<T>{
	using MatrixBase<T>::MatrixBase;
	public:
		MatrixInv Inverse();
		MatrixInv(const MatrixBase<T>& matrix):MatrixBase<T>(matrix){}

};

#endif