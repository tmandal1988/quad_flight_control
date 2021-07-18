#ifndef MATRIXFACT_H
#define MATRIXFACT_H
#endif

#include "matrix_base_class.h"
#include<array>

template <typename T>
class MatrixFact : public MatrixBase<T>{
	using MatrixBase<T>::MatrixBase;

	public:
		array<MatrixFact, 2> Lu();

};