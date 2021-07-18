#include <Matrix/matrix_factorization_class.h>


int main(){
	MatrixFact<float> m(3, 3);
	m(0, 0) = 1;
	m(0, 1) = 5;
	m(0, 2) = 3;
	m(1, 0) = 10;
	m(1, 1) = 7;
	m(1, 2) = 1;
	m(2, 0) = 2;
	m(2, 1) = 0;
	m(2, 2) = 13;

	MatrixFact<float> m1 = m;
	m1.PrintMatrix();

	MatrixFact<float> m2;
	m2 = m1;
	array<MatrixFact<float>, 2> lu_decomp = m2.Lu();
	lu_decomp[0].PrintMatrix();
	lu_decomp[1].PrintMatrix();

	MatrixBase<float> m3;
	m3 = lu_decomp[0]*lu_decomp[1];
	m3.PrintMatrix();
	return 0;
}