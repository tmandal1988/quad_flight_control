#include <Matrix/matrix_inv_class.h>


int main(){
	MatrixInv<float> m(3, 3);
	m(0, 0) = 1;
	m(0, 1) = 5;
	m(0, 2) = 3;
	m(1, 0) = 10;
	m(1, 1) = 7;
	m(1, 2) = 1;
	m(2, 0) = 2;
	m(2, 1) = 0;
	m(2, 2) = 13;
	m.PrintMatrix();

	MatrixInv<float> m1 = m;
	m1.PrintMatrix();

	MatrixInv<float> m2;
	m2 = m1;
	m2.Inverse().PrintMatrix();

	return 0;
}