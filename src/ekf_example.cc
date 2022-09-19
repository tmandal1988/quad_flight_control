#include <ekf_21dof_class.h>
#include <Matrix/matrix_factorization_class.h>


int main(){
/*	MatrixFact<float> m(3, 3);
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

	MatrixInv<float> m3;
	m3 = lu_decomp[0]*lu_decomp[1];
	m3.PrintMatrix();*/
	MatrixInv<float> initial_state(21, 1);
	MatrixInv<float> sensor_meas(9, 1);
	MatrixInv<float> state_sensor_val(6, 1);
	MatrixInv<float> process_noise_q(21, 21, "eye");
	MatrixInv<float> meas_noise_r(7, 7, "eye");

	Ekf21Dof<float> imu_gps_ekf(0.01, initial_state, process_noise_q, meas_noise_r);

	imu_gps_ekf.run(state_sensor_val, sensor_meas);


	return 0;
}