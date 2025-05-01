#include "external_mag.h"

#include <signal.h>
#include <chrono>

// To catch SIGINT
volatile sig_atomic_t sigint_flag = 0;

void sigint_handler(int sig){ // can be called asynstd::chronously
  sigint_flag = 1; // set flag
}

int main(int argc, char *argv[]){
	ExtMagDriver ext_mag;
	int status = ext_mag.TestConnection();
	if(status < 0){
		printf("Unable to initialize external mag\n");
	}

	int16_t mag_val[3];

	while(1){
		ext_mag.TestConnection();
		ext_mag.StartMeasurement();

		usleep(6500);

		printf("X_mag: %d, Y mag: %d, Z_mag: %d\n", mag_val[0], mag_val[1], mag_val[2]);

		ext_mag.GetMeasurement(mag_val);

		if(sigint_flag == 1)
   			break;
	}

	return 0;
}