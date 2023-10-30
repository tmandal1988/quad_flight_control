#ifndef PWMOUTPUTUTILS_H
#define PWMOUTPUTUTILS_H

#include "Navio/Navio2/PWM.h"
#include "Navio/Navio2/RCOutput_Navio2.h"
#include "Navio/Common/Util.h"
#include "Navio/Common/RCOutput.h"

#include <vector>
#include <unistd.h>
#include <memory>
#include <iostream>

using namespace std;

class PwmOutputHelper{

	#define MAX_NUM_PWM_OUTPUTS 14

	public:
		// Constructors
		PwmOutputHelper(int num_channels);

		// Initializer
		void InitializePwmOutput(float frequency = 400.0);

		// Set PWM duty cycle
		void SetPwmDutyCyle(vector<float> duty_cycle_ms);

		// Destructor
		~PwmOutputHelper();

	private:
		unique_ptr<RCOutput> rc_output_ptr_;
		int num_channels_;

};

#endif