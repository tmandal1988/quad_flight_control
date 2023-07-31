#include "pwm_output_utils.h"

//Constructor
PwmOutputHelper::PwmOutputHelper(int num_channels):
								num_channels_(num_channels){
	rc_output_ptr_ = unique_ptr <RCOutput>{ new RCOutput_Navio2() };
	num_channels_  = min(num_channels_, MAX_NUM_PWM_OUTPUTS);
}

//Destructor
PwmOutputHelper::~PwmOutputHelper(){
}

// Initializer
void PwmOutputHelper::InitializePwmOutput(float frequency){
	// Initialize the pwm channels from 0 to num_channels - 1
	for(int idx = 0; idx < num_channels_; idx++){
		if( !(rc_output_ptr_->initialize(idx)) ) {
            cerr << "Unable to initialize PWM channel: "<<idx<<endl;
        }
        rc_output_ptr_->set_frequency(idx, frequency);	

        if ( !(rc_output_ptr_->enable(idx)) ) {
	    	cerr << "Unable to enable PWM channel: "<<idx<<endl;
		}	
	}
}

// Set PWM Duty Cycle
void PwmOutputHelper::SetPwmDutyCyle(vector<float> duty_cycle_ms){
	int num_channels_to_write = min((int)duty_cycle_ms.size(), num_channels_);
	for(size_t idx = 0; idx < num_channels_to_write; idx++){
		rc_output_ptr_->set_duty_cycle(idx, duty_cycle_ms[idx]);
	}
}