#define MAX_BUFF_SIZE 1024

// Create two flags to indicate if the above buffers are full or not
bool is_data_buff1_full = false;
bool is_data_buff2_full = false;

// Struct defining what data to save
struct data_fields{
	// Time it took for the EKF to run
	float dt_s;
	// GPS Time if available
	float time_of_week_ms;
	// Raw Sensor Data
	float accel_mps2[3];
	float gyro_radps[3];
	// Corrected mag data for alignment and offset error
	float mag_t[3];
	// GPS NED position and velocity if available;
	float ned_pos_m[3];
	float ned_vel_mps[3];
	// Current EKF State
	float ekf_current_state[15];

};

// Create two buffers to store the data to save before writing
data_fields data_to_save1[MAX_BUFF_SIZE];
data_fields data_to_save2[MAX_BUFF_SIZE];