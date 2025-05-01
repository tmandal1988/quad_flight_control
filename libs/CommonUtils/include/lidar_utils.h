#ifndef LIDARUTILS_H
#define LIDARUTILS_H

#include <termios.h>
#include <unistd.h> // For usleep
#include <mutex>
#include <thread>
#include <pthread.h>
#include <errno.h>
#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <atomic>
#include <future>

#include "lw20api.h"

using namespace std;

class LidarHelper{
	public:
		// Constructors
		LidarHelper();

		// Destructor
		~LidarHelper();

		// To initialize the Lidar
		bool InitializeLidar();

		// Useful function that can be used to read Lidar data in a loop at the configured rate
		// This function can be passed to a thread to update range data
		void LidarReadLoop();

		// Start the Lidar running thread
		void CreateLidarThread();

		// Returns lidar range and validity
		bool GetLidarRange(float& lidar_range);

		// Stops the above loop
		void StopLidarReadLoop(){
			// Stop the Lidar reading loop
			stop_lidar_read_loop_.store(true);
		}

	private:
		struct lwSerialPort
		{
			int fd;
			bool connected;
		};

		struct lwSensorContext
		{
			lwLW20			lw20;
			lwSerialPort 	serialPort;
			uint8_t			inputBuffer[128];
			int32_t			inputBufferSize;
		};

		// Serial Port Name
		static string port_name_;

		// Create a thread
		thread lidar_thread_;
		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;
		// Create an LW20 context object to read Lidar data
		lwSensorContext lidar_context_ = {};

		// Create an LW20 Lidar service context
		lwServiceContext lidar_service_context_ = {};

		// Lidar range
		float lidar_range_ = 0.0f;
	    // Flag to indicate if lidar data has been updated
		mutable bool lidar_updated_ = false;

	    // Mutex to guard resource access between threads while running LidarReadLoop()
	    mutable mutex lidar_mutex_;

	    // Variable to indicate LidarReadLoop() to stop
	    atomic<bool> stop_lidar_read_loop_;

		//-------------------------------------------------------------------------
		// Platform Specific Functions.
		//-------------------------------------------------------------------------
	    inline int64_t PlatformGetMicrosecond();
		inline int32_t PlatformGetMS();

		// Function to connect/disconnect to the serial port
		static bool SerialDisconnect(lwSerialPort* com_port);	    
	    static bool SerialConnect(lwSerialPort* com_port, int bit_rate);
	    static int SerialWrite(lwSerialPort* com_port, char *buffer, int32_t buffer_size);
	    static int32_t SerialRead(lwSerialPort* com_port, char *buffer, int32_t buffer_size);

	    // Function to send a serial packet to the lidar
	    static bool SendPacket(lwLW20* Lw20, lwCmdPacket* Packet);
	    // Function to get a serial packet from the lidar
	    static bool GetPacket(lwLW20* Lw20, lwResponsePacket* Packet);
	    // Function to wait between measurements
	    static bool Sleep(lwLW20* Lw20, int32_t TimeMS);
	    // Function to stream data from the Lidar
	    static bool StreamResponse(lwLW20* Lw20, lwResponsePacket* Packet);
};

#endif