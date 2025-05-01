#include "lidar_utils.h"

#define LW20_API_IMPLEMENTATION
#include "lw20api.h"

string LidarHelper::port_name_ = "/dev/ttyAMA0";

LidarHelper::LidarHelper(){
    // Variable to indicate LidarReadLoop() to stop
	stop_lidar_read_loop_.store(false);

	lidar_context_.lw20 = lw20CreateLW20();
	lidar_context_.lw20.userData = &lidar_context_;
}

//-------------------------------------------------------------------------
// Platform Specific Functions.
//-------------------------------------------------------------------------
inline int64_t LidarHelper::PlatformGetMicrosecond()
{
	timespec time;
	clock_gettime(CLOCK_REALTIME, &time);

	return time.tv_sec * 1000000 + time.tv_nsec / 1000;
}

inline int32_t LidarHelper::PlatformGetMS()
{
	return (PlatformGetMicrosecond() / 1000);
}

//-------------------------------------------------------------------------
// Com Port Implementation.
//-------------------------------------------------------------------------
bool LidarHelper::SerialDisconnect(lwSerialPort* com_port)
{
	if (com_port != 0 && com_port->fd >= 0)
	{
		close(com_port->fd);
	}

	com_port->fd = -1;
	com_port->connected = false;

	return true;
}

bool LidarHelper::SerialConnect(lwSerialPort* com_port, int bit_rate)
{
	int fd = -1;
	printf("Attempting connection on: %s\n", port_name_.c_str());
		
	fd = open(port_name_.c_str(), O_RDWR | O_NOCTTY | O_SYNC);
	
	if (fd < 0)
	{
		printf("Couldn't open serial port!\n");
		return false;
	}

	struct termios tty;
	memset(&tty, 0, sizeof(tty));
	if (tcgetattr(fd, &tty) != 0)
	{
		printf("Error from tcgetattr\n");
		return false;
	}

	cfsetospeed(&tty, bit_rate);
	cfsetispeed(&tty, bit_rate);

	tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8;
	tty.c_cflag |= (CLOCAL | CREAD);
	tty.c_cflag &= ~(PARENB | PARODD);
	tty.c_cflag |= 0;
	tty.c_cflag &= ~CSTOPB;
	tty.c_cflag &= ~CRTSCTS;
	tty.c_iflag &= ~IGNBRK;
	tty.c_iflag &= ~ICRNL;
	tty.c_iflag &= ~(IXON | IXOFF | IXANY);
	tty.c_lflag = 0;
	tty.c_oflag = 0;
	tty.c_cc[VMIN] = 0;
	tty.c_cc[VTIME] = 1;

	if (tcsetattr(fd, TCSANOW, &tty) != 0)
	{
		printf("Error from tcsetattr\n");
		return false;
	}

	com_port->fd = fd;
	com_port->connected = true;

	printf("%s opened for connection\n", port_name_.c_str());

	return true;
}

int LidarHelper::SerialWrite(lwSerialPort* com_port, char *buffer, int32_t buffer_size)
{
	if (!com_port)
	{
		printf("Can't write to null coms\n");
		return -1;
	}

	if (!com_port->connected)
	{
		printf("Can't write to non connected coms\n");
		return -1;
	}

	int writtenBytes = write(com_port->fd, buffer, buffer_size);

	if (writtenBytes != buffer_size)
	{
		printf("Could not send all bytes!\n");
		return -1;
	}

	return writtenBytes;
}

int32_t LidarHelper::SerialRead(lwSerialPort* com_port, char *buffer, int32_t buffer_size)
{
	if (!com_port)
	{
		printf("Can't read from null coms\n");
		return -1;
	}

	if (!com_port->connected)
	{
		printf("Can't read from non connected coms\n");
		return -1;
	}

	errno = 0;
	int readBytes = read(com_port->fd, buffer, buffer_size);

	//if (readBytes == 0)
		//printf("No Data: %d Error %d (%s)\n", readBytes, errno, strerror(errno));

	return readBytes;
}

bool LidarHelper::InitializeLidar(){
	printf("*********************LIDAR INITIALIZATION********************\n");
	bool lidar_com_status = SerialConnect(&lidar_context_.serialPort, B115200);

	if(lidar_com_status){
		lidar_service_context_.sendPacketCallback = SendPacket;
		lidar_service_context_.getPacketCallback = GetPacket;
		lidar_service_context_.sleepCallback = Sleep;
		lidar_service_context_.streamCallback = StreamResponse;

		// NOTE: Run event loop for first time init.
        bool event_status = runEventLoop(&lidar_context_.lw20, &lidar_service_context_);
        if(!event_status){
        	printf("*********************LIDAR INITIALIZATION FAILED*************\n");
        	return false;
        }

        lwProductInfo productInfo = executeCmd_GetProduct(&lidar_context_.lw20, &lidar_service_context_);
    	std::cout << "Product: " << lidar_context_.lw20.response.product.model
              	<< " Hardware: V" << lidar_context_.lw20.response.product.hardwareVersion
              	<< " Firmware V" << lidar_context_.lw20.response.product.firmwareVersion
              	<< "\n";

    	executeCmd_SetLaserMode(&lidar_context_.lw20, &lidar_service_context_, LWMS_48);
    	printf("*********************LIDAR INITIALIZATION********************\n");

	    return true;
	}

	return false;	
}

bool LidarHelper::SendPacket(lwLW20* Lw20, lwCmdPacket* Packet)
{
	//std::cout << "Send Packet " << Packet->length << "\n";
	lwSensorContext* lidar_context_ = (lwSensorContext*)Lw20->userData;
	if (SerialWrite(&lidar_context_->serialPort, (char*)Packet->buffer, Packet->length) != -1)
		return true;

	return false;
}

bool LidarHelper::GetPacket(lwLW20* Lw20, lwResponsePacket* Packet)
{
	int32_t timeout_ms = 5000;
	lwSensorContext* lidar_context_ = (lwSensorContext*)Lw20->userData;
	auto start_time = std::chrono::steady_clock::now();
	while (true)
	{
		if (lidar_context_->inputBufferSize == 0)
		{
			int bytesRead = 0;
			if ((bytesRead = SerialRead(&lidar_context_->serialPort, (char*)lidar_context_->inputBuffer, sizeof(lidar_context_->inputBuffer))) != -1)
				lidar_context_->inputBufferSize = bytesRead;
			else
				return false;
		}

		lwResolvePacketResult packetResolve = lw20ResolvePacket(&lidar_context_->lw20.response, lidar_context_->inputBuffer, lidar_context_->inputBufferSize);

		// NOTE: You can use a circular buffer or so to avoid the shuffle.
		if (packetResolve.bytesRead > 0)
		{
			int32_t remaining = lidar_context_->inputBufferSize - packetResolve.bytesRead;
			for (int i = 0; i < remaining; ++i)
				lidar_context_->inputBuffer[i] = lidar_context_->inputBuffer[packetResolve.bytesRead + i];

			lidar_context_->inputBufferSize = remaining;
		}
		
		if (packetResolve.status == LWRPS_COMPLETE)
		{
			return true;
		}

		// Check for timeout
        auto elapsed_time = std::chrono::steady_clock::now() - start_time;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count() >= timeout_ms) {
            break; // Exit the loop if timeout is reached
        }
	}

	return false;
}

bool LidarHelper::Sleep(lwLW20* Lw20, int32_t TimeMS)
{
	usleep(TimeMS * 1000);
	return true;
};

bool LidarHelper::StreamResponse(lwLW20* Lw20, lwResponsePacket* Packet)
{
	if (Packet->type == LWC_SERVO_SCAN)
	{
		std::cout << "Scan: " << Packet->scanSample.angle << " " << Packet->scanSample.firstPulse << " " << Packet->scanSample.lastPulse << "\n";
	}
	else if (Packet->type == LWC_SERVO_POSITION)
	{
		std::cout << "Pos: " << Packet->floatValue << "\n";
	}
	else if (Packet->type == LWC_LASER_TEMPERATURE)
	{
		std::cout << "Temp: " << Packet->floatValue << "\n";
	}

	return true;
};

void LidarHelper::LidarReadLoop(){
	while (1){
		lidar_range_ = executeCmd_GetLaserDistanceFirst(&lidar_context_.lw20, &lidar_service_context_);
		{
			unique_lock<mutex> lidar_data_lock(lidar_mutex_);
			lidar_updated_ = true;
		}
		if(stop_lidar_read_loop_.load()){
			break;
		}
	}
	SerialDisconnect(&lidar_context_.serialPort);
}// LidarReadLoop function

bool LidarHelper::GetLidarRange(float& lidar_range){
	if(lidar_updated_){
		{
			unique_lock<mutex> lidar_data_lock(lidar_mutex_);
			lidar_range = lidar_range_;
			lidar_updated_ = false;
		}
		return true;
	}else{
		return false;
	}
}

void LidarHelper::CreateLidarThread(){
	lidar_thread_ = thread(&LidarHelper::LidarReadLoop, this);

    CPU_ZERO(&cpuset_);
    CPU_SET(0, &cpuset_);

	pthread_getschedparam(lidar_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = 5;
	pthread_setschedparam(lidar_thread_.native_handle(), SCHED_FIFO, &sch_);
	int rc = pthread_setaffinity_np(lidar_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on GPS Thread: " << rc << "\n";
    }
}

LidarHelper::~LidarHelper(){
	// Making sure that GPS reading loop is stopped
	stop_lidar_read_loop_.store(true);
	if(lidar_thread_.joinable())
		lidar_thread_.join();
}