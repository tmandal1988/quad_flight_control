#ifndef UBXDRIVER_H
#define UBXDRIVER_H

#include <string>
#include <vector>
#include <iostream>

#include "SPIdev.h"

#define UBX_SYNC1 0xb5
#define UBX_SYNC2 0x62

#define PREAMBLE_OFFSET 2
#define CHECKSUM_OFFSET 2

#define RAD2DEG			57.2957795130823

#define PACKED __attribute__((__packed__))

static const int UBX_BUFFER_LENGTH = 1024;

class UbxDriver{
	public:
		UbxDriver(std::string name = "/dev/spidev0.0");
		~UbxDriver(){
		};
		int TestConnection();

		struct PACKED NavPvtData{
    		std::uint32_t iTow_ms;
    		std::uint16_t year;
    		std::uint8_t month;
    		std::uint8_t day;
    		std::uint8_t hour;
    		std::uint8_t min;
    		std::uint8_t sec;
    		std::uint8_t valid;
    		std::uint32_t tAcc_ns;
    		std::int32_t sFrac_ns;
    		std::uint8_t fixType;
    		std::uint8_t flags;
    		std::uint8_t flags2;
    		std::uint8_t numSV;
    		double lon_rad;
    		double lat_rad;
    		float hElpsd_m;
    		float hMsl_m;
    		float hAcc_m;
    		float vAcc_m;
    		float velN_mps;
    		float velE_mps;
    		float velD_mps;
    		float gSpeed_mps;
    		float headMot_rad;
    		float sAcc_mps;
    		float headAcc_rad;
    		float pDop;
    		std::uint16_t flags3;
    		float headVeh_rad;
    		float magDec_rad;
    		float magAcc_rad;

    	};

    	int GetNavPvtData(NavPvtData &_nav_pvt_data);

    	int ConfigureSolutionRate(std::uint16_t meas_rate_ms,
                                      std::uint16_t nav_rate = 1,
                                      std::uint16_t time_ref = 0);

	private:
		std::string spi_device_name_;

		struct PACKED UbxHeader {
	        std::uint8_t sync1;
	        std::uint8_t sync2;
	        std::uint8_t msg_class;
	        std::uint8_t msg_id;
	        std::uint16_t length;
    	};

    	struct PACKED CheckSum {
	        std::uint8_t CK_A;
	        std::uint8_t CK_B;
    	};

		// Class and ID enums
		 enum UbxClassAndIdx{
	        CLASS_CFG = 0x06,
	        MSG_CFG_RATE = 0x01,
	        RST_CFG = 0x09,
	        SOFT_RST_CFG = 0x04,
	        NAV5_CFG = 0x24,
	        PRT_CFG = 0x00,
	        HNR_RATE = 0x5C,
	        GNSS_CFG = 0x3E,

	        CLASS_NAV = 0x01,
	        MSG_HP_LLH = 0x14,
	        MSG_LLH = 0x02,
	        MSG_SOL_STATUS = 0x03,
	        MSG_VEL_NED = 0x12,
	        MSG_NAV_PVT = 0x07,

	        CLASS_HNR_NAV = 0x28,
	        MSG_NAV_HNR_PVT = 0x00,

	        CFG_SAMPLE_RATE = 0x08
    	};

		// Data scanning state
		enum State
	    {
	        SYNC1,
	        SYNC2,
	        CLASS,
	        ID,
	        LEN1,
	        LEN2,
	        PAYLOAD,
	        CK_A,
	        CK_B,
	        DONE,
	        INVALID
	    };

	    State state_{SYNC1};

		std::uint8_t ReadSingleByte();
		void UpdateUbxBuffer(std::uint8_t data);	 
		void ResetBuffParser();
		int SpliceMemory(std::uint8_t *dest, const void * const src, std::size_t size, int dest_offset = 0);
		int SendUbxMsg(std::uint8_t msg_class, std::uint8_t msg_id, void *msg, std::uint16_t size);
		std::uint16_t DecodeSingleGenericMessage(std::vector<std::uint8_t>& data);
		UbxDriver::CheckSum CalculateCheckSum(std::uint8_t *msg_buff, std::size_t size);
		int GetUbxAck(std::uint8_t msg_class, std::uint8_t msg_id);



		int ResetConfig();
		int ConfigureUbloxSpiPort();
		int ConfigureNavEngine();
		int SaveConfig();

		int EnableNavPvt();

		std::uint32_t GetUnsigned32BitData(const std::vector<std::uint8_t> &data, size_t start_idx);
		std::int32_t Get32BitData(const std::vector<std::uint8_t> &data, size_t start_idx);
		std::uint16_t GetUnsigned16BitData(const std::vector<std::uint8_t> &data, size_t start_idx);
		std::int16_t Get16BitData(const std::vector<std::uint8_t> &data, size_t start_idx);

		std::uint8_t message_buff_[UBX_BUFFER_LENGTH];   // Buffer for UBX message 
		unsigned int buff_idx_{0};
		unsigned int payload_length_{0};
		unsigned int message_length_{0};


		// UBX Message Structs
		// Reset the device
		struct PACKED ResetUblox {
	        std::uint16_t nav_bbr_mask;
	        std::uint8_t reset_mode;
	        std::uint8_t reserved;
    	};

    	// Reset the configuration
		struct PACKED ResetCfgUblox {
	        std::uint32_t clear_mask;
	        std::uint32_t save_mask;
	        std::uint32_t load_mask;
    	};	

    	// Configure the SPI port
    	struct PACKED CfgPrt{
	        std::uint8_t port_id;
	        std::uint8_t reserved1;
	        std::uint16_t tx_ready;
	        std::uint32_t spi_mode;
	        std::uint8_t reserved2;
	        std::uint8_t reserved3;
	        std::uint8_t reserved4;
	        std::uint8_t reserved5;
	        std::uint16_t in_proto_mask;
	        std::uint16_t out_proto_mask;
	        std::uint16_t flags;
	        std::uint8_t reserved6;
	        std::uint8_t reserved7;
    	};

    	struct PACKED CfgNavEng {
	        std::uint16_t set_mask;
	        std::uint8_t dyn_model;
	        std::uint8_t fix_mode;
	        std::int32_t fixed_alt;
	        std::uint32_t fixed_alt_var;
	        std::int8_t min_elev;
	        std::uint8_t dr_limit;
	        std::uint16_t p_dop;
	        std::uint16_t t_dop;
	        std::uint16_t p_acc;
	        std::uint16_t t_acc;
	        std::uint8_t static_hold_threshold;
	        std::uint8_t dgnss_timeout;
	        std::uint8_t cno_thresh_num_svs;
	        std::uint8_t cno_thresh;
	        std::uint8_t reserved1;
	        std::uint8_t reserved2;
	        std::uint16_t static_hold_max_dist;
	        std::uint8_t utc_standard;
	        std::uint8_t reserved3;
	        std::uint8_t reserved4;
	        std::uint8_t reserved5;
	        std::uint8_t reserved6;
	        std::uint8_t reserved7;
    	};

    	struct PACKED CfgMeasrate{
	        std::uint8_t msg_class;
	        std::uint8_t msg_id;
	        std::uint8_t msg_rate;
    	};

    	struct PACKED CfgNavRate {
	        std::uint16_t measure_rate;
	        std::uint16_t nav_rate;
	        std::uint16_t timeref;
    	};
};


#endif

