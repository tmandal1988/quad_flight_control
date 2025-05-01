#include "UbloxDriver.h"

UbxDriver::UbxDriver(std::string name):
spi_device_name_(name){

}

uint8_t UbxDriver::ReadSingleByte(){
	uint8_t to_gps_data = 0x00, from_gps_data = 0x00;
	// Per Ublox Interface document "Back-To-Back Read and Write Access"
    SPIdev::transfer(spi_device_name_.c_str(), &to_gps_data, &from_gps_data, 1, 2500000);
    return from_gps_data;
}

int UbxDriver::SendUbxMsg(uint8_t msg_class, uint8_t msg_id, void *msg, std::uint16_t size){
	uint8_t buffer[UBX_BUFFER_LENGTH];

    UbxHeader header;
    header.sync1 = UBX_SYNC1;
    header.sync2 = UBX_SYNC2;
    header.msg_class = msg_class;
    header.msg_id    = msg_id;
    header.length    = size;

    int offset = SpliceMemory(buffer, &header, sizeof(UbxHeader));
    offset = SpliceMemory(buffer, msg, size, offset);

    auto checksum = CalculateCheckSum(buffer, offset);
    offset = SpliceMemory(buffer, &checksum, sizeof(CheckSum), offset);

    // for(size_t idx = 0; idx < offset; idx++){
    //     printf("Msg Sent [%d]: 0x%02x\n", idx, buffer[idx]);
    // }

    return SPIdev::transfer(spi_device_name_.c_str(), buffer, nullptr, offset);
}

int UbxDriver::SpliceMemory(uint8_t *dest, const void * const src, std::size_t size, int dest_offset){
    std::memmove(dest + dest_offset, src, size);
    return dest_offset + size;
}

UbxDriver::CheckSum UbxDriver::CalculateCheckSum(uint8_t *msg_buff, std::size_t size) {
    CheckSum checksum;
    checksum.CK_A = checksum.CK_B = 0;

    for (size_t idx = PREAMBLE_OFFSET; idx < size; idx++) {
        checksum.CK_A += msg_buff[idx];
        checksum.CK_B += checksum.CK_A;
    }
    return checksum;
}

void UbxDriver::ResetBuffParser(){
	buff_idx_ = 0;
	state_ = SYNC1;
}

void UbxDriver::UpdateUbxBuffer(uint8_t data){
    if (state_ != DONE){
        message_buff_[buff_idx_] = data;
        buff_idx_++;
    }

    // If Buffer is already full reset and startover, any old data or partial
    // data will be discarded
    if(buff_idx_ == UBX_BUFFER_LENGTH)
    	 ResetBuffParser();

    switch (state_)
    {
	    case SYNC1:
	        if (data == UBX_SYNC1)
	            state_ = SYNC2;
	        else
	            ResetBuffParser();
	        break;

	    case SYNC2:
	        if (data == UBX_SYNC2)
	            state_ = CLASS;
	        else
	            ResetBuffParser();
	        break;

	    case CLASS:
	        state_ = ID;
	        break;

	    case ID:
	        state_ = LEN1;
	        break;

	    case LEN1:
	        payload_length_ = data;
	        state_ = LEN2;
	        break;

	    case LEN2:
	        payload_length_ += data << 8;
	        state_ = PAYLOAD;
	        break;

	    case PAYLOAD:
	        if (buff_idx_ == payload_length_ + 6)
	            state_ = CK_A;
	        break;

	    case CK_A:
	        state_ = CK_B;
	        break;

	    case CK_B:
	        message_length_ = 6 + payload_length_ + 2;
	        state_ = DONE;
	        break;

	    default:
	        break;
    }

}

std::uint16_t UbxDriver::DecodeSingleGenericMessage(std::vector<uint8_t>& data){
	uint8_t from_gps_data;
    size_t byte_rec_count = 0;

    size_t msg_start_idx;

    std::uint16_t msg_id = 0x00;

    // Clear the data
    data.clear();
    ResetBuffParser();

    bool data_frame_flag = true;

    while (byte_rec_count < UBX_BUFFER_LENGTH/2)
    {
        // From now on, we will send zeroes to the receiver, which it will ignore
        // However, we are simultaneously getting useful information from it
        from_gps_data = ReadSingleByte();
        // Scanner checks the message structure with every byte received
        UpdateUbxBuffer(from_gps_data);
        if (state_ == DONE)
        {
            // Once we have a full message we decode it and reset the scanner, making it look for another message
            // in the data stream, coming over SPI
            msg_start_idx = buff_idx_ - message_length_; // count the message start position

            // All UBX messages start with 2 sync chars: 0xb5 and 0x62
            if (message_buff_[msg_start_idx] != UBX_SYNC1)
            	data_frame_flag = false;
            if (message_buff_[msg_start_idx + 1]!= UBX_SYNC2)
            	data_frame_flag = false;

            // Count the checksum
            CheckSum checksum = CalculateCheckSum(&message_buff_[msg_start_idx], message_length_ - 2);

            if (checksum.CK_A != message_buff_[msg_start_idx + message_length_ - 2])
            	data_frame_flag = false;
            if (checksum.CK_B != message_buff_[msg_start_idx + message_length_ - 1])
            	data_frame_flag = false;

            // printf("Computed CheckSum:0x%02x, 0x%02x, Received CheckSum: 0x%02x, 0x%02x\n", checksum.CK_A, checksum.CK_B, message_buff_[msg_start_idx + message_length_ - 2], message_buff_[msg_start_idx + message_length_ - 1]);
            
            if(data_frame_flag){
                // If we got everything right, then return the raw data
                //data.clear();
                data.push_back(message_buff_[msg_start_idx]);
                data.push_back(message_buff_[msg_start_idx + 1]);
                data.push_back(message_buff_[msg_start_idx + 2]);
                data.push_back(message_buff_[msg_start_idx + 3]);
                data.push_back(message_buff_[msg_start_idx + 4]);
                data.push_back(message_buff_[msg_start_idx + 5]);

                for(size_t idx = 6; idx < message_length_; idx++)
                    data.push_back(message_buff_[msg_start_idx + idx]);

                msg_id = message_buff_[msg_start_idx + 2] << 8 | message_buff_[msg_start_idx + 3]; // ID is a two-byte number with little endianness

                // Reset data buffer
                ResetBuffParser();
                return msg_id;
            }else{
            	ResetBuffParser();
            	data_frame_flag = true;
            }      
        }

        byte_rec_count++;
    }

    return msg_id;
}

int UbxDriver::GetUbxAck(uint8_t msg_class, uint8_t msg_id){
	std::vector<uint8_t> raw_data;
    int ack_received = -1;
    size_t count = 0;

    if(DecodeSingleGenericMessage(raw_data) > 0)
    {
        if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == msg_class && raw_data[7] == msg_id){
            ack_received = 1;
        }
    }

    return ack_received;
}

/******************************************************************************************************************************************
 *********************************************************    RECEIVER CONFIG    **********************************************************
 ******************************************************************************************************************************************/

int UbxDriver::TestConnection(){

    // reset Ublox config
    if(ResetConfig() < 0){
        return -1;
    }

    uint8_t spi_config_count = 0;
    while(spi_config_count < 3){
        if(ConfigureUbloxSpiPort() < 0){
            printf("Retrying Ublox SPI config\n");
            usleep(2000);
            spi_config_count++;
        }else{
            break;
        }
    }

    if(spi_config_count > 3){
        printf("Unable to config Ublox SPI after 3 attempts, aborting\n");
        return -1;
    }


    if(ConfigureNavEngine() < 0){
        return -1;
    }

    if(ConfigureSolutionRate(50) < 0){
        return -1;
    }

    if(EnableNavPvt() < 0){
        return -1;
    }

    if(SaveConfig() < 0){
        return -1;
    }

    return 1;
}


int UbxDriver::ResetConfig(){
    ResetUblox msg_rst;
    msg_rst.nav_bbr_mask = 1; //Warm Start
    msg_rst.reset_mode = 1; // Controlled Software Reset
    msg_rst.reserved = 0;
    printf("Resetting Ublox Hardware....\n");
    int send_status = SendUbxMsg(CLASS_CFG, SOFT_RST_CFG, &msg_rst, sizeof(ResetUblox));

    if (send_status < 0)
        return -1;
    // Just wait and don't check for ACK from ublox as ACK after reset is not reliable
    usleep(500000);
    printf("Ublox Hardware Reset Complete\n");

    ResetCfgUblox msg_cf_rst;
    msg_cf_rst.clear_mask   = 0x0000001F;

    msg_cf_rst.save_mask    = 0x00000000;

    msg_cf_rst.load_mask    = 0x0000001F;

    printf("Loading default Ublox config....\n");
    send_status = SendUbxMsg(CLASS_CFG, RST_CFG, &msg_cf_rst, sizeof(ResetCfgUblox));
    if(GetUbxAck(CLASS_CFG, RST_CFG) < 0){
    	std::cerr << "Could not load default Ublox config over SPI\n";
        return -1;
    }
    printf("Default Ublox config loaded\n");

    // printf("Reading the loaded config...\n");
    // send_status = SendUbxMsg(CLASS_CFG, RST_CFG, nullptr, 0);
    // std::vector<std::uint8_t> config_data;
    // DecodeSingleGenericMessage(config_data);

    // for (const auto & element : config_data) {
    //     printf("config_data: 0x%02x\n", element);
    // }

    return 1;
}

int UbxDriver::SaveConfig(){
    ResetCfgUblox msg_cf_rst;
    msg_cf_rst.clear_mask   = 0x00000000;

    msg_cf_rst.save_mask    = 0x0000001F;

    msg_cf_rst.load_mask    = 0x00000000;

    printf("Saving modified Ublox config....\n");
    int send_status = SendUbxMsg(CLASS_CFG, RST_CFG, &msg_cf_rst, sizeof(ResetCfgUblox));
    if(GetUbxAck(CLASS_CFG, RST_CFG) < 0){
        std::cerr << "Could not save modified Ublox config over SPI\n";
        return -1;
    }
    printf("Modified Ublox config saved to non-volatile memory\n");

    // printf("Reading the loaded config...\n");
    // send_status = SendUbxMsg(CLASS_CFG, RST_CFG, nullptr, 0);
    // std::vector<std::uint8_t> config_data;
    // DecodeSingleGenericMessage(config_data);

    // for (const auto & element : config_data) {
    //     printf("config_data: 0x%02x\n", element);
    // }

    return 1;
}

int UbxDriver::ConfigureUbloxSpiPort(){
    CfgPrt msg_prt;
    msg_prt.port_id = 4;
    msg_prt.reserved1 = 0;
    msg_prt.tx_ready = 0;
    msg_prt.spi_mode = 0;
    msg_prt.reserved2 = 0;
    msg_prt.reserved3 = 0;
    msg_prt.reserved4 = 0;
    msg_prt.reserved5 = 0;
    msg_prt.in_proto_mask = 1;
    msg_prt.out_proto_mask = 1;
    msg_prt.flags = 0;
    msg_prt.reserved6 = 0;
    msg_prt.reserved7 = 0;

    printf("Configuring Ublox SPI port for UBX protocol only ....\n");
    int send_status = SendUbxMsg(CLASS_CFG, PRT_CFG, &msg_prt, sizeof(CfgPrt));

    if(GetUbxAck(CLASS_CFG, PRT_CFG) < 0){
        std::cerr << "Could not configure Ublox SPI port for UBX only\n";
        return -1;
    }
    printf("Successfully configured Ublox SPI port for UBX only\n");

    return 1;
}

 int UbxDriver::ConfigureNavEngine(){
    CfgNavEng msg_nav_eng;

    msg_nav_eng.set_mask = 5;
    msg_nav_eng.dyn_model = 7;
    msg_nav_eng.fix_mode = 2;;
    msg_nav_eng.fixed_alt = 0;
    msg_nav_eng.fixed_alt_var = 0;
    msg_nav_eng.min_elev = 0;
    msg_nav_eng.dr_limit = 0;
    msg_nav_eng.p_dop = 0;
    msg_nav_eng.t_dop = 0;
    msg_nav_eng.p_acc = 0;
    msg_nav_eng.t_acc = 0;
    msg_nav_eng.static_hold_threshold = 0;
    msg_nav_eng.dgnss_timeout = 0;
    msg_nav_eng.cno_thresh_num_svs = 0;
    msg_nav_eng.cno_thresh = 0;
    msg_nav_eng.reserved1 = 0;
    msg_nav_eng.reserved2 = 0;
    msg_nav_eng.static_hold_max_dist = 0;
    msg_nav_eng.utc_standard = 0;
    msg_nav_eng.reserved3 = 0;
    msg_nav_eng.reserved4 = 0;
    msg_nav_eng.reserved5 = 0;
    msg_nav_eng.reserved6 = 0;
    msg_nav_eng.reserved7 = 0;

    printf("Configuring Ublox NAV Engine ....\n");
    int send_status = SendUbxMsg(CLASS_CFG, NAV5_CFG, &msg_nav_eng, sizeof(CfgNavEng));

    if(GetUbxAck(CLASS_CFG, NAV5_CFG) < 0){
        std::cerr << "Could not configure Ublox NAV Engine\n";
        return -1;
    }
    printf("Successfully configured Ublox NAV Engine\n");

    return 1;
 }

 int UbxDriver::EnableNavPvt(){
    CfgMeasrate meas_msg_rate;
    meas_msg_rate.msg_class = CLASS_NAV;
    meas_msg_rate.msg_id = MSG_NAV_PVT;
    meas_msg_rate.msg_rate = 0x01;

    int send_status = SendUbxMsg(CLASS_CFG, MSG_CFG_RATE, &meas_msg_rate, sizeof(CfgMeasrate));

    printf("Enabling NAV PVT messages...\n");
    if(GetUbxAck(CLASS_CFG, MSG_CFG_RATE) < 0){
        std::cerr << "Could not enable NAV PVT messages\n";
        return -1;
    }
    printf("Successfully enabled NAV PVT messages...\n");

    return 1;
 }

 int UbxDriver::ConfigureSolutionRate(std::uint16_t meas_rate_ms,
                                      std::uint16_t nav_rate,
                                      std::uint16_t time_ref){
    CfgNavRate msg_nav_rate;
    msg_nav_rate.measure_rate = meas_rate_ms;
    msg_nav_rate.nav_rate     = nav_rate;
    msg_nav_rate.timeref      = time_ref;

    printf("Configuring sampling rate to %g Hz, nav rate to %d and setting time reference to %d...\n", 
            (1.0/(float)meas_rate_ms)*1000, nav_rate, time_ref);

    //Try 10 times to set desired rate
    int rate_set_count = 0;
    while(rate_set_count < 10){
        int send_status = SendUbxMsg(CLASS_CFG, CFG_SAMPLE_RATE, &msg_nav_rate, sizeof(CfgNavRate));
        if(GetUbxAck(CLASS_CFG, CFG_SAMPLE_RATE) < 0){
            printf("Failed to configure sampling rate to %g Hz, nav rate to %d and setting time reference to %d, try no: %d\n", 
                    (1.0/(float)meas_rate_ms)*1000, nav_rate, time_ref, rate_set_count + 1);
        }else{
            printf("Successfully configured sampling rate to %g Hz, nav rate to %d and setting time reference to %d, try no: %d\n", 
            (1.0/(float)meas_rate_ms)*1000, nav_rate, time_ref, rate_set_count + 1);
            return 1;
        }
        usleep(2000);
        rate_set_count++;
    }

    if(rate_set_count == 10){
        printf("Failed to configure sampling rate to %g Hz, nav rate to %d and setting time reference to %d, after %d tries\n", 
                    (1.0/(float)meas_rate_ms)*1000, nav_rate, time_ref, rate_set_count);
        return -1;
    }
    
    return 1;
}

/******************************************************************************************************************************************
 *********************************************************    NAV DATA   ******************************************************************
 ******************************************************************************************************************************************/

std::uint32_t UbxDriver::GetUnsigned32BitData(const std::vector<std::uint8_t> &data, size_t start_idx){
    return (std::uint32_t)((data[start_idx + 3] << 24) | (data[start_idx + 2] << 16) | (data[start_idx + 1] << 8) | (data[start_idx]));
}

std::int32_t UbxDriver::Get32BitData(const std::vector<std::uint8_t> &data, size_t start_idx){
    return ((data[start_idx + 3] << 24) | (data[start_idx + 2] << 16) | (data[start_idx + 1] << 8) | (data[start_idx]));
}

std::uint16_t UbxDriver::GetUnsigned16BitData(const std::vector<std::uint8_t> &data, size_t start_idx){
    return (std::uint16_t)((data[start_idx + 1] << 8) | (data[start_idx]));
}

std::int16_t UbxDriver::Get16BitData(const std::vector<std::uint8_t> &data, size_t start_idx){
    return ((data[start_idx + 1] << 8) | (data[start_idx]));
}

int UbxDriver::GetNavPvtData(NavPvtData &nav_pvt_data){

    std::vector<std::uint8_t> data;

    // Read frames till correct data is avalable
    size_t count_frames = 0;
    while(count_frames < 20000){
        if(DecodeSingleGenericMessage(data) == ((CLASS_NAV <<8) | MSG_NAV_PVT)){
            nav_pvt_data.iTow_ms = GetUnsigned32BitData(data, 6);
            nav_pvt_data.year = GetUnsigned16BitData(data, 10);
            nav_pvt_data.month = data[12];
            nav_pvt_data.day = data[13];
            nav_pvt_data.hour = data[14];
            nav_pvt_data.min = data[15];
            nav_pvt_data.sec = data[16];
            nav_pvt_data.valid = data[17];
            nav_pvt_data.tAcc_ns = GetUnsigned32BitData(data, 18);
            nav_pvt_data.sFrac_ns = Get32BitData(data, 22);
            nav_pvt_data.fixType = data[26];
            nav_pvt_data.flags = data[27];
            nav_pvt_data.flags2 = data[28];
            nav_pvt_data.numSV = data[29];
            nav_pvt_data.lon_rad = ((double)Get32BitData(data, 30) * 1e-7)/RAD2DEG;
            nav_pvt_data.lat_rad = ((double)Get32BitData(data, 34) * 1e-7)/RAD2DEG;
            nav_pvt_data.hElpsd_m = ((double)Get32BitData(data, 38))/1000;            
            nav_pvt_data.hMsl_m = ((double)Get32BitData(data, 42))/1000;  
            nav_pvt_data.hAcc_m = ((double)GetUnsigned32BitData(data, 46))/1000;  
            nav_pvt_data.vAcc_m = ((double)GetUnsigned32BitData(data, 50))/1000; 
            nav_pvt_data.velN_mps = ((double)Get32BitData(data, 54))/1000; 
            nav_pvt_data.velE_mps = ((double)Get32BitData(data, 58))/1000;
            nav_pvt_data.velD_mps = ((double)Get32BitData(data, 62))/1000;
            nav_pvt_data.gSpeed_mps = ((double)Get32BitData(data, 66))/1000;
            nav_pvt_data.headMot_rad = ((double)Get32BitData(data, 70)*1e-5)/RAD2DEG;
            nav_pvt_data.sAcc_mps = ((double)GetUnsigned32BitData(data, 74))/1000;
            nav_pvt_data.headAcc_rad = ((double)GetUnsigned32BitData(data, 78) * 1e-5)/RAD2DEG;
            nav_pvt_data.pDop = ((double)GetUnsigned16BitData(data, 82) * 0.01);
            nav_pvt_data.flags3 = GetUnsigned16BitData(data, 84);
            nav_pvt_data.headVeh_rad = ((double)Get32BitData(data, 90) * 1e-5)/RAD2DEG;
            nav_pvt_data.magDec_rad = ((double)Get16BitData(data, 94) * 1e-2)/RAD2DEG;
            nav_pvt_data.magAcc_rad =  ((double)GetUnsigned16BitData(data, 96) * 1e-2)/RAD2DEG;

            // printf("################################################\n");
            // printf("iTow: %d, min: %d, sec: %d, fixType: %d, fix validity: %d, lon: %.10f, lat: %.10f\n", nav_pvt_data.iTow_ms, 
            //     nav_pvt_data.min, nav_pvt_data.sec, nav_pvt_data.fixType, nav_pvt_data.flags & (1 << 0), nav_pvt_data.lon_rad, 
            //     nav_pvt_data.lat_rad);
            // printf("Vn: %g, Ve: %g, Vd: %g, H Ellpsoid: %g, H Msl: %g, Hz Acc: %g, vt Acc: %g\n", nav_pvt_data.velN_mps, 
            //     nav_pvt_data.velE_mps, nav_pvt_data.velD_mps, nav_pvt_data.hElpsd_m, nav_pvt_data.hMsl_m, nav_pvt_data.hAcc_m, 
            //     nav_pvt_data.vAcc_m);
            // printf("Mag Dec: %g, pDop: %g\n", nav_pvt_data.magDec_rad*RAD2DEG, nav_pvt_data.pDop);

            // int vd = nav_pvt_data.valid & (1 << 0);
            // int vt = nav_pvt_data.valid & (1 << 1);
            // int fr = nav_pvt_data.valid & (1 << 2);
            // int vm = nav_pvt_data.valid & (1 << 3);
            // printf("Valid Date: %d, Valid Time: %d, Fully Resolved: %d, Valid Mag Dec: %d\n", vd, vt, fr, vm);

            // int ca = nav_pvt_data.flags2 & (1 << 5);
            // int cd = nav_pvt_data.flags2 & (1 << 6);
            // int ct = nav_pvt_data.flags2 & (1 << 7);
            // printf("Conf Avail: %d, Conf Date: %d, Conf Time: %d\n", ca, cd, ct);

            return 1;
            break;
        }
        usleep(20000);
        count_frames++;
    }

    return -1;
}





