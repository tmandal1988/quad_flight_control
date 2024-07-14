/*
GPS driver dcode is placed under the BSD license.
Written by Egor Fedorov (egor.fedorov@emlid.com)
Copyright (c) 2014, Emlid Limited
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Emlid Limited nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL EMLID LIMITED BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "Ublox.h"

#define PREAMBLE_OFFSET 2
// class UBXScanner

UBXScanner::UBXScanner()
{
    reset();
}

unsigned char* UBXScanner::getMessage()
{
    return message;
}

unsigned int UBXScanner::getMessageLength()
{
    return message_length;
}

unsigned int UBXScanner::getPosition()
{
    return position;
}

void UBXScanner::reset()
{
    message_length = 0;
    position = 0;
    state = Sync1;
}

int UBXScanner::update(unsigned char data)
{
    if (state != Done)
        message[position++] = data;

    switch (state)
    {
    case Sync1:
        if (data == 0xb5)
            state = Sync2;
        else
            reset();
        break;

    case Sync2:
        if (data == 0x62)
            state = Class;
        else
            if (data == 0xb5)
                state = Sync1;
            else
                reset();
        break;

    case Class:
        state = ID;
        break;

    case ID:
        state = Length1;
        break;

    case Length1:
        payload_length = data;
        state = Length2;
        break;

    case Length2:
        payload_length += data << 8;
        state = Payload;
        break;

    case Payload:
        if (position == payload_length + 6)
            state = CK_A;
        if (position >= 1022)
            reset();
        break;

    case CK_A:
        state = CK_B;
        break;

    case CK_B:
        message_length = 6 + payload_length + 2;
        state = Done;
        break;

    case Done:
    default:
        break;
    }

    return state;
}

// class UBXParser


UBXParser::UBXParser(UBXScanner* ubxsc) : scanner(ubxsc), message(ubxsc->getMessage())
{

}

// This function updates message length and end position in the scanner's buffer

void UBXParser::updateMessageData(){
    length = scanner->getMessageLength();
    position = scanner->getPosition();
}

// Function decodeMessage() returns 1 in case of a successful message verification.
// It looks through the buffer, checks 2 ubx protocol sync chars (0xb5 and 0x62),
// counts the checksum and if everything is alright, defines the message type by id.
// In this example we only decode two messages: Nav-Status and Nav-Posllh. It is possible to add more
// id cases

int UBXParser::decodeMessageGeneric(std::string spi_device_name, std::vector<unsigned char>& data)
{
    unsigned char to_gps_data = 0x00, from_gps_data = 0x00;
    int status;
    int count = 0;

    int flag=1;
    int pos;
    uint16_t id;
    uint8_t CK_A=0, CK_B=0;

    data.clear();

    while (count < UBX_BUFFER_LENGTH/2)
    {
        // From now on, we will send zeroes to the receiver, which it will ignore
        // However, we are simultaneously getting useful information from it
        SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

        // Scanner checks the message structure with every byte received
        status = scanner->update(from_gps_data);

        if (status == UBXScanner::Done)
        {
            // Once we have a full message we decode it and reset the scanner, making it look for another message
            // in the data stream, coming over SPI
            updateMessageData(); // get the length and end message coordinate from UBX scanner
            pos = position-length; // count the message start position

            // All UBX messages start with 2 sync chars: 0xb5 and 0x62

            if (*(message+pos)!=0xb5) flag=0;
            if (*(message+pos+1)!=0x62) flag=0;

            // Count the checksum

            for (unsigned int i=2;i<(length-2);i++){
                CK_A += *(message+pos+i);
                CK_B += CK_A;
            }

            if (CK_A != *(message+pos+length-2)) flag=0;
            if (CK_B != *(message+pos+length-1)) flag=0;

            // If the sync chars, or the checksum is wrong, we should not be doing this anymore

            if (flag==0){
                scanner->reset();
            }else{
                // If we got everything right, then return the raw data
                data.clear();
                data.push_back(*message + pos);
                data.push_back(*(message +pos + 1));
                data.push_back(*(message + pos + 2));
                data.push_back(*(message + pos + 3));
                data.push_back(*(message + pos + 4));
                data.push_back(*(message + pos + 5));

                for(size_t idx = 6; idx < length; idx++)
                    data.push_back(*(message + pos + idx));

                id = (*(message+pos+2)) << 8 | (*(message+pos+3)); // ID is a two-byte number with little endianness

                flag = id; // will return message id, if we got the info decoded
                scanner->reset();
                return flag;

            }         
        }

        count++;
    }

    return -1;    
}

// Function decodeMessage() returns 1 in case of a successful message verification.
// It looks through the buffer, checks 2 ubx protocol sync chars (0xb5 and 0x62),
// counts the checksum and if everything is alright, defines the message type by id.
// In this example we only decode two messages: Nav-Status and Nav-Posllh. It is possible to add more
// id cases

int UBXParser::decodeMessage(std::vector<double>& data)
{
    int flag=1;
    int pos;
    uint16_t id;
    uint8_t CK_A=0, CK_B=0;

    updateMessageData(); // get the length and end message coordinate from UBX scanner

    pos = position-length; // count the message start position

    // All UBX messages start with 2 sync chars: 0xb5 and 0x62

    if (*(message+pos)!=0xb5) flag=0;
    if (*(message+pos+1)!=0x62) flag=0;

    // Count the checksum

    for (unsigned int i=2;i<(length-2);i++){
        CK_A += *(message+pos+i);
        CK_B += CK_A;
    }

    if (CK_A != *(message+pos+length-2)) flag=0;
    if (CK_B != *(message+pos+length-1)) flag=0;

    // If the sync chars, or the checksum is wrong, we should not be doing this anymore

    if (flag==0) return 0;

    // If we got everything right, then it's time to decide, what type of message this is

    id = (*(message+pos+2)) << 8 | (*(message+pos+3)); // ID is a two-byte number with little endianness

    flag = id; // will return message id, if we got the info decoded

    switch(id){
        case 258:
                // ID for Nav-Posllh messages is 0x0102 == 258
                // In this example we extract 7 variables - longitude, latitude,
                // height above ellipsoid and mean sea level, horizontal and vertical
                // accuracy estimate and iTOW - GPS Millisecond Time of Week

                // All the needed parameters are 4-byte numbers with little endianness.
                // We know the current message and we want to update the info in the data vector.
                // First we clear the old data:

                data.clear();

                // Second, we extract the needed data from the message buffer and save it to the vector.

                //iTOW
                data.push_back ((unsigned)((*(message+pos+9) << 24) | (*(message+pos+8) << 16) | (*(message+pos+7) << 8) | (*(message+pos+6))));
                //Longitude
                data.push_back ((*(message+pos+13) << 24) | (*(message+pos+12) << 16) | (*(message+pos+11) << 8) | (*(message+pos+10)));
                //Latitude
                data.push_back ((*(message+pos+17) << 24) | (*(message+pos+16) << 16) | (*(message+pos+15) << 8) | (*(message+pos+14)));
                //Height above Ellipsoid
                data.push_back ((*(message+pos+21) << 24) | (*(message+pos+20) << 16) | (*(message+pos+19) << 8) | (*(message+pos+18)));
                //Height above mean sea level
                data.push_back ((*(message+pos+25) << 24) | (*(message+pos+24) << 16) | (*(message+pos+23) << 8) | (*(message+pos+22)));
                //Horizontal Accuracy Estateimate
                data.push_back ((unsigned)((*(message+pos+29) << 24) | (*(message+pos+28) << 16) | (*(message+pos+27) << 8) | (*(message+pos+26))));
                //Vertical Accuracy Estateimate
                data.push_back ((unsigned)((*(message+pos+33) << 24) | (*(message+pos+32) << 16) | (*(message+pos+31) << 8) | (*(message+pos+30))));
                break;

        case 276:
                // ID for Nav-HpPosllh messages is 0x0114 == 276
                // In this example we extract 11 variables - longitude, latitude,
                // height above ellipsoid and mean sea level, HP lon, HP Lat, HP Alt,
                // HP MSL, horizontal and vertical
                // accuracy estimate and iTOW - GPS Millisecond Time of Week

                // All the needed parameters are 4-byte numbers with little endianness.
                // We know the current message and we want to update the info in the data vector.
                // First we clear the old data:

                data.clear();

                // Second, we extract the needed data from the message buffer and save it to the vector.

                //iTOW
                data.push_back ((unsigned)((*(message+pos+13) << 24) | (*(message+pos+12) << 16) | (*(message+pos+11) << 8) | (*(message+pos+10))));
                //Longitude
                data.push_back ((*(message+pos+17) << 24) | (*(message+pos+16) << 16) | (*(message+pos+15) << 8) | (*(message+pos+14)));
                //Latitude
                data.push_back ((*(message+pos+21) << 24) | (*(message+pos+20) << 16) | (*(message+pos+19) << 8) | (*(message+pos+18)));
                //Height above Ellipsoid
                data.push_back ((*(message+pos+25) << 24) | (*(message+pos+24) << 16) | (*(message+pos+23) << 8) | (*(message+pos+22)));
                //Height above mean sea level
                data.push_back ((*(message+pos+29) << 24) | (*(message+pos+28) << 16) | (*(message+pos+27) << 8) | (*(message+pos+26)));
                //HP Lon
                data.push_back (*(message+pos+30));
                //HP Lat
                data.push_back (*(message+pos+31));
                //HP Height
                data.push_back (*(message+pos+32));
                //HP MSL
                data.push_back (*(message+pos+33));
                //Horizontal Accuracy Estateimate
                data.push_back ((unsigned)((*(message+pos+37) << 24) | (*(message+pos+36) << 16) | (*(message+pos+35) << 8) | (*(message+pos+34))));
                //Vertical Accuracy Estateimate
                data.push_back ((unsigned)((*(message+pos+41) << 24) | (*(message+pos+40) << 16) | (*(message+pos+39) << 8) | (*(message+pos+38))));
                break;

        case 259:
                // ID for Nav-Status messages is 0x0103 == 259
                // This message contains a lot of information, but the reason we use it the GPS fix status
                // This status is a one byte flag variable, with the offset of 10: first 6 bytes of the captured
                // message contain message header, id, payload length, next 4 are the first payload variable - iTOW.
                // We are not interested in it, so we just proceed to the fifth byte - gpsFix flag.
                // We are also interested in checking the gpsFixOk flag in the next byte/
                data.clear();
                // gpsFix
                data.push_back ((unsigned)(*(message+pos+10)));
                // flags, which contain gpsFixOk
                data.push_back ((unsigned)(*(message+pos+10)));

                break;


        case 274:
                // ID for Nav-Velned messages is 0x0112 == 274
                // In this example we extract 9 variables
                // 1. iTOW - GPS Millisecond Time of Week
                // 2. North velocity
                // 3. East velocity
                // 4. Down velocity
                // 5. 3D speed
                // 6. 2D speed - ground speed
                // 7. 2D heading
                // 8. Speed accuracy estimate
                // 9. Course/Heading accuracy estimate

                // All the needed parameters are 4-byte numbers with little endianness.
                // We know the current message and we want to update the info in the data vector.
                // First we clear the old data:

                data.clear();

                //iTOW
                data.push_back ((unsigned)((*(message+pos+9) << 24) | (*(message+pos+8) << 16) | (*(message+pos+7) << 8) | (*(message+pos+6))));
                //North velocity
                data.push_back ((*(message+pos+13) << 24) | (*(message+pos+12) << 16) | (*(message+pos+11) << 8) | (*(message+pos+10)));
                //East velocity
                data.push_back ((*(message+pos+17) << 24) | (*(message+pos+16) << 16) | (*(message+pos+15) << 8) | (*(message+pos+14)));
                //Down velocity
                data.push_back ((*(message+pos+21) << 24) | (*(message+pos+20) << 16) | (*(message+pos+19) << 8) | (*(message+pos+18)));
                //3D speed
                data.push_back ((*(message+pos+25) << 24) | (*(message+pos+24) << 16) | (*(message+pos+23) << 8) | (*(message+pos+22)));
                //2D speed - ground speed
                data.push_back ((*(message+pos+29) << 24) | (*(message+pos+28) << 16) | (*(message+pos+27) << 8) | (*(message+pos+26)));
                //2D heading
                data.push_back ((*(message+pos+33) << 24) | (*(message+pos+32) << 16) | (*(message+pos+31) << 8) | (*(message+pos+30)));
                //Speed accuracy estimate
                data.push_back ((unsigned)((*(message+pos+37) << 24) | (*(message+pos+36) << 16) | (*(message+pos+35) << 8) | (*(message+pos+34))));
                //Course/Heading accuracy estimate
                data.push_back ((unsigned)((*(message+pos+41) << 24) | (*(message+pos+40) << 16) | (*(message+pos+39) << 8) | (*(message+pos+38))));
                break;


        default:
                // In case we don't want to decode the received message

                flag = 0;

                break;
    }

    return flag;
}

// Function checkMessage() returns 1 if the message, currently stored in the buffer is valid.
// Validity means, that the necessary sync chars are present and the checksum test is passed

int UBXParser::checkMessage()
{
    int flag=1;
    int pos;
    uint8_t CK_A=0, CK_B=0;

    updateMessageData(); // get the length and end message coordinate from UBX scanner

    pos = position-length; // count the message start position

    // All UBX messages start with 2 sync chars: 0xb5 and 0x62

    if (*(message+pos)!=0xb5) flag=0;
    if (*(message+pos+1)!=0x62) flag=0;

    // Count the checksum

    for (unsigned int i=2;i<(length-2);i++){
        CK_A += *(message+pos+i);
        CK_B += CK_A;
    }

    if (CK_A != *(message+pos+length-2)) flag=0;
    if (CK_B != *(message+pos+length-1)) flag=0;

    return flag;
}

// class Ublox

Ublox::Ublox(std::string name) : spi_device_name(name), scanner(new UBXScanner()), parser(new UBXParser(scanner))
{

}

Ublox::Ublox(std::string name, UBXScanner* scan, UBXParser* pars) : spi_device_name(name), scanner(scan), parser(pars)
{

}

int Ublox::enableNAV_HPPOSLLH()
{
    CfgMeasrate msg_meas_rate;
    msg_meas_rate.msg_class = CLASS_NAV;
    msg_meas_rate.msg_id = MSG_HP_LLH;
    msg_meas_rate.msg_rate = 0x01;

    int send_status = _sendMessage(CLASS_CFG, MSG_CFG_RATE, &msg_meas_rate, sizeof(CfgMeasrate));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            // printf("Data at idx 2: %#04x, Data at idx 3: %#04x\n", raw_data[2], raw_data[3]);
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == MSG_CFG_RATE)
                ack_received = 1;
        }
        count++;
    }

    if (send_status > 0 && ack_received > 0)
        return 1;

    return -1;
}

int Ublox::enableNAV_POSLLH()
{
    CfgMeasrate msg_meas_rate;
    msg_meas_rate.msg_class = CLASS_NAV;
    msg_meas_rate.msg_id = MSG_LLH;
    msg_meas_rate.msg_rate = 0x01;

    int send_status = _sendMessage(CLASS_CFG, MSG_CFG_RATE, &msg_meas_rate, sizeof(CfgMeasrate));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == MSG_CFG_RATE)
                ack_received = 1;
        }
        count++;
    }

    if (send_status > 0 && ack_received > 0)
        return 1;

    return -1;
}

int Ublox::enableNAV_STATUS()
{
    CfgMeasrate msg_meas_rate;
    msg_meas_rate.msg_class = CLASS_NAV;
    msg_meas_rate.msg_id = MSG_SOL_STATUS;
    msg_meas_rate.msg_rate = 0x01;

    int send_status = _sendMessage(CLASS_CFG, MSG_CFG_RATE, &msg_meas_rate, sizeof(CfgMeasrate));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == MSG_CFG_RATE)
                ack_received = 1;
        }
        count++;
    }

    if (send_status > 0 && ack_received > 0)
        return 1;

    return -1;
}

int Ublox::enableNAV_VELNED()
{
    CfgMeasrate msg_meas_rate;
    msg_meas_rate.msg_class = CLASS_NAV;
    msg_meas_rate.msg_id = MSG_VEL_NED;
    msg_meas_rate.msg_rate = 0x01;

    int send_status = _sendMessage(CLASS_CFG, MSG_CFG_RATE, &msg_meas_rate, sizeof(CfgMeasrate));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == MSG_CFG_RATE)
                ack_received = 1;
        }
        count++;
    }

    if (send_status > 0 && ack_received > 0)
        return 1;

    return -1;
}

int Ublox::resetConfig(){
    ResetUblox msg_rst;
    msg_rst.nav_bbr_mask = 0; //Warm Start
    msg_rst.reset_mode = 1; // Controller Software Reset
    msg_rst.reserved = 0;
    int send_status1 = _sendMessage(CLASS_CFG, SOFT_RST_CFG, &msg_rst, sizeof(ResetUblox));

    if (send_status1 < 0)
        return -1;
    // Just wait and don't check for ACK from ublox as ACK after reset is not reliable
    usleep(500000);

    ResetCfgUblox msg_cf_rst;
    msg_cf_rst.clear_mask1     = 15;
    msg_cf_rst.clear_mask2     = 1;
    msg_cf_rst.clear_mask3     = 0;
    msg_cf_rst.clear_mask4     = 0;

    msg_cf_rst.save_mask1      = 0;
    msg_cf_rst.save_mask2      = 0;
    msg_cf_rst.save_mask3      = 0;
    msg_cf_rst.save_mask4      = 0;

    msg_cf_rst.load_mask1      = 15;
    msg_cf_rst.load_mask2      = 1;
    msg_cf_rst.load_mask3      = 0;
    msg_cf_rst.load_mask4      = 0;

    int send_status2 = _sendMessage(CLASS_CFG, RST_CFG, &msg_cf_rst, sizeof(ResetCfgUblox));

    std::vector<unsigned char> raw_data;
    int ack_received  = -1;
    if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
    {
        if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == RST_CFG)
            ack_received = 1;
    }
    

    if(send_status2 > 0 && ack_received > 0){
        return 1;
    }

    return -1;
}

int Ublox::testConnection()
{
    int status;
    int count = 0;
    unsigned char to_gps_data = 0x00, from_gps_data = 0x00;

    // reset Ublox config
    if(resetConfig() < 0)
    {
        std::cerr << "Could not reset ublox config over SPI\n";
        return 0;
    }

    // we do this, so that at least one ubx message is enabled

    // if (enableNAV_HPPOSLLH()<0)//Not supported
    // {
    //     std::cerr << "Could not configure HPPOSLLH for ublox over SPI\n";
    // }

    if(configureUbloxSpiPort()<0){
        std::cerr << "Could not configure SPI Port for ublox over SPI\n";
        return 0;
    }

    if (enableNAV_POSLLH()<0)
    {
        std::cerr << "Could not configure POSLLH for ublox over SPI\n";
        return 0;
    }

    if (enableNAV_STATUS()<0)
    {
        std::cerr << "Could not configure NAV_STATUS for ublox over SPI\n";
        return 0;
    }

    if (enableNAV_VELNED()<0)
    {
        std::cerr << "Could not configure VELNED for ublox over SPI\n";
        return 0;
    }

    if (configureNavEngine() < 0){
        std::cerr << "Could not configure NAV5 Engine for ublox over SPI\n";
        return 0;
    }

    // printf("Here\n");

    // while (count < UBX_BUFFER_LENGTH/2)
    // {
    //     // From now on, we will send zeroes to the receiver, which it will ignore
    //     // However, we are simultaneously getting useful information from it
    //     SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

    //     // Scanner checks the message structure with every byte received
    //     status = scanner->update(from_gps_data);

    //     if (status == UBXScanner::Done)
    //     {
    //         // Once we have a full message we decode it and reset the scanner, making it look for another message
    //         // in the data stream, coming over SPI

    //         // If we find at least one valid message in the buffer, we consider connection to be established
    //         if(parser->checkMessage()==1)
    //         {
    //             scanner->reset();
    //             return 1;
    //         }

    //         scanner->reset();
    //     }

    //     count++;
    // }

    return 1;
}

int Ublox::configureSolutionRate(std::uint16_t meas_rate,
                           std::uint16_t nav_rate,
                           std::uint16_t timeref)
{
    CfgNavRate msg;
    msg.measure_rate = meas_rate;
    msg.nav_rate     = nav_rate;
    msg.timeref      = timeref;

    int send_status = _sendMessage(CLASS_CFG, CFG_SAMPLE_RATE, &msg, sizeof(CfgNavRate));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == CFG_SAMPLE_RATE)
                ack_received = 1;
        }
        count++;
    }

    if(send_status > 0 && ack_received > 0){
        return 1;
    }

    return 0;
}

int Ublox::configureNavEngine(){
    CfgNavEng msg;
    msg.set_mask = 5;
    msg.dyn_model = 7;
    msg.fix_mode = 2;;
    msg.fixed_alt = 0;
    msg.fixed_alt_var = 0;
    msg.min_elev = 0;
    msg.dr_limit = 0;
    msg.p_dop = 0;
    msg.t_dop = 0;
    msg.p_acc = 0;
    msg.t_acc = 0;
    msg.static_hold_threshold = 0;
    msg.dgnss_timeout = 0;
    msg.cno_thresh_num_svs = 0;
    msg.cno_thresh = 0;
    msg.reserved1 = 0;
    msg.reserved2 = 0;
    msg.static_hold_max_dist = 0;
    msg.utc_standard = 0;
    msg.reserved3 = 0;
    msg.reserved4 = 0;
    msg.reserved5 = 0;
    msg.reserved6 = 0;
    msg.reserved7 = 0;

    int send_status = _sendMessage(CLASS_CFG, NAV5_CFG, &msg, sizeof(CfgNavEng));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == NAV5_CFG)
                ack_received = 1;
        }
        count++;
    }

    if(send_status > 0 && ack_received > 0){
        return 1;
    }

    return -1;
}

int Ublox::configureUbloxSpiPort(){
    CfgPrt msg;
    msg.port_id = 4;
    msg.reserved1 = 0;
    msg.tx_ready = 0;
    msg.spi_mode = 0;
    msg.reserved2 = 0;
    msg.reserved3 = 0;
    msg.reserved4 = 0;
    msg.reserved5 = 0;
    msg.in_proto_mask = 1;
    msg.out_proto_mask = 1;
    msg.flags = 0;
    msg.reserved6 = 0;
    msg.reserved7 = 0;

    int send_status = _sendMessage(CLASS_CFG, PRT_CFG, &msg, sizeof(CfgPrt));

    std::vector<unsigned char> raw_data;
    int ack_received = -1;
    int count = 0;
    while(ack_received != 1 && count < 100){
        if(parser->decodeMessageGeneric(spi_device_name, raw_data) > 0)
        {
            if(raw_data[2] == 0x05 && raw_data[3] == 0x01 && raw_data[6] == CLASS_CFG && raw_data[7] == PRT_CFG)
                ack_received = 1;
        }
        count++;
    }

    if(send_status > 0 && ack_received > 0){
        return 1;
    }

    return -1;
}

int Ublox::_sendMessage(std::uint8_t msg_class, std::uint8_t msg_id, void *msg, std::uint16_t size)
{
    unsigned char buffer[UBX_BUFFER_LENGTH];

    UbxHeader header;
    header.preamble1 = PREAMBLE1;
    header.preamble2 = PREAMBLE2;
    header.msg_class = msg_class;
    header.msg_id    = msg_id;
    header.length    = size;

    int offset = _spliceMemory(buffer, &header, sizeof(UbxHeader));
    offset = _spliceMemory(buffer, msg, size, offset);

    auto checksum = _calculateCheckSum(buffer, offset);
    offset = _spliceMemory(buffer, &checksum, sizeof(CheckSum), offset);

    return SPIdev::transfer(spi_device_name.c_str(), buffer, nullptr, offset);
}

int Ublox::_spliceMemory(unsigned char *dest, const void * const src, std::size_t size, int dest_offset)
{
    std::memmove(dest + dest_offset, src, size);
    return dest_offset + size;
}

Ublox::CheckSum Ublox::_calculateCheckSum(unsigned char *message, std::size_t size) {
    CheckSum checkSum;
    checkSum.CK_A = checkSum.CK_B = 0;

    for (int i = PREAMBLE_OFFSET; i < size; i++) {
        checkSum.CK_A += message[i];
        checkSum.CK_B += checkSum.CK_A;
    }
    return checkSum;
}

int Ublox::decodeMessages()
{
    int status;
    unsigned char to_gps_data = 0x00, from_gps_data = 0x00;
    std::vector<double> position_data;

    // if (enableNAV_HPPOSLLH()<0)
    // {
    //     std::cerr << "Could not configure HPPOSLLH in ublox over SPI\n";
    // }

    if (enableNAV_POSLLH()<0)
    {
        std::cerr << "Could not configure ublox over SPI\n";
    }

    if (enableNAV_STATUS()<0)
    {
        std::cerr << "Could not configure ublox over SPI\n";
    }

    while (true)
    {
        // From now on, we will send zeroes to the receiver, which it will ignore
        // However, we are simultaneously getting useful information from it
        SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

        // Scanner checks the message structure with every byte received
        status = scanner->update(from_gps_data);

        if (status == UBXScanner::Done)
        {
            // Once we have a full message we decode it and reset the scanner, making it look for another message
            // in the data stream, coming over SPI
            if (parser->decodeMessage(position_data) > 0)
            {
                // do something with the obtained data
                //
                // in case if NAV-POSLLH messages we can do this:
                // printf("decodeMessages(): \nCurrent location data:\n");
                // printf("GPS Millisecond Time of Week: %lf\n", position_data[0]/1000);
                // printf("Longitude: %lf\n", position_data[1]/10000000);
                // printf("Latitude: %lf\n", position_data[2]/10000000);
                // printf("Height above Ellipsoid: %.3lf m\n", pos_data[3]/1000);
                // printf("Height above mean sea level: %.3lf m\n", pos_data[4]/1000);
                // printf("Horizontal Accuracy Estateimate: %.3lf m\n", pos_data[5]/1000);
                // printf("Vertical Accuracy Estateimate: %.3lf m\n", pos_data[6]/1000);
                //
                // in case of NAV-STATUS messages we can do this:
                //
                // printf("Current GPS status:\n");
                // printf("gpsFixOk: %d\n", ((int)pos_data[1] & 0x01));
                //
                // printf("gps Fix status: ");
                // switch((int)pos_data[0]){
                //     case 0x00:
                //         printf("no fix\n");
                //         break;
                //
                //     case 0x01:
                //         printf("dead reckoning only\n");
                //         break;
                //
                //     case 0x02:
                //         printf("2D-fix\n");
                //         break;
                //
                //     case 0x03:
                //         printf("3D-fix\n");
                //         break;
                //
                //     case 0x04:
                //         printf("GPS + dead reckoning combined\n");
                //         break;
                //
                //     case 0x05:
                //         printf("Time only fix\n");
                //         break;
                //
                //     default:
                //         printf("Reserved value. Current state unknown\n");
                //         break;
                //
                // }
                //
                // printf("\n");

            }

            scanner->reset();
        }

        usleep(200);

    }

    return 0 ;
}

int Ublox::decodeSingleMessage(message_t msg, std::vector<double>& position_data)
{
    switch(msg){
        case NAV_HPPOSLLH:
            {
                uint16_t id = 0x0114;
                int status;
                int count = 0;
                unsigned char to_gps_data = 0x00, from_gps_data = 0x00;

                while (count < UBX_BUFFER_LENGTH/2)
                {
                    // From now on, we will send zeroes to the receiver, which it will ignore
                    // However, we are simultaneously getting useful information from it
                    SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

                    // Scanner checks the message structure with every byte received
                    status = scanner->update(from_gps_data);

                    if (status == UBXScanner::Done)
                    {
                        // Once we have a full message we decode it and reset the scanner, making it look for another message
                        // in the data stream, coming over SPI
                        if(parser->decodeMessage(position_data) == id)
                        {
                            // Now let's do something with the extracted information
                            // in case of NAV-POSLLH messages we can print the information like this:
                            // printf("GPS Millisecond Time of Week: %lf\n", position_data[0]/1000);
                            // printf("Longitude: %lf\n", position_data[1]/10000000);
                            // printf("Latitude: %lf\n", position_data[2]/10000000);
                            // printf("Height above Ellipsoid: %.3lf m\n", pos_data[3]/1000);
                            // printf("Height above mean sea level: %.3lf m\n", pos_data[4]/1000);
                            // printf("Horizontal Accuracy Estateimate: %.3lf m\n", pos_data[5]/1000);
                            // printf("Vertical Accuracy Estateimate: %.3lf m\n", pos_data[6]/1000);


                            // You can see ubx message structure in ublox reference manual

                            scanner->reset();

                            return 1;
                        }

                        scanner->reset();
                    }

                    count++;
                }

                return 0;
            }

        break;

        case NAV_POSLLH:
            {
                uint16_t id = 0x0102;
                int status;
                int count = 0;
                unsigned char to_gps_data = 0x00, from_gps_data = 0x00;

                while (count < UBX_BUFFER_LENGTH/2)
                {
                    // From now on, we will send zeroes to the receiver, which it will ignore
                    // However, we are simultaneously getting useful information from it
                    SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

                    // Scanner checks the message structure with every byte received
                    status = scanner->update(from_gps_data);

                    if (status == UBXScanner::Done)
                    {
                        // Once we have a full message we decode it and reset the scanner, making it look for another message
                        // in the data stream, coming over SPI
                        if(parser->decodeMessage(position_data) == id)
                        {
                            // Now let's do something with the extracted information
                            // in case of NAV-POSLLH messages we can print the information like this:
                            // printf("GPS Millisecond Time of Week: %lf\n", position_data[0]/1000);
                            // printf("Longitude: %lf\n", position_data[1]/10000000);
                            // printf("Latitude: %lf\n", position_data[2]/10000000);
                            // printf("Height above Ellipsoid: %.3lf m\n", pos_data[3]/1000);
                            // printf("Height above mean sea level: %.3lf m\n", pos_data[4]/1000);
                            // printf("Horizontal Accuracy Estateimate: %.3lf m\n", pos_data[5]/1000);
                            // printf("Vertical Accuracy Estateimate: %.3lf m\n", pos_data[6]/1000);


                            // You can see ubx message structure in ublox reference manual

                            scanner->reset();

                            return 1;
                        }

                        scanner->reset();
                    }

                    count++;
                }

                return 0;
            }

        break;

        case NAV_STATUS:
            {
                uint16_t id = 0x0103;
                int status;
                int count = 0;
                unsigned char to_gps_data = 0x00, from_gps_data = 0x00;

                while (count < UBX_BUFFER_LENGTH/2)
                {
                    // From now on, we will send zeroes to the receiver, which it will ignore
                    // However, we are simultaneously getting useful information from it
                    SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

                    // Scanner checks the message structure with every byte received
                    status = scanner->update(from_gps_data);

                    if (status == UBXScanner::Done)
                    {
                        // Once we have a full message we decode it and reset the scanner, making it look for another message
                        // in the data stream, coming over SPI
                        if(parser->decodeMessage(position_data) == id)
                        {
                            // Now let's do something with the extracted information
                            // in case of NAV-STATUS messages we can do this:
                            //
                            // printf("Current GPS status:\n");
                            // printf("gpsFixOk: %d\n", ((int)pos_data[1] & 0x01));
                            //
                            // printf("gps Fix status: ");
                            // switch((int)pos_data[0]){
                            //     case 0x00:
                            //         printf("no fix\n");
                            //         break;
                            //
                            //     case 0x01:
                            //         printf("dead reckoning only\n");
                            //         break;
                            //
                            //     case 0x02:
                            //         printf("2D-fix\n");
                            //         break;
                            //
                            //     case 0x03:
                            //         printf("3D-fix\n");
                            //         break;
                            //
                            //     case 0x04:
                            //         printf("GPS + dead reckoning combined\n");
                            //         break;
                            //
                            //     case 0x05:
                            //         printf("Time only fix\n");
                            //         break;
                            //
                            //     default:
                            //         printf("Reserved value. Current state unknown\n");
                            //         break;
                            //
                            // }
                            //
                            // printf("\n");

                            // You can see ubx message structure in ublox reference manual

                            scanner->reset();

                            return 1;
                        }

                        scanner->reset();
                    }

                    count++;
                }

                return 0;
            }

        break;

        case NAV_VELNED:
            {
                uint16_t id = 0x0112;
                int status;
                int count = 0;
                unsigned char to_gps_data = 0x00, from_gps_data = 0x00;

                while (count < UBX_BUFFER_LENGTH/2)
                {
                    // From now on, we will send zeroes to the receiver, which it will ignore
                    // However, we are simultaneously getting useful information from it
                    SPIdev::transfer(spi_device_name.c_str(), &to_gps_data, &from_gps_data, 1);

                    // Scanner checks the message structure with every byte received
                    status = scanner->update(from_gps_data);

                    if (status == UBXScanner::Done)
                    {
                        // Once we have a full message we decode it and reset the scanner, making it look for another message
                        // in the data stream, coming over SPI
                        if(parser->decodeMessage(position_data) == id)
                        {

                            // You can see ubx message structure in ublox reference manual

                            scanner->reset();

                            return 1;
                        }

                        scanner->reset();
                    }

                    count++;
                }

                return 0;
            }

        break;


        // add your ubx message type here!

        default:
            return 0;

        break;
    }
}
