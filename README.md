# Custom Flight Controller For Quadrotor
Instead of using open source flight controller stacks, this repo contains a custom flight controller for a quadrotor built from ground up. This is to have better control and understanding of the flight code and easy extensibility for personal hobby and research projects in the future.

## Hardware
 1. [Navio2 + Raspberry Pi](https://docs.emlid.com/navio2/) -  With SMP PREEMPT Emlid Raspberry Pi OS Buster image.
2.  4 x # SunnySky X Series V3 X2216 V3 880KV motors
3. 1 x FrSky Taranis Q X7 Digital Telemetry Radio System 2.4GHz ACCS
4. 1 x FrSky X8R Receiver
5. 4 x any 60A ESC with high update rate
6. 1 x 4000mAh 4S 30C Lipo
7. PM02D Power Module or any other DC-DC stepdown converter with 5.2V and 3A output. (See Navio2 doc for various options to power the board and the Pi)

## Setup
X-Frame quadrotor configuration

<img width="195" alt="Screenshot 2024-05-11 at 9 09 55 AM" src="https://github.com/tmandal1988/quad_flight_control/assets/7689397/38f93127-c2c4-4f71-9851-1927e8ad933a">

1. Connect the SBUS output from the receiver to the pin labeled PPM/SB on Navio2 board.
2. Tx-Rx Channel Mapping
    1. Chn 1 -> Roll Stick
    2. Chn 2 -> Pitch Stick
    3. Chn 3 -> Throttle Stick
    4. Chn 4 -> Rudder Stick
    5. Chn 5 -> 3 Position Toggle Switch. Position 1 < 1100, 1100 <= Position 2 < 1700, Position 3 >= 1700
3. Connect Mtr 1, 2, 3, 4 ESCs to their corresponding PWM port numbers on Navio2.

## Build and Run
1. ssh into raspberry pi
2. Clone the repository and `cd` into it
3. `mkdir build`
4. `cmake -S . -B build`
5. `make -C build/`
6. Use tmux and in a new session run `sudo ./build/QUADFLIGHTCONTROL`
7. It will take some time to initialize the IMU and the EKF, if there is no GPS lock it will initialize Mahony Filter for attitude estimation instead of EKF.
8. Better to log out of raspberry pi before flying
9. To stop the code ssh back into raspberry pi and rejoin the tmux session and use `Ctrl + C`
10. Flight data is saved in binary file `data_file.dat` which can be decoded using the .m script in extras folder. fcsDebug field names can be found in `libs\CommonUtils\src\fcs_out_data.cc`

## Flying
1. Power up everything.
2. Set all the trim tabs on the transmitter to zero.
3. Hold Throttle stick, rudder stick, roll stick and pitch stick at the positions where each of their PWM values are less than 1000 for 3 sec.
4. If flight is not commenced within 30 seconds, repeat step 3 again.
5. Toggle Switch Position 1 -> Stabilize Mode, Roll and Pitch Hold and Yaw Rate Control
6. Toggle Switch Position 2 -> Altitude Hold Mode, attitude control, altitude hold when Throttle stick is appoximately centered and heading hold when Yaw stick is appoximately centered, 
    1. Transition from Stabilize Mode to Altitude Hold Mode and vice versa while in flight can cause abrupt change in altitude. Altitude Hold Mode flights are carried out from takeoff to landing without mode change to Stabilize Mode.
    2. In Altitude Hold Mode moving Throttle stick below mid point will cause quadrotor to descend and moving Throttle stick above mid point will cause quadrotor to ascend.
7. Toggle Switch Position 3 -> Position Hold Mode, Position hold when all sticks (except Yaw stick) are centered. **Not Tuned Yet**
    1. Tested Procedure -> Takeoff in Altitude Holde Mode and when in stable flight switch to Position Hold Mode with all stick centered.
    2. Moving Sticks away from center position will move the quadrotor in the same direction.
