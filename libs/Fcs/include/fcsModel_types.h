//
// File: fcsModel_types.h
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.48
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Mon May  1 19:34:21 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM 7
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. Traceability
// Validation result: All passed
//
#ifndef RTW_HEADER_fcsModel_types_h_
#define RTW_HEADER_fcsModel_types_h_
#include "rtwtypes.h"
#include <array>

// Model Code Variants
#ifndef DEFINED_TYPEDEF_FOR_enumStateMachine_
#define DEFINED_TYPEDEF_FOR_enumStateMachine_

// Defines the current state of the state machine
enum class enumStateMachine
  : int32_T {
  INACTIVE = 0,                        // Default value
  MTR_ARMED,
  INFLIGHT
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_enumFlightMode_
#define DEFINED_TYPEDEF_FOR_enumFlightMode_

// Defines the flight mode in use
enum class enumFlightMode
  : int32_T {
  STABILIZE = 0,                       // Default value
  VEL_CONTROL
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busRcInCmds_
#define DEFINED_TYPEDEF_FOR_busRcInCmds_

// Radio controller commands
struct busRcInCmds
{
  // Non dimensional throttle command from RC transmitter
  uint16_T throttleCmd_nd;

  // Non dimensional joystick Y command from RC transmitter
  uint16_T joystickYCmd_nd;

  // Non dimensional joystick X command from RC transmitter
  uint16_T joystickXCmd_nd;

  // Non dimensional joystick Z command from RC transmitter
  uint16_T joystickZCmd_nd;

  // Switch 1 from RC transmitter
  uint16_T rcSwitch1_nd;

  // Switch 2 from RC transmitter
  uint16_T rcSwitch2_nd;

  // Switch 3 from RC transmitter
  uint16_T rcSwitch3_nd;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busGeodeticPos_
#define DEFINED_TYPEDEF_FOR_busGeodeticPos_

// Bus containing position in lat, lon and altitude
struct busGeodeticPos
{
  // Lattitude in radian
  real_T lat_rad;

  // Longitude in radian
  real_T lon_rad;

  // Height above the ellipsoid in meter
  real_T alt_m;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busStateEstimate_
#define DEFINED_TYPEDEF_FOR_busStateEstimate_

// Bus containing process measurement data to be used by the controllers
struct busStateEstimate
{
  std::array<real_T, 3> attitude_rad;
  std::array<real_T, 3> bodyAngRates_radps;
  busGeodeticPos geodeticPos;
  std::array<real_T, 3> nedVel_mps;
  std::array<real_T, 9> ned2BodyDcm_nd;
  std::array<real_T, 9> ned2FepDcm_nd;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busPidParams_
#define DEFINED_TYPEDEF_FOR_busPidParams_

// Bus with all the params needed by the PID controller
struct busPidParams
{
  real_T Kp;
  real_T Ki;
  real_T Kd;
  real_T Kb;
  real_T Kt;
  real_T filterBandwidth_radps;
  std::array<real_T, 2> outputLimits;
  std::array<real_T, 2> outputRateLimits;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busSecondOrderFilterParam_
#define DEFINED_TYPEDEF_FOR_busSecondOrderFilterParam_

struct busSecondOrderFilterParam
{
  // Bandwidth of the second order filter
  real_T filterBandwidth_radps;

  // Damping ration of the filter
  real_T dampingRatio_nd;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busSignalConditioningParams_
#define DEFINED_TYPEDEF_FOR_busSignalConditioningParams_

// Bus to contain all the params required by the signal conditioning block
struct busSignalConditioningParams
{
  std::array<real_T, 2> filteredInputLimits;
  std::array<real_T, 2> filteredInputRateLimits;
  std::array<real_T, 2> filteredInputAccelLimits;
  std::array<real_T, 2> filteredInputJerkLimits;
  busSecondOrderFilterParam filterParams;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busAttCtrlParams_
#define DEFINED_TYPEDEF_FOR_busAttCtrlParams_

// Bus containing all the parameters necessary for attitude controller
struct busAttCtrlParams
{
  std::array<busPidParams, 3> ctrlParamsArray;
  std::array<busSignalConditioningParams, 3> cmdSignalConditioningParamsArray;
  std::array<busSignalConditioningParams, 3> measSignalConditioningParamsArray;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busAngRateCtrlParams_
#define DEFINED_TYPEDEF_FOR_busAngRateCtrlParams_

// Bus containing all the parameters necessary for angular rate controller
struct busAngRateCtrlParams
{
  std::array<busPidParams, 3> ctrlParamsArray;
  std::array<busSignalConditioningParams, 3> cmdSignalConditioningParamsArray;
  std::array<busSignalConditioningParams, 3> measSignalConditioningParamsArray;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busInnerLoopCtrlParams_
#define DEFINED_TYPEDEF_FOR_busInnerLoopCtrlParams_

// Bus containing parameters for inner loop controllers
struct busInnerLoopCtrlParams
{
  busAttCtrlParams attCtrlParams;
  busAngRateCtrlParams angRateCtrlParams;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busInnerLoopToAlloc_
#define DEFINED_TYPEDEF_FOR_busInnerLoopToAlloc_

// Bus containing data from inner loop to allocation block
struct busInnerLoopToAlloc
{
  // Commanded thrust
  real_T thrustCmd_N;

  // X moment command
  real_T xMomCmd_Nm;

  // Y moment command
  real_T yMomCmd_Nm;

  // Z moment command
  real_T zMomCmd_Nm;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busPidDebug_
#define DEFINED_TYPEDEF_FOR_busPidDebug_

// Bus with PID controller debug signals for post run analysis
struct busPidDebug
{
  // Final output of the PID controller
  real_T output;

  // Output portion due to propertional gain
  real_T proportionalOutput;

  // Output portion due to integral gain
  real_T integralOutput;

  // Output portion due to derivative gain
  real_T derivativeOutput;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busInnerLoopCtrlDebug_
#define DEFINED_TYPEDEF_FOR_busInnerLoopCtrlDebug_

// Bus containing inner loop controller debug data
struct busInnerLoopCtrlDebug
{
  std::array<busPidDebug, 3> angRateCtrlDebug;
  std::array<busPidDebug, 3> attCtrlDebug;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busFcsDebug_
#define DEFINED_TYPEDEF_FOR_busFcsDebug_

// Bus containing debug data from flight control system
struct busFcsDebug
{
  busInnerLoopCtrlDebug innerLoopCtrldebug;
  enumStateMachine state;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busOuterLoopCmds_
#define DEFINED_TYPEDEF_FOR_busOuterLoopCmds_

// Bus containing outer loop controller outputs
struct busOuterLoopCmds
{
  // Thrust command generated by the outer loop controller
  real_T thrustCmd_N;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busCtrlInputs_
#define DEFINED_TYPEDEF_FOR_busCtrlInputs_

// A generic bus to contains inputs to a controller
struct busCtrlInputs
{
  // Feed forward command to be fed into the controller
  real_T feedForwardCmd;

  // Controller commanded set point
  real_T cmd;

  // measurement data
  real_T meas;

  // Command that the controller should track in tracking mode
  real_T trackingCtrlCmd;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busAttCtrlInputs_
#define DEFINED_TYPEDEF_FOR_busAttCtrlInputs_

// Bus containing all the inputs necessary for attitude controllers
struct busAttCtrlInputs
{
  std::array<busCtrlInputs, 3> ctrlInputsArray;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busOuterLoopToInnerLoop_
#define DEFINED_TYPEDEF_FOR_busOuterLoopToInnerLoop_

// Bus containing all the data from outer loop to inner loop
struct busOuterLoopToInnerLoop
{
  busOuterLoopCmds outerLoopCmds;
  busAttCtrlInputs attCtrlInputs;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busAngRateCtrlInputs_
#define DEFINED_TYPEDEF_FOR_busAngRateCtrlInputs_

// Bus containing all the inputs necessary for angular rate controllers
struct busAngRateCtrlInputs
{
  std::array<busCtrlInputs, 3> ctrlInputsArray;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_busRcOutCmds_
#define DEFINED_TYPEDEF_FOR_busRcOutCmds_

// Bus containing data from RC Interpreter Library
struct busRcOutCmds
{
  // Dimensionalized command computed using throttle stick input.
  // Interpreted unit depends on the flight mode.
  real_T throttleStick;

  // Dimensionalized command computed using roll stick input.
  // Interpreted unit depends on the flight mode.
  real_T rollStick;

  // Dimensionalized command computed using pitch stick input.
  // Interpreted unit depends on the flight mode.
  real_T pitchStick;

  // Dimensionalized command computed using yaw stick input.
  // Interpreted unit depends on the flight mode.
  real_T yawStick;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_A1FBxl6bVqC2yteQsTJPYF_
#define DEFINED_TYPEDEF_FOR_struct_A1FBxl6bVqC2yteQsTJPYF_

struct struct_A1FBxl6bVqC2yteQsTJPYF
{
  std::array<real_T, 2> pitchRate_radps;
  std::array<real_T, 2> rollRate_radps;
  std::array<real_T, 2> yawRate_radps;
  std::array<real_T, 2> zForce_N;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_RgPxAEXAj72sD2xo4DfKOE_
#define DEFINED_TYPEDEF_FOR_struct_RgPxAEXAj72sD2xo4DfKOE_

struct struct_RgPxAEXAj72sD2xo4DfKOE
{
  std::array<real_T, 4> roll_nd;
  std::array<real_T, 4> pitch_nd;
  std::array<real_T, 4> yaw_nd;
  std::array<real_T, 4> throttle_nd;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_ek64hRkWPGOrPAEJD7wTRH_
#define DEFINED_TYPEDEF_FOR_struct_ek64hRkWPGOrPAEJD7wTRH_

struct struct_ek64hRkWPGOrPAEJD7wTRH
{
  struct_A1FBxl6bVqC2yteQsTJPYF cmdLimits;
  std::array<real_T, 2> pwmLimits;
  std::array<real_T, 2> lowPwmThreshold;
  struct_RgPxAEXAj72sD2xo4DfKOE coeffs;
};

#endif
#endif                                 // RTW_HEADER_fcsModel_types_h_

//
// File trailer for generated code.
//
// [EOF]
//
