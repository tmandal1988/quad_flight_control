//
// File: stateEstimator.h
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM Cortex-M
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. ROM efficiency
// Validation result: Not run
//
#ifndef RTW_HEADER_stateEstimator_h_
#define RTW_HEADER_stateEstimator_h_
#include "rtwtypes.h"
#include "stateEstimator_types.h"
#include <cstring>

// Class declaration for model stateEstimator
class stateEstimator final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for model 'stateEstimator'
  struct DW_stateEstimator_T {
    real_T latLonAltOut[3];            // '<Root>/estimatorStateMachine'
    real_T refLatLonAlt[3];            // '<Root>/estimatorStateMachine'
    real_T gpsIdx;                     // '<Root>/estimatorStateMachine'
    real32_T Product[3];               // '<S4>/Product'
    real32_T Divide[3];                // '<S9>/Divide'
    real32_T TmpSignalConversionAtSFunctionI[3];// '<Root>/estimatorStateMachine' 
    real32_T initialStates[23];        // '<Root>/estimatorStateMachine'
    real32_T initialDcmBodyToNed[9];   // '<Root>/estimatorStateMachine'
    real32_T bodyRatesOut_radps[3];    // '<Root>/estimatorStateMachine'
    real32_T normMagVecOut_nd[3];      // '<Root>/estimatorStateMachine'
    real32_T UnitDelay_DSTATE[23];     // '<Root>/Unit Delay'
    real32_T XAxis_states[2];          // '<S23>/X Axis'
    real32_T XAxis1_states[2];         // '<S23>/X Axis1'
    real32_T XAxis2_states[2];         // '<S23>/X Axis2'
    real32_T XAxis_states_b[2];        // '<S24>/X Axis'
    real32_T XAxis1_states_m[2];       // '<S24>/X Axis1'
    real32_T XAxis2_states_h[2];       // '<S24>/X Axis2'
    real32_T Delay_DSTATE[23];         // '<S1>/Delay'
    real32_T Delay2_DSTATE[9];         // '<S1>/Delay2'
    real32_T UnitDelay1_DSTATE[3];     // '<Root>/Unit Delay1'
    real32_T gyroBias_radps[3];        // '<Root>/estimatorStateMachine'
    real32_T initialQuat[4];           // '<Root>/estimatorStateMachine'
    real32_T imuM2[6];                 // '<Root>/estimatorStateMachine'
    real32_T imuMean[6];               // '<Root>/estimatorStateMachine'
    real32_T magM2[3];                 // '<Root>/estimatorStateMachine'
    real32_T magMean[3];               // '<Root>/estimatorStateMachine'
    real32_T accelBias_mps2[3];        // '<Root>/estimatorStateMachine'
    real32_T magBias_nd[3];            // '<Root>/estimatorStateMachine'
    real32_T nedMagVecNorm_nd[3];      // '<Root>/estimatorStateMachine'
    real32_T covP[324];                // '<S1>/EKF NO GPS'
    real32_T covP_a[529];              // '<S1>/EKF'
    real32_T Divide1;                  // '<S10>/Divide1'
    real32_T stateEstInitPct;          // '<Root>/estimatorStateMachine'
    real32_T UnitDelay_DSTATE_c;       // '<S7>/Unit Delay'
    real32_T XAxis_tmp;                // '<S23>/X Axis'
    real32_T XAxis1_tmp;               // '<S23>/X Axis1'
    real32_T XAxis2_tmp;               // '<S23>/X Axis2'
    real32_T XAxis_tmp_k;              // '<S24>/X Axis'
    real32_T XAxis1_tmp_m;             // '<S24>/X Axis1'
    real32_T XAxis2_tmp_g;             // '<S24>/X Axis2'
    real32_T imuIdx;                   // '<Root>/estimatorStateMachine'
    real32_T magIdx;                   // '<Root>/estimatorStateMachine'
    real32_T baroIdx;                  // '<Root>/estimatorStateMachine'
    real32_T baroBias_m;               // '<Root>/estimatorStateMachine'
    real32_T baroInitAltM2;            // '<Root>/estimatorStateMachine'
    real32_T baroInitAltMean;          // '<Root>/estimatorStateMachine'
    int32_T durationCounter_1;         // '<Root>/estimatorStateMachine'
    int32_T durationCounter_1_c;       // '<Root>/estimatorStateMachine'
    int32_T durationCounter_1_p;       // '<Root>/estimatorStateMachine'
    uint16_T gpsValidCount;            // '<Root>/estimatorStateMachine'
    uint8_T is_active_c3_stateEstimator;// '<Root>/estimatorStateMachine'
    uint8_T is_c3_stateEstimator;      // '<Root>/estimatorStateMachine'
    boolean_T resetStates;             // '<Root>/estimatorStateMachine'
    boolean_T DelayInput1_DSTATE;      // '<S26>/Delay Input1'
    boolean_T icLoad;                  // '<S1>/Delay'
    boolean_T icLoad_j;                // '<S1>/Delay2'
    boolean_T isAttInitialized;        // '<Root>/estimatorStateMachine'
    boolean_T isBaroInitialized;       // '<Root>/estimatorStateMachine'
    boolean_T isPosInitialized;        // '<Root>/estimatorStateMachine'
    boolean_T covP_not_empty;          // '<S1>/EKF NO GPS'
    boolean_T covP_not_empty_p;        // '<S1>/EKF'
  };

  // Initial conditions function
  void init(busStateEstimatorDebug *rty_stateEstimatorDebug);

  // Copy Constructor
  stateEstimator(stateEstimator const&) = delete;

  // Assignment Operator
  stateEstimator& operator= (stateEstimator const&) & = delete;

  // Move Constructor
  stateEstimator(stateEstimator &&) = delete;

  // Move Assignment Operator
  stateEstimator& operator= (stateEstimator &&) = delete;

  // model step function
  void step(const busImuData *rtu_imuData, const busMagData *rtu_magData, const
            busGpsData *rtu_gpsData, const busBaroData *rtu_baroData, const
            busLidarData *rtu_lidarData, const busImuNtchFiltParams
            *rtu_imuNotchFiltParams, const busAccelParams *rtu_accelParams,
            const busMagParams *rtu_magParams, const busLidarParams
            *rtu_lidarParams, const busStateEstSmParams *rtu_stateEstSmParams,
            const real32_T rtu_processNoiseQ[529], const real32_T
            rtu_measNoiseR[121], const real32_T rtu_initCovP[529], const
            real32_T rtu_processNoiseNoGpsQ[324], const real32_T
            rtu_measNoiseNoGpsR[64], const real32_T rtu_initCovNoGpsP[324],
            const real32_T *rtu_gEarth_mps2, real32_T rty_states[23], real32_T
            rty_eulAng_rad[3], real32_T rty_dcmNedToBody[9], real32_T
            rty_dcmNedToFep[9], real32_T rty_bodyAccels_mps2[3],
            busStateEstimatorDebug *rty_stateEstimatorDebug);

  // Constructor
  stateEstimator();

  // Destructor
  ~stateEstimator();

  // private data and function members
 private:
  // Block states
  DW_stateEstimator_T stateEstimator_DW;

  // private member function(s) for subsystem '<Root>/TmpModelReferenceSubsystem'
  void stateEstimator_INITIALIZE(enumStateEstimateMode *mode, boolean_T
    *isMagValid, boolean_T *isGpsValid, real32_T *baroAltOut_m, boolean_T
    *isBaroValid, real32_T bodyAccelsOut_mps2[3], const busMagData *rtu_magData,
    const busGpsData *rtu_gpsData, const busBaroData *rtu_baroData, const
    busStateEstSmParams *rtu_stateEstSmParams);
  void state_enter_atomic_RUN_INIT_GPS(enumStateEstimateMode *mode, boolean_T
    *isMagValid, boolean_T *isGpsValid, real32_T *baroAltOut_m, boolean_T
    *isBaroValid, real32_T bodyAccelsOut_mps2[3], const busMagData *rtu_magData,
    const busGpsData *rtu_gpsData, const busBaroData *rtu_baroData);
};

//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'stateEstimator'
//  '<S1>'   : 'stateEstimator/EKF'
//  '<S2>'   : 'stateEstimator/Quaternions to Rotation Angles'
//  '<S3>'   : 'stateEstimator/Subsystem Reference'
//  '<S4>'   : 'stateEstimator/accelCorrection'
//  '<S5>'   : 'stateEstimator/estimatorStateMachine'
//  '<S6>'   : 'stateEstimator/eulToDcm'
//  '<S7>'   : 'stateEstimator/latLonAltToNedPos'
//  '<S8>'   : 'stateEstimator/lidarRangeToAgl'
//  '<S9>'   : 'stateEstimator/magCorrection'
//  '<S10>'  : 'stateEstimator/pressureToAlt'
//  '<S11>'  : 'stateEstimator/EKF/EKF'
//  '<S12>'  : 'stateEstimator/EKF/EKF NO GPS'
//  '<S13>'  : 'stateEstimator/EKF/State Selector'
//  '<S14>'  : 'stateEstimator/Quaternions to Rotation Angles/Angle Calculation'
//  '<S15>'  : 'stateEstimator/Quaternions to Rotation Angles/Quaternion Normalize'
//  '<S16>'  : 'stateEstimator/Quaternions to Rotation Angles/Angle Calculation/Protect asincos input'
//  '<S17>'  : 'stateEstimator/Quaternions to Rotation Angles/Angle Calculation/Protect asincos input/If Action Subsystem'
//  '<S18>'  : 'stateEstimator/Quaternions to Rotation Angles/Angle Calculation/Protect asincos input/If Action Subsystem1'
//  '<S19>'  : 'stateEstimator/Quaternions to Rotation Angles/Angle Calculation/Protect asincos input/If Action Subsystem2'
//  '<S20>'  : 'stateEstimator/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus'
//  '<S21>'  : 'stateEstimator/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
//  '<S22>'  : 'stateEstimator/Subsystem Reference/IMU Filters'
//  '<S23>'  : 'stateEstimator/Subsystem Reference/IMU Filters/Accel Notch Filters'
//  '<S24>'  : 'stateEstimator/Subsystem Reference/IMU Filters/Gyro Notch Filters'
//  '<S25>'  : 'stateEstimator/latLonAltToNedPos/Compare To Constant'
//  '<S26>'  : 'stateEstimator/latLonAltToNedPos/Detect Rise Positive'
//  '<S27>'  : 'stateEstimator/latLonAltToNedPos/convertLlhToNedPos'
//  '<S28>'  : 'stateEstimator/latLonAltToNedPos/Detect Rise Positive/Positive'


//-
//  Requirements for '<Root>': stateEstimator

#endif                                 // RTW_HEADER_stateEstimator_h_

//
// File trailer for generated code.
//
// [EOF]
//
