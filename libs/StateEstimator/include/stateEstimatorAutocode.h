//
// File: stateEstimatorAutocode.h
//
// Code generated for Simulink model 'stateEstimatorAutocode'.
//
// Model version                  : 1.45
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:54:26 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM Cortex-M
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. ROM efficiency
// Validation result: Not run
//
#ifndef RTW_HEADER_stateEstimatorAutocode_h_
#define RTW_HEADER_stateEstimatorAutocode_h_
#include "rtwtypes.h"
#include "stateEstimatorAutocode_types.h"
#include "stateEstimator.h"

// Class declaration for model stateEstimatorAutocode
class stateEstimatorAutocode final
{
  // public data and function members
 public:
  // Constant parameters (default storage)
  struct ConstP_stateEstimatorAutocode_T {
    // Expression: imuNtchFiltParams
    //  Referenced by: '<Root>/imuNtchFilterParams'

    busImuNtchFiltParams imuNtchFilterParams_Value;

    // Expression: accelParams
    //  Referenced by: '<Root>/accelParams'

    busAccelParams accelParams_Value;

    // Expression: magParams
    //  Referenced by: '<Root>/magParams'

    busMagParams magParams_Value;

    // Computed Parameter: stateEstSmParams_Value
    //  Referenced by: '<Root>/stateEstSmParams'

    busStateEstSmParams stateEstSmParams_Value;

    // Expression: lidarParams
    //  Referenced by: '<Root>/lidarParams'

    busLidarParams lidarParams_Value;
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_stateEstimatorAutocode_T {
    busImuData imuData;                // '<Root>/imuData'
    busMagData magData;                // '<Root>/magData'
    busGpsData gpsData;                // '<Root>/gpsData'
    busBaroData baroData;              // '<Root>/baroData'
    busLidarData lidarData;            // '<Root>/lidarData'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_stateEstimatorAutocode_T {
    real32_T states[23];               // '<Root>/states'
    real32_T eulAng_rad[3];            // '<Root>/eulAng_rad'
    real32_T dcmNedToBody[9];          // '<Root>/dcmNedToBody'
    real32_T dcmNedToFep[9];           // '<Root>/dcmNedToFep'
    real32_T bodyAccels_mps2[3];       // '<Root>/bodyAccels_mps2'
    busStateEstimatorDebug stateEstimatorDebug;// '<Root>/stateEstimatorDebug'
  };

  // Copy Constructor
  stateEstimatorAutocode(stateEstimatorAutocode const&) = delete;

  // Assignment Operator
  stateEstimatorAutocode& operator= (stateEstimatorAutocode const&) & = delete;

  // Move Constructor
  stateEstimatorAutocode(stateEstimatorAutocode &&) = delete;

  // Move Assignment Operator
  stateEstimatorAutocode& operator= (stateEstimatorAutocode &&) = delete;

  // Root inports set method
  void setExternalInputs(const ExtU_stateEstimatorAutocode_T
    *pExtU_stateEstimatorAutocode_T);

  // Root outports get method
  const ExtY_stateEstimatorAutocode_T &getExternalOutputs() const;

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  stateEstimatorAutocode();

  // Destructor
  ~stateEstimatorAutocode();

  // private data and function members
 private:
  // External inputs
  ExtU_stateEstimatorAutocode_T stateEstimatorAutocode_U;

  // External outputs
  ExtY_stateEstimatorAutocode_T stateEstimatorAutocode_Y;

  // model instance variable for '<Root>/State Estimator'
  stateEstimator State_EstimatorMDLOBJ1;
};

// Constant parameters (default storage)
extern const stateEstimatorAutocode::ConstP_stateEstimatorAutocode_T
  stateEstimatorAutocode_ConstP;

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
//  '<Root>' : 'stateEstimatorAutocode'


//-
//  Requirements for '<Root>': stateEstimatorAutocode

#endif                                 // RTW_HEADER_stateEstimatorAutocode_h_

//
// File trailer for generated code.
//
// [EOF]
//
