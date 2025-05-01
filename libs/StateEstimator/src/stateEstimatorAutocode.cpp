//
// File: stateEstimatorAutocode.cpp
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
#include "stateEstimatorAutocode.h"
#include "stateEstimatorAutocode_types.h"
#include "rtwtypes.h"
#include "stateEstimatorAutocode_private.h"
#include "stateEstimator.h"

// Model step function
void stateEstimatorAutocode::step()
{
  // ModelReference: '<Root>/State Estimator' incorporates:
  //   Constant: '<Root>/accelParams'
  //   Constant: '<Root>/gEarth_mps2'
  //   Constant: '<Root>/imuNtchFilterParams'
  //   Constant: '<Root>/initCovNoGpsP'
  //   Constant: '<Root>/initCovP'
  //   Constant: '<Root>/lidarParams'
  //   Constant: '<Root>/magParams'
  //   Constant: '<Root>/measNoiseNoGpsR'
  //   Constant: '<Root>/measNoiseR'
  //   Constant: '<Root>/processNoiseNoGpsQ'
  //   Constant: '<Root>/processNoiseQ'
  //   Constant: '<Root>/stateEstSmParams'
  //   Inport: '<Root>/baroData'
  //   Inport: '<Root>/gpsData'
  //   Inport: '<Root>/imuData'
  //   Inport: '<Root>/lidarData'
  //   Inport: '<Root>/magData'
  //   Outport: '<Root>/bodyAccels_mps2'
  //   Outport: '<Root>/dcmNedToBody'
  //   Outport: '<Root>/dcmNedToFep'
  //   Outport: '<Root>/eulAng_rad'
  //   Outport: '<Root>/stateEstimatorDebug'
  //   Outport: '<Root>/states'

  State_EstimatorMDLOBJ1.step(&stateEstimatorAutocode_U.imuData,
    &stateEstimatorAutocode_U.magData, &stateEstimatorAutocode_U.gpsData,
    &stateEstimatorAutocode_U.baroData, &stateEstimatorAutocode_U.lidarData,
    &stateEstimatorAutocode_ConstP.imuNtchFilterParams_Value,
    &stateEstimatorAutocode_ConstP.accelParams_Value,
    &stateEstimatorAutocode_ConstP.magParams_Value,
    &stateEstimatorAutocode_ConstP.lidarParams_Value,
    &stateEstimatorAutocode_ConstP.stateEstSmParams_Value,
    &rtCP_processNoiseQ_Value[0], &rtCP_measNoiseR_Value[0],
    &rtCP_initCovP_Value[0], &rtCP_processNoiseNoGpsQ_Value[0],
    &rtCP_measNoiseNoGpsR_Value[0], &rtCP_initCovNoGpsP_Value[0],
    &rtCP_gEarth_mps2_Value, &stateEstimatorAutocode_Y.states[0],
    &stateEstimatorAutocode_Y.eulAng_rad[0],
    &stateEstimatorAutocode_Y.dcmNedToBody[0],
    &stateEstimatorAutocode_Y.dcmNedToFep[0],
    &stateEstimatorAutocode_Y.bodyAccels_mps2[0],
    &stateEstimatorAutocode_Y.stateEstimatorDebug);
}

// Model initialize function
void stateEstimatorAutocode::initialize()
{
  // SystemInitialize for ModelReference: '<Root>/State Estimator' incorporates:
  //   Outport: '<Root>/stateEstimatorDebug'

  State_EstimatorMDLOBJ1.init(&stateEstimatorAutocode_Y.stateEstimatorDebug);
}

// Model terminate function
void stateEstimatorAutocode::terminate()
{
  // (no terminate code required)
}

// Root inports set method
void stateEstimatorAutocode::setExternalInputs(const stateEstimatorAutocode::
  ExtU_stateEstimatorAutocode_T *pExtU_stateEstimatorAutocode_T)
{
  stateEstimatorAutocode_U = *pExtU_stateEstimatorAutocode_T;
}

// Root outports get method
const stateEstimatorAutocode::ExtY_stateEstimatorAutocode_T
  &stateEstimatorAutocode::getExternalOutputs() const
{
  return stateEstimatorAutocode_Y;
}

// Constructor
stateEstimatorAutocode::stateEstimatorAutocode():
  stateEstimatorAutocode_U(),
  stateEstimatorAutocode_Y()
{
  // Currently there is no constructor body generated.
}

// Destructor
stateEstimatorAutocode::~stateEstimatorAutocode()
{
  // Currently there is no destructor body generated.
}

//
// File trailer for generated code.
//
// [EOF]
//
