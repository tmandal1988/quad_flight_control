//
// File: stateEstimatorAutocode_private.h
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
#ifndef RTW_HEADER_stateEstimatorAutocode_private_h_
#define RTW_HEADER_stateEstimatorAutocode_private_h_
#include "rtwtypes.h"

extern const real32_T rtCP_pooled_KMmQGEgMZJGF[529];
extern const real32_T rtCP_pooled_RMz2mfhVMd4n[529];
extern const real32_T rtCP_pooled_h5N04J6mrCDP[324];
extern const real32_T rtCP_pooled_u8rExSs9SlNG[324];
extern const real32_T rtCP_pooled_LNk8l36ZWAmg[121];
extern const real32_T rtCP_pooled_HTJQPte3ACu1[64];
extern const real32_T rtCP_pooled_FQUpkqF8e3zU;

#define rtCP_processNoiseQ_Value       rtCP_pooled_KMmQGEgMZJGF  // Expression: processNoiseQ
                                                                 //  Referenced by: '<Root>/processNoiseQ'

#define rtCP_initCovP_Value            rtCP_pooled_RMz2mfhVMd4n  // Expression: initCovP
                                                                 //  Referenced by: '<Root>/initCovP'

#define rtCP_processNoiseNoGpsQ_Value  rtCP_pooled_h5N04J6mrCDP  // Expression: processNoiseNoGpsQ
                                                                 //  Referenced by: '<Root>/processNoiseNoGpsQ'

#define rtCP_initCovNoGpsP_Value       rtCP_pooled_u8rExSs9SlNG  // Expression: initCovNoGpsP
                                                                 //  Referenced by: '<Root>/initCovNoGpsP'

#define rtCP_measNoiseR_Value          rtCP_pooled_LNk8l36ZWAmg  // Expression: measNoiseR
                                                                 //  Referenced by: '<Root>/measNoiseR'

#define rtCP_measNoiseNoGpsR_Value     rtCP_pooled_HTJQPte3ACu1  // Expression: measNoiseNoGpsR
                                                                 //  Referenced by: '<Root>/measNoiseNoGpsR'

#define rtCP_gEarth_mps2_Value         rtCP_pooled_FQUpkqF8e3zU  // Computed Parameter: rtCP_gEarth_mps2_Value
                                                                 //  Referenced by: '<Root>/gEarth_mps2'

#endif                          // RTW_HEADER_stateEstimatorAutocode_private_h_

//
// File trailer for generated code.
//
// [EOF]
//
