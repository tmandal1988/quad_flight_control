//
// File: fcsModel_data.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.100
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Mon Apr  8 11:00:11 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM 7
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. Traceability
// Validation result: All passed
//
#include "fcsModel.h"

// Constant parameters (default storage)
const fcsModel::ConstP_fcsModel_T fcsModel_ConstP{
  // Pooled Parameter (Mixed Expressions)
  //  Referenced by:
  //    '<S3>/Constant'
  //    '<S155>/Constant'

  {
    {
      0.0                              // thrustCmd_N
    },                                 // outerLoopCmds

    {
      { {
          {
            0.0,                       // feedForwardCmd
            0.0,                       // cmd
            0.0,                       // meas
            false,                     // integratorReset
            0.0                        // trackingCtrlCmd
          }, {
            0.0,                       // feedForwardCmd
            0.0,                       // cmd
            0.0,                       // meas
            false,                     // integratorReset
            0.0                        // trackingCtrlCmd
          }, {
            0.0,                       // feedForwardCmd
            0.0,                       // cmd
            0.0,                       // meas
            false,                     // integratorReset
            0.0                        // trackingCtrlCmd
          } } }
      // ctrlInputsArray
    }                                  // attCtrlInputs
  },

  // Expression: allocationDataStruct.allocationMatrix
  //  Referenced by: '<S1>/Constant'

  { { -18260.93139809632, -18260.931398096309, -18260.931398096312,
      -18260.931398096312, 89077.714137205781, -89077.71413720581,
      -89077.714137205781, 89077.7141372058, 80444.631708577974,
      80444.631708577959, -80444.631708577988, -80444.631708577974,
      -152414.0730882342, 152414.0730882342, -152414.0730882342,
      152414.07308823423 } },

  // Expression: vehicleConstants.inertia_kgm2
  //  Referenced by: '<S2>/Constant'

  { { 0.02, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.03 } },

  // Computed Parameter: Constant_Value_e
  //  Referenced by: '<S9>/Constant'

  { { 1U, 2U, 3U } }
};

//
// File trailer for generated code.
//
// [EOF]
//
