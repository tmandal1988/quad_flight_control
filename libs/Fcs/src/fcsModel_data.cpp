//
// File: fcsModel_data.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.51
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Fri Sep 22 15:02:49 2023
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
  // Computed Parameter: Constant_Value
  //  Referenced by: '<S3>/Constant'

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

  { { -9334.7739593088463, -9334.77395930884, -9334.773959308839,
      -9334.7739593088445, 45535.4827282571, -45535.4827282571,
      -45535.4827282571, 45535.4827282571, 41122.352243582813,
      41122.352243582813, -41122.352243582813, -41122.352243582813,
      -61251.797633172762, 61251.797633172762, -61251.797633172762,
      61251.797633172762 } },

  // Expression: vehicleConstants.inertia_kgm2
  //  Referenced by: '<S2>/Constant'

  { { 0.03, 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, 0.06 } }
};

//
// File trailer for generated code.
//
// [EOF]
//
