//
// File: fcsModel_data.cpp
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
#include "fcsModel.h"

// Constant parameters (default storage)
const fcsModel::ConstP_fcsModel_T fcsModel_ConstP{
  // Expression: rcParamsStruct
  //  Referenced by: '<S4>/State Machine'

  {
    {
      { { -0.52359877559829882, 0.52359877559829882 } },

      { { -0.52359877559829882, 0.52359877559829882 } },

      { { -0.52359877559829882, 0.52359877559829882 } },

      { { 0.0, 50.0 } }
    },

    { { 990.0, 2006.0 } },

    { { 900.0, 990.0 } },

    {
      { { 3.9045645325900015e-9, -1.7538875049548508e-5, 0.027221770655502103,
          -14.548300722291762 } },

      { { 3.9045645325900015e-9, -1.7538875049548508e-5, 0.027221770655502103,
          -14.548300722291762 } },

      { { 3.9045645325900015e-9, -1.7538875049548508e-5, 0.027221770655502103,
          -14.548300722291762 } },

      { { 1.8913975371836569e-9, -8.4651009094729818e-6, 0.013109462896716882,
          -6.5014878863088388 } }
    }
  },

  // Computed Parameter: Constant_Value_h
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
            0.0                        // trackingCtrlCmd
          }, {
            0.0,                       // feedForwardCmd
            0.0,                       // cmd
            0.0,                       // meas
            0.0                        // trackingCtrlCmd
          }, {
            0.0,                       // feedForwardCmd
            0.0,                       // cmd
            0.0,                       // meas
            0.0                        // trackingCtrlCmd
          } } }
      // ctrlInputsArray
    }                                  // attCtrlInputs
  },

  // Expression: allocationDataStruct.allocationMatrix
  //  Referenced by: '<S1>/Constant'

  { { -933.47739593970618, -933.47739593970618, -933.47739593970618,
      -933.47739593970618, -5280.5455741165442, -5280.5455741165442,
      5280.5455741165442, 5280.5455741165442, -5280.5455741165442,
      5280.5455741165442, 5280.5455741165442, -5280.5455741165442,
      -30625.898816586381, 30625.898816586381, -30625.898816586381,
      30625.898816586381 } },

  // Expression: vehicleConstants.inertia_kgm2
  //  Referenced by: '<S2>/Constant'

  { { 0.03, 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.0, 0.06 } }
};

//
// File trailer for generated code.
//
// [EOF]
//
