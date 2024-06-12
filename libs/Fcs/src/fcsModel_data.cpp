//
// File: fcsModel_data.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.103
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Mon Jun 10 17:10:25 2024
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
  //    '<S156>/Constant'

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

  // Expression: rpmToPwmLut(:, 1)
  //  Referenced by: '<S7>/Prelookup'

  { { 2250.0, 3138.0, 3950.0, 4635.0, 5170.0, 5660.0, 6200.0, 6580.0, 7010.0,
      7450.0, 7880.0, 8300.0, 8680.0, 9000.0, 9325.0 } },

  // Expression: rpmToPwmLut(:, 2)
  //  Referenced by: '<S7>/Interpolation Using Prelookup'

  { { 1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0, 1550.0, 1600.0,
      1650.0, 1700.0, 1750.0, 1800.0, 1850.0, 1900.0 } },

  // Expression: allocationDataStruct.allocationMatrix
  //  Referenced by: '<S1>/Constant'

  { { -18260.932199467716, -18260.932199467723, -18260.932199467705,
      -18260.932199467716, 89077.718047137911, -89077.71804713794,
      -89077.718047137882, 89077.718047137911, 80444.635240298026,
      80444.635240298012, -80444.635240298012, -80444.635240298026,
      -152414.08813250146, 152414.08813250143, -152414.08813250143,
      152414.08813250146 } },

  // Expression: vehicleConstants.inertia_kgm2
  //  Referenced by: '<S2>/Constant'

  { { 0.02, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.03 } },

  // Computed Parameter: Constant_Value_e
  //  Referenced by: '<S10>/Constant'

  { { 1U, 2U, 3U } }
};

//
// File trailer for generated code.
//
// [EOF]
//
