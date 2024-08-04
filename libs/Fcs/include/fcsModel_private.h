//
// File: fcsModel_private.h
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.112
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Sat Aug  3 00:09:24 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM 7
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. Traceability
// Validation result: All passed
//
#ifndef RTW_HEADER_fcsModel_private_h_
#define RTW_HEADER_fcsModel_private_h_
#include "rtwtypes.h"
#include "fcsModel.h"

extern real_T rt_urand_Upu32_Yd_f_pw(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw(uint32_T *u);
extern uint32_T plook_bincpag(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, uint32_T *prevIndex);
extern real_T intrp1d_la(uint32_T bpIndex, real_T frac, const real_T table[],
  uint32_T maxIndex);
extern uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T
  startIndex, uint32_T maxIndex);

#endif                                 // RTW_HEADER_fcsModel_private_h_

//
// File trailer for generated code.
//
// [EOF]
//
