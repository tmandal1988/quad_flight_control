//
// File: computeMagMeasJacNoGps_mcxOPfaA.cpp
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
#include "rtwtypes.h"
#include "computeMagMeasJacNoGps_mcxOPfaA.h"
#include <cstring>

//
// Function for MATLAB Function: '<S1>/EKF NO GPS'
// function magMeasJac = computeMagMeasJacNoGps(states)
// COMPUTEMAGMEASJACNOGPS computes the measurement jacobian for magnetometer
// states in the estimator
//
// Inputs:
// states:                EKF states
//
// Outputs:
// magMeasJac:            3x18 Mag Meas Jacobian
//
void computeMagMeasJacNoGps_mcxOPfaA(const real32_T states[18], real32_T
  magMeasJac[54])
{
  real32_T tmp1;
  real32_T tmp10;
  real32_T tmp11;
  real32_T tmp12;
  real32_T tmp13;
  real32_T tmp14;
  real32_T tmp15;
  real32_T tmp16;
  real32_T tmp17;
  real32_T tmp18;
  real32_T tmp19;
  real32_T tmp2;
  real32_T tmp20;
  real32_T tmp21;
  real32_T tmp22;
  real32_T tmp25;
  real32_T tmp3;
  real32_T tmp4;
  real32_T tmp5;
  real32_T tmp6;
  real32_T tmp7;
  real32_T tmp8;

  // 'computeMagMeasJacNoGps:11' magMeasJac = zeros(3, 18, 'single');
  std::memset(&magMeasJac[0], 0, 54U * sizeof(real32_T));

  // Extract quat states
  // 'computeMagMeasJacNoGps:14' q0 = states(1);
  // 'computeMagMeasJacNoGps:15' q1 = states(2);
  // 'computeMagMeasJacNoGps:16' q2 = states(3);
  // 'computeMagMeasJacNoGps:17' q3 = states(4);
  // Extract Est NED Mag
  // 'computeMagMeasJacNoGps:19' magN = states(12);
  // 'computeMagMeasJacNoGps:20' magE = states(13);
  // 'computeMagMeasJacNoGps:21' magD = states(14);
  // 'computeMagMeasJacNoGps:23' tmp1 = 2*magD;
  tmp1 = 2.0F * states[13];

  // 'computeMagMeasJacNoGps:24' tmp2 = q2*tmp1;
  tmp2 = states[2] * tmp1;

  // 'computeMagMeasJacNoGps:25' tmp3 = q3*tmp1;
  tmp3 = states[3] * tmp1;

  // 'computeMagMeasJacNoGps:26' tmp4 = 2*magE;
  tmp4 = 2.0F * states[12];

  // 'computeMagMeasJacNoGps:27' tmp5 = q2*tmp4;
  tmp5 = states[2] * tmp4;

  // 'computeMagMeasJacNoGps:28' tmp6 = q0*tmp1;
  tmp6 = states[0] * tmp1;

  // 'computeMagMeasJacNoGps:29' tmp7 = q1*tmp4;
  tmp7 = states[1] * tmp4;

  // 'computeMagMeasJacNoGps:30' tmp8 = 4*magN;
  tmp8 = 4.0F * states[11];

  // 'computeMagMeasJacNoGps:31' tmp9 = q1*tmp1;
  tmp1 *= states[1];

  // 'computeMagMeasJacNoGps:32' tmp10 = q0*tmp4;
  tmp10 = states[0] * tmp4;

  // 'computeMagMeasJacNoGps:33' tmp11 = 2*q2^2;
  tmp11 = states[2] * states[2] * 2.0F;

  // 'computeMagMeasJacNoGps:34' tmp12 = 2*q3^2 - 1;
  tmp12 = states[3] * states[3] * 2.0F - 1.0F;

  // 'computeMagMeasJacNoGps:35' tmp13 = 2*q0;
  tmp13 = 2.0F * states[0];

  // 'computeMagMeasJacNoGps:36' tmp14 = q3*tmp13;
  tmp14 = states[3] * tmp13;

  // 'computeMagMeasJacNoGps:37' tmp15 = 2*q1;
  tmp15 = 2.0F * states[1];

  // 'computeMagMeasJacNoGps:38' tmp16 = q2*tmp13;
  tmp16 = states[2] * tmp13;

  // 'computeMagMeasJacNoGps:39' tmp17 = 2*magN;
  tmp17 = 2.0F * states[11];

  // 'computeMagMeasJacNoGps:40' tmp18 = -q3*tmp17;
  tmp18 = -states[3] * tmp17;

  // 'computeMagMeasJacNoGps:41' tmp19 = 4*magE;
  tmp19 = 4.0F * states[12];

  // 'computeMagMeasJacNoGps:42' tmp20 = magN*tmp15;
  tmp20 = states[11] * tmp15;

  // 'computeMagMeasJacNoGps:43' tmp21 = magN*tmp13;
  tmp21 = states[11] * tmp13;

  // 'computeMagMeasJacNoGps:44' tmp22 = 2*q1^2;
  tmp22 = states[1] * states[1] * 2.0F;

  // 'computeMagMeasJacNoGps:45' tmp23 = q1*tmp13;
  tmp13 *= states[1];

  // 'computeMagMeasJacNoGps:46' tmp24 = 2*q2;
  // 'computeMagMeasJacNoGps:47' tmp25 = 4*magD;
  tmp25 = 4.0F * states[13];

  // 'computeMagMeasJacNoGps:49' magMeasJac(1, 1) = 2*magE*q3 - tmp2;
  magMeasJac[0] = 2.0F * states[12] * states[3] - tmp2;

  // 'computeMagMeasJacNoGps:50' magMeasJac(1, 2) = tmp3 + tmp5;
  magMeasJac[3] = tmp3 + tmp5;

  // 'computeMagMeasJacNoGps:51' magMeasJac(1, 3) = -q2*tmp8 - tmp6 + tmp7;
  magMeasJac[6] = (-states[2] * tmp8 - tmp6) + tmp7;

  // 'computeMagMeasJacNoGps:52' magMeasJac(1, 4) = -q3*tmp8 + tmp10 + tmp9;
  magMeasJac[9] = (-states[3] * tmp8 + tmp10) + tmp1;

  // 'computeMagMeasJacNoGps:53' magMeasJac(1, 12) = -tmp11 - tmp12;
  magMeasJac[33] = -tmp11 - tmp12;

  // 'computeMagMeasJacNoGps:54' magMeasJac(1, 13) = q2*tmp15 + tmp14;
  magMeasJac[36] = states[2] * tmp15 + tmp14;

  // 'computeMagMeasJacNoGps:55' magMeasJac(1, 14) = 2*q1*q3 - tmp16;
  magMeasJac[39] = 2.0F * states[1] * states[3] - tmp16;

  // 'computeMagMeasJacNoGps:56' magMeasJac(1, 15) = 1;
  magMeasJac[42] = 1.0F;

  // 'computeMagMeasJacNoGps:58' magMeasJac(2, 1) = tmp18 + tmp9;
  magMeasJac[1] = tmp18 + tmp1;

  // 'computeMagMeasJacNoGps:59' magMeasJac(2, 2) = -q1*tmp19 + q2*tmp17 + tmp6; 
  magMeasJac[4] = (-states[1] * tmp19 + states[2] * tmp17) + tmp6;

  // 'computeMagMeasJacNoGps:60' magMeasJac(2, 3) = tmp20 + tmp3;
  magMeasJac[7] = tmp20 + tmp3;

  // 'computeMagMeasJacNoGps:61' magMeasJac(2, 4) = -q3*tmp19 + tmp2 - tmp21;
  magMeasJac[10] = (-states[3] * tmp19 + tmp2) - tmp21;

  // 'computeMagMeasJacNoGps:62' magMeasJac(2, 12) = 2*q1*q2 - tmp14;
  magMeasJac[34] = 2.0F * states[1] * states[2] - tmp14;

  // 'computeMagMeasJacNoGps:63' magMeasJac(2, 13) = -tmp12 - tmp22;
  magMeasJac[37] = -tmp12 - tmp22;

  // 'computeMagMeasJacNoGps:64' magMeasJac(2, 14) = q3*tmp24 + tmp23;
  tmp1 = 2.0F * states[2] * states[3];
  magMeasJac[40] = tmp1 + tmp13;

  // 'computeMagMeasJacNoGps:65' magMeasJac(2, 16) = 1;
  magMeasJac[46] = 1.0F;

  // 'computeMagMeasJacNoGps:67' magMeasJac(3, 1) = 2*magN*q2 - tmp7;
  magMeasJac[2] = 2.0F * states[11] * states[2] - tmp7;

  // 'computeMagMeasJacNoGps:68' magMeasJac(3, 2) = -q1*tmp25 - tmp10 - tmp18;
  magMeasJac[5] = (-states[1] * tmp25 - tmp10) - tmp18;

  // 'computeMagMeasJacNoGps:69' magMeasJac(3, 3) = -q2*tmp25 + q3*tmp4 + tmp21; 
  magMeasJac[8] = (-states[2] * tmp25 + states[3] * tmp4) + tmp21;

  // 'computeMagMeasJacNoGps:70' magMeasJac(3, 4) = tmp20 + tmp5;
  magMeasJac[11] = tmp20 + tmp5;

  // 'computeMagMeasJacNoGps:71' magMeasJac(3, 12) = q3*tmp15 + tmp16;
  magMeasJac[35] = states[3] * tmp15 + tmp16;

  // 'computeMagMeasJacNoGps:72' magMeasJac(3, 13) = 2*q2*q3 - tmp23;
  magMeasJac[38] = tmp1 - tmp13;

  // 'computeMagMeasJacNoGps:73' magMeasJac(3, 14) = -tmp11 - tmp22 + 1;
  magMeasJac[41] = (-tmp11 - tmp22) + 1.0F;

  // 'computeMagMeasJacNoGps:74' magMeasJac(3, 17) = 1;
  magMeasJac[50] = 1.0F;
}

//
// File trailer for generated code.
//
// [EOF]
//
