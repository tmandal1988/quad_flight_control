//
// File: computeStateJac_YJuzBPAK.cpp
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
#include "rtwtypes.h"
#include "computeStateJac_YJuzBPAK.h"

//
// Function for MATLAB Function: '<S1>/EKF'
// function  stateJac = computeStateJac(states, bodyAccels_mps2, bodyRates_radps, sampleTime_s)
// COMPUTESTATEJAC Computes the EKF state Jacobian
//
// Inputs:
// states:             Time propogated states
// bodyAccels_mps2:    Measured body accels
// bodyRates_radps:    Measured body rates
// sampleTime_s:       Sample time
//
// Outputs:
// stateJac:           State Jacobian
//
void computeStateJac_YJuzBPAK(const real32_T states[23], const real32_T
  bodyAccels_mps2[3], const real32_T bodyRates_radps[3], real32_T b_sampleTime_s,
  real32_T stateJac[58])
{
  real32_T t1;
  real32_T t2;
  real32_T t3;
  real32_T tmp1;
  real32_T tmp10;
  real32_T tmp12;
  real32_T tmp13;
  real32_T tmp15;
  real32_T tmp16;
  real32_T tmp17;
  real32_T tmp18;
  real32_T tmp19;
  real32_T tmp2;
  real32_T tmp20;
  real32_T tmp21;
  real32_T tmp22;
  real32_T tmp23;
  real32_T tmp24;
  real32_T tmp25;
  real32_T tmp28;
  real32_T tmp29;
  real32_T tmp3;
  real32_T tmp30;
  real32_T tmp31;
  real32_T tmp32;
  real32_T tmp35;
  real32_T tmp36;
  real32_T tmp37;
  real32_T tmp39;
  real32_T tmp4;
  real32_T tmp40;
  real32_T tmp5;
  real32_T tmp6;
  real32_T tmp7;

  // Initialize the state jacobian to zero
  //  stateJac = zeros(23, 23, 'single');
  // 'computeStateJac:15' stateJac = zeros(58, 1, 'single');
  // Extract quat states
  // 'computeStateJac:18' q0 = states(1);
  // 'computeStateJac:19' q1 = states(2);
  // 'computeStateJac:20' q2 = states(3);
  // 'computeStateJac:21' q3 = states(4);
  // Extract gyro biases from states
  // 'computeStateJac:23' bwx = states(11);
  // 'computeStateJac:24' bwy = states(12);
  // 'computeStateJac:25' bwz = states(13);
  // Extract accel biases from states
  // 'computeStateJac:27' bax = states(14);
  // 'computeStateJac:28' bay = states(15);
  // 'computeStateJac:29' baz = states(16);
  // Extract gyro inputs
  // 'computeStateJac:31' wx = bodyRates_radps(1);
  // 'computeStateJac:32' wy = bodyRates_radps(2);
  // 'computeStateJac:33' wz = bodyRates_radps(3);
  // Extract accel inputs
  // 'computeStateJac:35' ax = bodyAccels_mps2(1);
  // 'computeStateJac:36' ay = bodyAccels_mps2(2);
  // 'computeStateJac:37' az = bodyAccels_mps2(3);
  // 'computeStateJac:39' t1 = (bwx - wx);
  t1 = states[10] - bodyRates_radps[0];

  // 'computeStateJac:40' t2 = (bwy - wy);
  t2 = states[11] - bodyRates_radps[1];

  // 'computeStateJac:41' t3 = (bwz - wz);
  t3 = states[12] - bodyRates_radps[2];

  // 'computeStateJac:42' tmp1 = sampleTime_s*t1*0.5;
  tmp1 = b_sampleTime_s * t1 * 0.5F;

  // 'computeStateJac:43' tmp2 = sampleTime_s*t2*0.5;
  tmp2 = b_sampleTime_s * t2 * 0.5F;

  // 'computeStateJac:44' tmp3 = sampleTime_s*t3*0.5;
  tmp3 = b_sampleTime_s * t3 * 0.5F;

  // 'computeStateJac:45' tmp4 = sampleTime_s/2;
  tmp4 = b_sampleTime_s / 2.0F;

  // 'computeStateJac:46' tmp5 = q1*tmp4;
  tmp5 = states[1] * tmp4;

  // 'computeStateJac:47' tmp6 = q2*tmp4;
  tmp6 = states[2] * tmp4;

  // 'computeStateJac:48' tmp7 = q3*tmp4;
  tmp7 = states[3] * tmp4;

  // 'computeStateJac:49' tmp8 = -tmp4*t1;
  t1 *= -tmp4;

  // 'computeStateJac:50' tmp9 = -tmp4*t3;
  t3 *= -tmp4;

  // 'computeStateJac:51' tmp10 = -q0*tmp4;
  tmp10 = -states[0] * tmp4;

  // 'computeStateJac:52' tmp11 = -tmp4*t2;
  t2 *= -tmp4;

  // 'computeStateJac:53' tmp12 = az - baz;
  tmp12 = bodyAccels_mps2[2] - states[15];

  // 'computeStateJac:54' tmp13 = 2*q2;
  tmp13 = 2.0F * states[2];

  // 'computeStateJac:55' tmp14 = tmp12*tmp13;
  tmp4 = tmp12 * tmp13;

  // 'computeStateJac:56' tmp15 = ay - bay;
  tmp15 = bodyAccels_mps2[1] - states[14];

  // 'computeStateJac:57' tmp16 = 2*q3;
  tmp16 = 2.0F * states[3];

  // 'computeStateJac:58' tmp17 = -tmp15*tmp16;
  tmp17 = -tmp15 * tmp16;

  // 'computeStateJac:59' tmp18 = tmp13*tmp15;
  tmp18 = tmp13 * tmp15;

  // 'computeStateJac:60' tmp19 = tmp12*tmp16;
  tmp19 = tmp12 * tmp16;

  // 'computeStateJac:61' tmp20 = 2*tmp12;
  tmp20 = 2.0F * tmp12;

  // 'computeStateJac:62' tmp21 = q0*tmp20;
  tmp21 = states[0] * tmp20;

  // 'computeStateJac:63' tmp22 = 2*tmp15;
  tmp22 = 2.0F * tmp15;

  // 'computeStateJac:64' tmp23 = q1*tmp22;
  tmp23 = states[1] * tmp22;

  // 'computeStateJac:65' tmp24 = ax - bax;
  tmp24 = bodyAccels_mps2[0] - states[13];

  // 'computeStateJac:66' tmp25 = 4*tmp24;
  tmp25 = 4.0F * tmp24;

  // 'computeStateJac:67' tmp26 = q0*tmp22;
  tmp22 *= states[0];

  // 'computeStateJac:68' tmp27 = q1*tmp20;
  tmp20 *= states[1];

  // 'computeStateJac:69' tmp28 = 2*q2*q2;
  tmp28 = 2.0F * states[2] * states[2];

  // 'computeStateJac:70' tmp29 = 2*q3*q3 - 1;
  tmp29 = 2.0F * states[3] * states[3] - 1.0F;

  // 'computeStateJac:71' tmp30 = q0*tmp16;
  tmp30 = states[0] * tmp16;

  // 'computeStateJac:72' tmp31 = q1*tmp13;
  tmp31 = states[1] * tmp13;

  // 'computeStateJac:73' tmp32 = q0*tmp13;
  tmp32 = states[0] * tmp13;

  // 'computeStateJac:74' tmp33 = q1*tmp16;
  tmp16 *= states[1];

  // 'computeStateJac:75' tmp34 = 4*tmp15;
  tmp15 *= 4.0F;

  // 'computeStateJac:76' tmp35 = -tmp13*tmp24;
  tmp35 = -tmp13 * tmp24;

  // 'computeStateJac:77' tmp36 = 2*tmp24;
  tmp36 = 2.0F * tmp24;

  // 'computeStateJac:78' tmp37 = q1*tmp36;
  tmp37 = states[1] * tmp36;

  // 'computeStateJac:79' tmp38 = q0*tmp36;
  tmp36 *= states[0];

  // 'computeStateJac:80' tmp39 = 2*q1*q1;
  tmp39 = 2.0F * states[1] * states[1];

  // 'computeStateJac:81' tmp40 = 2*q0*q1;
  tmp40 = 2.0F * states[0] * states[1];

  // 'computeStateJac:82' tmp41 = q3*tmp13;
  tmp13 *= states[3];

  // 'computeStateJac:83' tmp42 = 4*tmp12;
  tmp12 *= 4.0F;

  // stateJac(1, 1) = 1;
  // 'computeStateJac:86' stateJac(1) = 1;
  stateJac[0] = 1.0F;

  // stateJac(1, 2) = tmp1;
  // 'computeStateJac:89' stateJac(2) = tmp1;
  stateJac[1] = tmp1;

  // stateJac(1, 3) = tmp2;
  // 'computeStateJac:92' stateJac(3) = tmp2;
  stateJac[2] = tmp2;

  // stateJac(1, 4) = tmp3;
  // 'computeStateJac:95' stateJac(4) = tmp3;
  stateJac[3] = tmp3;

  // stateJac(1, 11) = tmp5;
  // 'computeStateJac:98' stateJac(5) = tmp5;
  stateJac[4] = tmp5;

  // stateJac(1, 12) = tmp6;
  // 'computeStateJac:101' stateJac(6) = tmp6;
  stateJac[5] = tmp6;

  // stateJac(1, 13) = tmp7;
  // 'computeStateJac:104' stateJac(7) = tmp7;
  stateJac[6] = tmp7;

  // stateJac(2, 1) = tmp8;
  // 'computeStateJac:107' stateJac(8) = tmp8;
  stateJac[7] = t1;

  // stateJac(2, 2) = 1;
  // 'computeStateJac:110' stateJac(9) = 1;
  stateJac[8] = 1.0F;

  // stateJac(2, 3) = tmp9;
  // 'computeStateJac:113' stateJac(10) = tmp9;
  stateJac[9] = t3;

  // stateJac(2, 4) = tmp2;
  // 'computeStateJac:116' stateJac(11) = tmp2;
  stateJac[10] = tmp2;

  // stateJac(2, 11) = tmp10;
  // 'computeStateJac:119' stateJac(12) = tmp10;
  stateJac[11] = tmp10;

  // stateJac(2, 12) = tmp7;
  // 'computeStateJac:122' stateJac(13) = tmp7;
  stateJac[12] = tmp7;

  // stateJac(2, 13) = -tmp6;
  // 'computeStateJac:125' stateJac(14) = -tmp6;
  stateJac[13] = -tmp6;

  // stateJac(3, 1) = tmp11;
  // 'computeStateJac:128' stateJac(15) = tmp11;
  stateJac[14] = t2;

  // stateJac(3, 2) = tmp3;
  // 'computeStateJac:131' stateJac(16) = tmp3;
  stateJac[15] = tmp3;

  // stateJac(3, 3) = 1;
  // 'computeStateJac:134' stateJac(17) = 1;
  stateJac[16] = 1.0F;

  // stateJac(3, 4) = tmp8;
  // 'computeStateJac:137' stateJac(18) = tmp8;
  stateJac[17] = t1;

  // stateJac(3, 11) = -tmp7;
  // 'computeStateJac:140' stateJac(19) = -tmp7;
  stateJac[18] = -tmp7;

  // stateJac(3, 12) = tmp10;
  // 'computeStateJac:143' stateJac(20) = tmp10;
  stateJac[19] = tmp10;

  // stateJac(3, 13) = tmp5;
  // 'computeStateJac:146' stateJac(21) = tmp5;
  stateJac[20] = tmp5;

  // stateJac(4, 1) = tmp9;
  // 'computeStateJac:149' stateJac(22) = tmp9;
  stateJac[21] = t3;

  // stateJac(4, 2) = tmp11;
  // 'computeStateJac:152' stateJac(23) = tmp11;
  stateJac[22] = t2;

  // stateJac(4, 3) = tmp1;
  // 'computeStateJac:155' stateJac(24) = tmp1;
  stateJac[23] = tmp1;

  // stateJac(4, 4) = 1;
  // 'computeStateJac:158' stateJac(25) = 1;
  stateJac[24] = 1.0F;

  // stateJac(4, 11) = tmp6;
  // 'computeStateJac:161' stateJac(26) = tmp6;
  stateJac[25] = tmp6;

  // stateJac(4, 12) = -tmp5;
  // 'computeStateJac:164' stateJac(27) = -tmp5;
  stateJac[26] = -tmp5;

  // stateJac(4, 13) = tmp10;
  // 'computeStateJac:167' stateJac(28) = tmp10;
  stateJac[27] = tmp10;

  // stateJac(5, 5) = 1;
  // 'computeStateJac:170' stateJac(29) = 1;
  stateJac[28] = 1.0F;

  // stateJac(5, 8) = sampleTime_s;
  // 'computeStateJac:173' stateJac(30) = sampleTime_s;
  stateJac[29] = b_sampleTime_s;

  // stateJac(6, 6) = 1;
  // 'computeStateJac:176' stateJac(31) = 1;
  stateJac[30] = 1.0F;

  // stateJac(6, 9) = sampleTime_s;
  // 'computeStateJac:179' stateJac(32) = sampleTime_s;
  stateJac[31] = b_sampleTime_s;

  // stateJac(7, 7) = 1;
  // 'computeStateJac:182' stateJac(33) = 1;
  stateJac[32] = 1.0F;

  // stateJac(7, 10) = sampleTime_s;
  // 'computeStateJac:185' stateJac(34) = sampleTime_s;
  stateJac[33] = b_sampleTime_s;

  // stateJac(8, 1) = -sampleTime_s*(-tmp14 - tmp17);
  // 'computeStateJac:188' stateJac(35) = -sampleTime_s*(-tmp14 - tmp17);
  stateJac[34] = (-tmp4 - tmp17) * -b_sampleTime_s;

  // stateJac(8, 2) = sampleTime_s*(tmp18 + tmp19);
  // 'computeStateJac:191' stateJac(36) = sampleTime_s*(tmp18 + tmp19);
  stateJac[35] = (tmp18 + tmp19) * b_sampleTime_s;

  // stateJac(8, 3) = sampleTime_s*(-q2*tmp25 + tmp21 + tmp23);
  // 'computeStateJac:194' stateJac(37) = sampleTime_s*(-q2*tmp25 + tmp21 + tmp23); 
  stateJac[36] = ((-states[2] * tmp25 + tmp21) + tmp23) * b_sampleTime_s;

  // stateJac(8, 4) = -sampleTime_s*(q3*tmp25 + tmp26 - tmp27);
  // 'computeStateJac:197' stateJac(38) = -sampleTime_s*(q3*tmp25 + tmp26 - tmp27); 
  stateJac[37] = ((states[3] * tmp25 + tmp22) - tmp20) * -b_sampleTime_s;

  // stateJac(8, 8) = 1;
  // 'computeStateJac:200' stateJac(39) = 1;
  stateJac[38] = 1.0F;

  // stateJac(8, 14) = sampleTime_s*(tmp28 + tmp29);
  // 'computeStateJac:203' stateJac(40) = sampleTime_s*(tmp28 + tmp29);
  stateJac[39] = (tmp28 + tmp29) * b_sampleTime_s;

  // stateJac(8, 15) = sampleTime_s*(tmp30 - tmp31);
  // 'computeStateJac:206' stateJac(41) = sampleTime_s*(tmp30 - tmp31);
  stateJac[40] = (tmp30 - tmp31) * b_sampleTime_s;

  // stateJac(8, 16) = -sampleTime_s*(tmp32 + tmp33);
  // 'computeStateJac:209' stateJac(42) = -sampleTime_s*(tmp32 + tmp33);
  stateJac[41] = (tmp32 + tmp16) * -b_sampleTime_s;

  // stateJac(9, 1) = sampleTime_s*(2*q3*tmp24 - tmp27);
  // 'computeStateJac:212' stateJac(43) = sampleTime_s*(2*q3*tmp24 - tmp27);
  t1 = 2.0F * states[3] * tmp24;
  stateJac[42] = (t1 - tmp20) * b_sampleTime_s;

  // stateJac(9, 2) = -sampleTime_s*(q1*tmp34 + tmp21 + tmp35);
  // 'computeStateJac:215' stateJac(44) = -sampleTime_s*(q1*tmp34 + tmp21 + tmp35); 
  stateJac[43] = ((states[1] * tmp15 + tmp21) + tmp35) * -b_sampleTime_s;

  // stateJac(9, 3) = sampleTime_s*(tmp19 + tmp37);
  // 'computeStateJac:218' stateJac(45) = sampleTime_s*(tmp19 + tmp37);
  stateJac[44] = (tmp19 + tmp37) * b_sampleTime_s;

  // stateJac(9, 4) = sampleTime_s*(-q3*tmp34 + tmp14 + tmp38);
  // 'computeStateJac:221' stateJac(46) = sampleTime_s*(-q3*tmp34 + tmp14 + tmp38); 
  stateJac[45] = ((-states[3] * tmp15 + tmp4) + tmp36) * b_sampleTime_s;

  // stateJac(9, 9) = 1;
  // 'computeStateJac:224' stateJac(47) = 1;
  stateJac[46] = 1.0F;

  // stateJac(9, 14) = -sampleTime_s*(tmp30 + tmp31);
  // 'computeStateJac:227' stateJac(48) = -sampleTime_s*(tmp30 + tmp31);
  stateJac[47] = (tmp30 + tmp31) * -b_sampleTime_s;

  // stateJac(9, 15) = sampleTime_s*(tmp29 + tmp39);
  // 'computeStateJac:230' stateJac(49) = sampleTime_s*(tmp29 + tmp39);
  stateJac[48] = (tmp29 + tmp39) * b_sampleTime_s;

  // stateJac(9, 16) = sampleTime_s*(tmp40 - tmp41);
  // 'computeStateJac:233' stateJac(50) = sampleTime_s*(tmp40 - tmp41);
  stateJac[49] = (tmp40 - tmp13) * b_sampleTime_s;

  // stateJac(10, 1) = -sampleTime_s*(-tmp23 - tmp35);
  // 'computeStateJac:236' stateJac(51) = -sampleTime_s*(-tmp23 - tmp35);
  stateJac[50] = (-tmp23 - tmp35) * -b_sampleTime_s;

  // stateJac(10, 2) = sampleTime_s*(-q1*tmp42 + tmp16*tmp24 + tmp26);
  // 'computeStateJac:239' stateJac(52) = sampleTime_s*(-q1*tmp42 + tmp16*tmp24 + tmp26); 
  stateJac[51] = ((-states[1] * tmp12 + t1) + tmp22) * b_sampleTime_s;

  // stateJac(10, 3) = -sampleTime_s*(q2*tmp42 + tmp17 + tmp38);
  // 'computeStateJac:242' stateJac(53) = -sampleTime_s*(q2*tmp42 + tmp17 + tmp38); 
  stateJac[52] = ((states[2] * tmp12 + tmp17) + tmp36) * -b_sampleTime_s;

  // stateJac(10, 4) = sampleTime_s*(tmp18 + tmp37);
  // 'computeStateJac:245' stateJac(54) = sampleTime_s*(tmp18 + tmp37);
  stateJac[53] = (tmp18 + tmp37) * b_sampleTime_s;

  // stateJac(10, 10) = 1;
  // 'computeStateJac:248' stateJac(55) = 1;
  stateJac[54] = 1.0F;

  // stateJac(10, 14) = sampleTime_s*(tmp32 - tmp33);
  // 'computeStateJac:251' stateJac(56) = sampleTime_s*(tmp32 - tmp33);
  stateJac[55] = (tmp32 - tmp16) * b_sampleTime_s;

  // stateJac(10, 15) = -sampleTime_s*(tmp40 + tmp41);
  // 'computeStateJac:254' stateJac(57) = -sampleTime_s*(tmp40 + tmp41);
  stateJac[56] = (tmp40 + tmp13) * -b_sampleTime_s;

  // stateJac(10, 16) = sampleTime_s*(tmp28 + tmp39 - 1);
  // 'computeStateJac:257' stateJac(58) = sampleTime_s*(tmp28 + tmp39 - 1);
  stateJac[57] = ((tmp28 + tmp39) - 1.0F) * b_sampleTime_s;

  //  %stateJac(11, 11) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(12, 12) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(13, 13) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(14, 14) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(15, 15) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(16, 16) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(17, 17) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(18, 18) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(19, 19) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(20, 20) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(21, 21) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(22, 22) = 1;
  //  stateJac(idx) = 1;
  //  idx = idx + 1;
  //
  //  %stateJac(23, 23) = 1;
  //  stateJac(idx) = 1;
}

//
// File trailer for generated code.
//
// [EOF]
//
