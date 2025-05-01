//
// File: updateCovPNoGps_vX01OT1j.cpp
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
#include "rtwtypes.h"
#include "updateCovPNoGps_vX01OT1j.h"

//
// Function for MATLAB Function: '<S1>/EKF NO GPS'
// function  covP = updateCovPNoGps(covP, stateJac, processNoiseQ)
// UPDATECOVTEST Computes the covariance for the prediction stage of the EKF
//
// Inputs:
// covP:           Covariance from previous time step
// stateJac:           State Jacobian using previous time estimate
// processNoiseQ:      Process noise Q
//
// Outputs:
// covP:               Updated covariance
//
void updateCovPNoGps_vX01OT1j(real32_T covP[324], const real32_T stateJac[28],
  const real32_T processNoiseQ[324])
{
  real32_T tmp1;
  real32_T tmp10;
  real32_T tmp106;
  real32_T tmp107;
  real32_T tmp108;
  real32_T tmp109;
  real32_T tmp11;
  real32_T tmp110;
  real32_T tmp111;
  real32_T tmp112;
  real32_T tmp12;
  real32_T tmp122;
  real32_T tmp123;
  real32_T tmp124;
  real32_T tmp125;
  real32_T tmp126;
  real32_T tmp127;
  real32_T tmp128;
  real32_T tmp13;
  real32_T tmp137;
  real32_T tmp138;
  real32_T tmp139;
  real32_T tmp14;
  real32_T tmp140;
  real32_T tmp141;
  real32_T tmp142;
  real32_T tmp143;
  real32_T tmp15;
  real32_T tmp151;
  real32_T tmp152;
  real32_T tmp153;
  real32_T tmp154;
  real32_T tmp155;
  real32_T tmp156;
  real32_T tmp157;
  real32_T tmp16;
  real32_T tmp164;
  real32_T tmp165;
  real32_T tmp166;
  real32_T tmp167;
  real32_T tmp168;
  real32_T tmp169;
  real32_T tmp17;
  real32_T tmp170;
  real32_T tmp176;
  real32_T tmp177;
  real32_T tmp178;
  real32_T tmp179;
  real32_T tmp18;
  real32_T tmp180;
  real32_T tmp181;
  real32_T tmp182;
  real32_T tmp187;
  real32_T tmp188;
  real32_T tmp189;
  real32_T tmp19;
  real32_T tmp190;
  real32_T tmp191;
  real32_T tmp192;
  real32_T tmp193;
  real32_T tmp197;
  real32_T tmp198;
  real32_T tmp199;
  real32_T tmp2;
  real32_T tmp20;
  real32_T tmp200;
  real32_T tmp201;
  real32_T tmp202;
  real32_T tmp203;
  real32_T tmp206;
  real32_T tmp207;
  real32_T tmp208;
  real32_T tmp209;
  real32_T tmp21;
  real32_T tmp210;
  real32_T tmp211;
  real32_T tmp212;
  real32_T tmp214;
  real32_T tmp215;
  real32_T tmp216;
  real32_T tmp217;
  real32_T tmp218;
  real32_T tmp219;
  real32_T tmp22;
  real32_T tmp220;
  real32_T tmp23;
  real32_T tmp24;
  real32_T tmp25;
  real32_T tmp26;
  real32_T tmp27;
  real32_T tmp28;
  real32_T tmp29;
  real32_T tmp3;
  real32_T tmp30;
  real32_T tmp31;
  real32_T tmp32;
  real32_T tmp33;
  real32_T tmp34;
  real32_T tmp35;
  real32_T tmp36;
  real32_T tmp37;
  real32_T tmp38;
  real32_T tmp39;
  real32_T tmp4;
  real32_T tmp40;
  real32_T tmp41;
  real32_T tmp42;
  real32_T tmp43;
  real32_T tmp44;
  real32_T tmp45;
  real32_T tmp46;
  real32_T tmp47;
  real32_T tmp5;
  real32_T tmp58;
  real32_T tmp59;
  real32_T tmp6;
  real32_T tmp60;
  real32_T tmp61;
  real32_T tmp62;
  real32_T tmp63;
  real32_T tmp7;
  real32_T tmp74;
  real32_T tmp75;
  real32_T tmp76;
  real32_T tmp77;
  real32_T tmp78;
  real32_T tmp79;
  real32_T tmp8;
  real32_T tmp9;
  real32_T tmp90;
  real32_T tmp91;
  real32_T tmp92;
  real32_T tmp93;
  real32_T tmp94;
  real32_T tmp95;

  //  sJ1_1 = stateJac_(1, 1)
  // 'updateCovPNoGps:13' sJ1_1 = stateJac(1);
  //  sJ1_2 = stateJac_(1, 2)
  // 'updateCovPNoGps:16' sJ1_2 = stateJac(2);
  //  sJ1_3 = stateJac_(1, 3)
  // 'updateCovPNoGps:19' sJ1_3 = stateJac(3);
  //  sJ1_4 = stateJac_(1, 4)
  // 'updateCovPNoGps:22' sJ1_4 = stateJac(4);
  //  sJ1_6 = stateJac_(1, 6)
  // 'updateCovPNoGps:25' sJ1_6 = stateJac(5);
  //  sJ1_7 = stateJac_(1, 7)
  // 'updateCovPNoGps:28' sJ1_7 = stateJac(6);
  //  sJ1_8 = stateJac_(1, 8)
  // 'updateCovPNoGps:31' sJ1_8 = stateJac(7);
  //  sJ2_1 = stateJac_(2, 1)
  // 'updateCovPNoGps:34' sJ2_1 = stateJac(8);
  //  sJ2_2 = stateJac_(2, 2)
  // 'updateCovPNoGps:37' sJ2_2 = stateJac(9);
  //  sJ2_3 = stateJac_(2, 3)
  // 'updateCovPNoGps:40' sJ2_3 = stateJac(10);
  //  sJ2_4 = stateJac_(2, 4)
  // 'updateCovPNoGps:43' sJ2_4 = stateJac(11);
  //  sJ2_6 = stateJac_(2, 6)
  // 'updateCovPNoGps:46' sJ2_6 = stateJac(12);
  //  sJ2_7 = stateJac_(2, 7)
  // 'updateCovPNoGps:49' sJ2_7 = stateJac(13);
  //  sJ2_8 = stateJac_(2, 8)
  // 'updateCovPNoGps:52' sJ2_8 = stateJac(14);
  //  sJ3_1 = stateJac_(3, 1)
  // 'updateCovPNoGps:55' sJ3_1 = stateJac(15);
  //  sJ3_2 = stateJac_(3, 2)
  // 'updateCovPNoGps:58' sJ3_2 = stateJac(16);
  //  sJ3_3 = stateJac_(3, 3)
  // 'updateCovPNoGps:61' sJ3_3 = stateJac(17);
  //  sJ3_4 = stateJac_(3, 4)
  // 'updateCovPNoGps:64' sJ3_4 = stateJac(18);
  //  sJ3_6 = stateJac_(3, 6)
  // 'updateCovPNoGps:67' sJ3_6 = stateJac(19);
  //  sJ3_7 = stateJac_(3, 7)
  // 'updateCovPNoGps:70' sJ3_7 = stateJac(20);
  //  sJ3_8 = stateJac_(3, 8)
  // 'updateCovPNoGps:73' sJ3_8 = stateJac(21);
  //  sJ4_1 = stateJac_(4, 1)
  // 'updateCovPNoGps:76' sJ4_1 = stateJac(22);
  //  sJ4_2 = stateJac_(4, 2)
  // 'updateCovPNoGps:79' sJ4_2 = stateJac(23);
  //  sJ4_3 = stateJac_(4, 3)
  // 'updateCovPNoGps:82' sJ4_3 = stateJac(24);
  //  sJ4_4 = stateJac_(4, 4)
  // 'updateCovPNoGps:85' sJ4_4 = stateJac(25);
  //  sJ4_6 = stateJac_(4, 6)
  // 'updateCovPNoGps:88' sJ4_6 = stateJac(26);
  //  sJ4_7 = stateJac_(4, 7)
  // 'updateCovPNoGps:91' sJ4_7 = stateJac(27);
  //  sJ4_8 = stateJac_(4, 8)
  // 'updateCovPNoGps:94' sJ4_8 = stateJac(28);
  // 'updateCovPNoGps:96' tmp1 = covP(1, 1)*sJ1_1 + covP(2, 1)*sJ1_2 + covP(3, 1)*sJ1_3 + covP(4, 1)*sJ1_4 + covP(6, 1)*sJ1_6 + covP(7, 1)*sJ1_7 + covP(8, 1)*sJ1_8; 
  tmp1 = (((((covP[0] * stateJac[0] + covP[1] * stateJac[1]) + covP[2] *
             stateJac[2]) + covP[3] * stateJac[3]) + stateJac[4] * covP[5]) +
          stateJac[5] * covP[6]) + stateJac[6] * covP[7];

  // 'updateCovPNoGps:97' tmp2 = covP(1, 2)*sJ1_1 + covP(2, 2)*sJ1_2 + covP(3, 2)*sJ1_3 + covP(4, 2)*sJ1_4 + covP(6, 2)*sJ1_6 + covP(7, 2)*sJ1_7 + covP(8, 2)*sJ1_8; 
  tmp2 = (((((stateJac[0] * covP[18] + stateJac[1] * covP[19]) + stateJac[2] *
             covP[20]) + stateJac[3] * covP[21]) + stateJac[4] * covP[23]) +
          stateJac[5] * covP[24]) + stateJac[6] * covP[25];

  // 'updateCovPNoGps:98' tmp3 = covP(1, 3)*sJ1_1 + covP(2, 3)*sJ1_2 + covP(3, 3)*sJ1_3 + covP(4, 3)*sJ1_4 + covP(6, 3)*sJ1_6 + covP(7, 3)*sJ1_7 + covP(8, 3)*sJ1_8; 
  tmp3 = (((((stateJac[0] * covP[36] + stateJac[1] * covP[37]) + stateJac[2] *
             covP[38]) + stateJac[3] * covP[39]) + stateJac[4] * covP[41]) +
          stateJac[5] * covP[42]) + stateJac[6] * covP[43];

  // 'updateCovPNoGps:99' tmp4 = covP(1, 4)*sJ1_1 + covP(2, 4)*sJ1_2 + covP(3, 4)*sJ1_3 + covP(4, 4)*sJ1_4 + covP(6, 4)*sJ1_6 + covP(7, 4)*sJ1_7 + covP(8, 4)*sJ1_8; 
  tmp4 = (((((stateJac[0] * covP[54] + stateJac[1] * covP[55]) + stateJac[2] *
             covP[56]) + stateJac[3] * covP[57]) + stateJac[4] * covP[59]) +
          stateJac[5] * covP[60]) + stateJac[6] * covP[61];

  // 'updateCovPNoGps:100' tmp5 = covP(6, 6)*sJ1_6;
  tmp5 = stateJac[4] * covP[95];

  // 'updateCovPNoGps:101' tmp6 = covP(1, 6)*sJ1_1 + covP(2, 6)*sJ1_2 + covP(3, 6)*sJ1_3 + covP(4, 6)*sJ1_4 + covP(7, 6)*sJ1_7 + covP(8, 6)*sJ1_8 + tmp5; 
  tmp6 = (((((stateJac[0] * covP[90] + stateJac[1] * covP[91]) + stateJac[2] *
             covP[92]) + stateJac[3] * covP[93]) + stateJac[5] * covP[96]) +
          stateJac[6] * covP[97]) + tmp5;

  // 'updateCovPNoGps:102' tmp7 = covP(7, 7)*sJ1_7;
  tmp7 = stateJac[5] * covP[114];

  // 'updateCovPNoGps:103' tmp8 = covP(1, 7)*sJ1_1 + covP(2, 7)*sJ1_2 + covP(3, 7)*sJ1_3 + covP(4, 7)*sJ1_4 + covP(6, 7)*sJ1_6 + covP(8, 7)*sJ1_8 + tmp7; 
  tmp8 = (((((stateJac[0] * covP[108] + stateJac[1] * covP[109]) + stateJac[2] *
             covP[110]) + stateJac[3] * covP[111]) + stateJac[4] * covP[113]) +
          stateJac[6] * covP[115]) + tmp7;

  // 'updateCovPNoGps:104' tmp9 = covP(8, 8)*sJ1_8;
  tmp9 = stateJac[6] * covP[133];

  // 'updateCovPNoGps:105' tmp10 = covP(1, 8)*sJ1_1 + covP(2, 8)*sJ1_2 + covP(3, 8)*sJ1_3 + covP(4, 8)*sJ1_4 + covP(6, 8)*sJ1_6 + covP(7, 8)*sJ1_7 + tmp9; 
  tmp10 = (((((stateJac[0] * covP[126] + stateJac[1] * covP[127]) + stateJac[2] *
              covP[128]) + stateJac[3] * covP[129]) + stateJac[4] * covP[131]) +
           stateJac[5] * covP[132]) + tmp9;

  // 'updateCovPNoGps:106' tmp11 = covP(1, 1)*sJ2_1 + covP(2, 1)*sJ2_2 + covP(3, 1)*sJ2_3 + covP(4, 1)*sJ2_4 + covP(6, 1)*sJ2_6 + covP(7, 1)*sJ2_7 + covP(8, 1)*sJ2_8; 
  tmp11 = (((((covP[0] * stateJac[7] + covP[1] * stateJac[8]) + covP[2] *
              stateJac[9]) + covP[3] * stateJac[10]) + covP[5] * stateJac[11]) +
           covP[6] * stateJac[12]) + covP[7] * stateJac[13];

  // 'updateCovPNoGps:107' tmp12 = covP(1, 2)*sJ2_1 + covP(2, 2)*sJ2_2 + covP(3, 2)*sJ2_3 + covP(4, 2)*sJ2_4 + covP(6, 2)*sJ2_6 + covP(7, 2)*sJ2_7 + covP(8, 2)*sJ2_8; 
  tmp12 = (((((stateJac[7] * covP[18] + stateJac[8] * covP[19]) + stateJac[9] *
              covP[20]) + stateJac[10] * covP[21]) + stateJac[11] * covP[23]) +
           stateJac[12] * covP[24]) + stateJac[13] * covP[25];

  // 'updateCovPNoGps:108' tmp13 = covP(1, 3)*sJ2_1 + covP(2, 3)*sJ2_2 + covP(3, 3)*sJ2_3 + covP(4, 3)*sJ2_4 + covP(6, 3)*sJ2_6 + covP(7, 3)*sJ2_7 + covP(8, 3)*sJ2_8; 
  tmp13 = (((((stateJac[7] * covP[36] + stateJac[8] * covP[37]) + stateJac[9] *
              covP[38]) + stateJac[10] * covP[39]) + stateJac[11] * covP[41]) +
           stateJac[12] * covP[42]) + stateJac[13] * covP[43];

  // 'updateCovPNoGps:109' tmp14 = covP(1, 4)*sJ2_1 + covP(2, 4)*sJ2_2 + covP(3, 4)*sJ2_3 + covP(4, 4)*sJ2_4 + covP(6, 4)*sJ2_6 + covP(7, 4)*sJ2_7 + covP(8, 4)*sJ2_8; 
  tmp14 = (((((stateJac[7] * covP[54] + stateJac[8] * covP[55]) + stateJac[9] *
              covP[56]) + stateJac[10] * covP[57]) + stateJac[11] * covP[59]) +
           stateJac[12] * covP[60]) + stateJac[13] * covP[61];

  // 'updateCovPNoGps:110' tmp15 = covP(6, 6)*sJ2_6;
  tmp15 = stateJac[11] * covP[95];

  // 'updateCovPNoGps:111' tmp16 = covP(1, 6)*sJ2_1 + covP(2, 6)*sJ2_2 + covP(3, 6)*sJ2_3 + covP(4, 6)*sJ2_4 + covP(7, 6)*sJ2_7 + covP(8, 6)*sJ2_8 + tmp15; 
  tmp16 = (((((stateJac[7] * covP[90] + stateJac[8] * covP[91]) + stateJac[9] *
              covP[92]) + stateJac[10] * covP[93]) + stateJac[12] * covP[96]) +
           stateJac[13] * covP[97]) + tmp15;

  // 'updateCovPNoGps:112' tmp17 = covP(7, 7)*sJ2_7;
  tmp17 = stateJac[12] * covP[114];

  // 'updateCovPNoGps:113' tmp18 = covP(1, 7)*sJ2_1 + covP(2, 7)*sJ2_2 + covP(3, 7)*sJ2_3 + covP(4, 7)*sJ2_4 + covP(6, 7)*sJ2_6 + covP(8, 7)*sJ2_8 + tmp17; 
  tmp18 = (((((stateJac[7] * covP[108] + stateJac[8] * covP[109]) + stateJac[9] *
              covP[110]) + stateJac[10] * covP[111]) + stateJac[11] * covP[113])
           + stateJac[13] * covP[115]) + tmp17;

  // 'updateCovPNoGps:114' tmp19 = covP(8, 8)*sJ2_8;
  tmp19 = stateJac[13] * covP[133];

  // 'updateCovPNoGps:115' tmp20 = covP(1, 8)*sJ2_1 + covP(2, 8)*sJ2_2 + covP(3, 8)*sJ2_3 + covP(4, 8)*sJ2_4 + covP(6, 8)*sJ2_6 + covP(7, 8)*sJ2_7 + tmp19; 
  tmp20 = (((((stateJac[7] * covP[126] + stateJac[8] * covP[127]) + stateJac[9] *
              covP[128]) + stateJac[10] * covP[129]) + stateJac[11] * covP[131])
           + stateJac[12] * covP[132]) + tmp19;

  // 'updateCovPNoGps:116' tmp21 = covP(1, 1)*sJ3_1 + covP(2, 1)*sJ3_2 + covP(3, 1)*sJ3_3 + covP(4, 1)*sJ3_4 + covP(6, 1)*sJ3_6 + covP(7, 1)*sJ3_7 + covP(8, 1)*sJ3_8; 
  tmp21 = (((((covP[0] * stateJac[14] + covP[1] * stateJac[15]) + covP[2] *
              stateJac[16]) + covP[3] * stateJac[17]) + covP[5] * stateJac[18])
           + covP[6] * stateJac[19]) + covP[7] * stateJac[20];

  // 'updateCovPNoGps:117' tmp22 = covP(1, 2)*sJ3_1 + covP(2, 2)*sJ3_2 + covP(3, 2)*sJ3_3 + covP(4, 2)*sJ3_4 + covP(6, 2)*sJ3_6 + covP(7, 2)*sJ3_7 + covP(8, 2)*sJ3_8; 
  tmp22 = (((((stateJac[14] * covP[18] + stateJac[15] * covP[19]) + stateJac[16]
              * covP[20]) + stateJac[17] * covP[21]) + stateJac[18] * covP[23])
           + stateJac[19] * covP[24]) + stateJac[20] * covP[25];

  // 'updateCovPNoGps:118' tmp23 = covP(1, 3)*sJ3_1 + covP(2, 3)*sJ3_2 + covP(3, 3)*sJ3_3 + covP(4, 3)*sJ3_4 + covP(6, 3)*sJ3_6 + covP(7, 3)*sJ3_7 + covP(8, 3)*sJ3_8; 
  tmp23 = (((((stateJac[14] * covP[36] + stateJac[15] * covP[37]) + stateJac[16]
              * covP[38]) + stateJac[17] * covP[39]) + stateJac[18] * covP[41])
           + stateJac[19] * covP[42]) + stateJac[20] * covP[43];

  // 'updateCovPNoGps:119' tmp24 = covP(1, 4)*sJ3_1 + covP(2, 4)*sJ3_2 + covP(3, 4)*sJ3_3 + covP(4, 4)*sJ3_4 + covP(6, 4)*sJ3_6 + covP(7, 4)*sJ3_7 + covP(8, 4)*sJ3_8; 
  tmp24 = (((((stateJac[14] * covP[54] + stateJac[15] * covP[55]) + stateJac[16]
              * covP[56]) + stateJac[17] * covP[57]) + stateJac[18] * covP[59])
           + stateJac[19] * covP[60]) + stateJac[20] * covP[61];

  // 'updateCovPNoGps:120' tmp25 = covP(6, 6)*sJ3_6;
  tmp25 = stateJac[18] * covP[95];

  // 'updateCovPNoGps:121' tmp26 = covP(1, 6)*sJ3_1 + covP(2, 6)*sJ3_2 + covP(3, 6)*sJ3_3 + covP(4, 6)*sJ3_4 + covP(7, 6)*sJ3_7 + covP(8, 6)*sJ3_8 + tmp25; 
  tmp26 = (((((stateJac[14] * covP[90] + stateJac[15] * covP[91]) + stateJac[16]
              * covP[92]) + stateJac[17] * covP[93]) + stateJac[19] * covP[96])
           + stateJac[20] * covP[97]) + tmp25;

  // 'updateCovPNoGps:122' tmp27 = covP(7, 7)*sJ3_7;
  tmp27 = stateJac[19] * covP[114];

  // 'updateCovPNoGps:123' tmp28 = covP(1, 7)*sJ3_1 + covP(2, 7)*sJ3_2 + covP(3, 7)*sJ3_3 + covP(4, 7)*sJ3_4 + covP(6, 7)*sJ3_6 + covP(8, 7)*sJ3_8 + tmp27; 
  tmp28 = (((((stateJac[14] * covP[108] + stateJac[15] * covP[109]) + stateJac
              [16] * covP[110]) + stateJac[17] * covP[111]) + stateJac[18] *
            covP[113]) + stateJac[20] * covP[115]) + tmp27;

  // 'updateCovPNoGps:124' tmp29 = covP(8, 8)*sJ3_8;
  tmp29 = stateJac[20] * covP[133];

  // 'updateCovPNoGps:125' tmp30 = covP(1, 8)*sJ3_1 + covP(2, 8)*sJ3_2 + covP(3, 8)*sJ3_3 + covP(4, 8)*sJ3_4 + covP(6, 8)*sJ3_6 + covP(7, 8)*sJ3_7 + tmp29; 
  tmp30 = (((((stateJac[14] * covP[126] + stateJac[15] * covP[127]) + stateJac
              [16] * covP[128]) + stateJac[17] * covP[129]) + stateJac[18] *
            covP[131]) + stateJac[19] * covP[132]) + tmp29;

  // 'updateCovPNoGps:126' tmp31 = covP(1, 1)*sJ4_1 + covP(2, 1)*sJ4_2 + covP(3, 1)*sJ4_3 + covP(4, 1)*sJ4_4 + covP(6, 1)*sJ4_6 + covP(7, 1)*sJ4_7 + covP(8, 1)*sJ4_8; 
  tmp31 = (((((covP[0] * stateJac[21] + covP[1] * stateJac[22]) + covP[2] *
              stateJac[23]) + covP[3] * stateJac[24]) + covP[5] * stateJac[25])
           + covP[6] * stateJac[26]) + covP[7] * stateJac[27];

  // 'updateCovPNoGps:127' tmp32 = covP(1, 2)*sJ4_1 + covP(2, 2)*sJ4_2 + covP(3, 2)*sJ4_3 + covP(4, 2)*sJ4_4 + covP(6, 2)*sJ4_6 + covP(7, 2)*sJ4_7 + covP(8, 2)*sJ4_8; 
  tmp32 = (((((covP[18] * stateJac[21] + covP[19] * stateJac[22]) + covP[20] *
              stateJac[23]) + covP[21] * stateJac[24]) + covP[23] * stateJac[25])
           + covP[24] * stateJac[26]) + covP[25] * stateJac[27];

  // 'updateCovPNoGps:128' tmp33 = covP(1, 3)*sJ4_1 + covP(2, 3)*sJ4_2 + covP(3, 3)*sJ4_3 + covP(4, 3)*sJ4_4 + covP(6, 3)*sJ4_6 + covP(7, 3)*sJ4_7 + covP(8, 3)*sJ4_8; 
  tmp33 = (((((stateJac[21] * covP[36] + stateJac[22] * covP[37]) + stateJac[23]
              * covP[38]) + stateJac[24] * covP[39]) + stateJac[25] * covP[41])
           + stateJac[26] * covP[42]) + stateJac[27] * covP[43];

  // 'updateCovPNoGps:129' tmp34 = covP(1, 4)*sJ4_1 + covP(2, 4)*sJ4_2 + covP(3, 4)*sJ4_3 + covP(4, 4)*sJ4_4 + covP(6, 4)*sJ4_6 + covP(7, 4)*sJ4_7 + covP(8, 4)*sJ4_8; 
  tmp34 = (((((stateJac[21] * covP[54] + stateJac[22] * covP[55]) + stateJac[23]
              * covP[56]) + stateJac[24] * covP[57]) + stateJac[25] * covP[59])
           + stateJac[26] * covP[60]) + stateJac[27] * covP[61];

  // 'updateCovPNoGps:130' tmp35 = covP(6, 6)*sJ4_6;
  tmp35 = stateJac[25] * covP[95];

  // 'updateCovPNoGps:131' tmp36 = covP(1, 6)*sJ4_1 + covP(2, 6)*sJ4_2 + covP(3, 6)*sJ4_3 + covP(4, 6)*sJ4_4 + covP(7, 6)*sJ4_7 + covP(8, 6)*sJ4_8 + tmp35; 
  tmp36 = (((((stateJac[21] * covP[90] + stateJac[22] * covP[91]) + stateJac[23]
              * covP[92]) + stateJac[24] * covP[93]) + stateJac[26] * covP[96])
           + stateJac[27] * covP[97]) + tmp35;

  // 'updateCovPNoGps:132' tmp37 = covP(7, 7)*sJ4_7;
  tmp37 = stateJac[26] * covP[114];

  // 'updateCovPNoGps:133' tmp38 = covP(1, 7)*sJ4_1 + covP(2, 7)*sJ4_2 + covP(3, 7)*sJ4_3 + covP(4, 7)*sJ4_4 + covP(6, 7)*sJ4_6 + covP(8, 7)*sJ4_8 + tmp37; 
  tmp38 = (((((stateJac[21] * covP[108] + stateJac[22] * covP[109]) + stateJac
              [23] * covP[110]) + stateJac[24] * covP[111]) + stateJac[25] *
            covP[113]) + stateJac[27] * covP[115]) + tmp37;

  // 'updateCovPNoGps:134' tmp39 = covP(8, 8)*sJ4_8;
  tmp39 = stateJac[27] * covP[133];

  // 'updateCovPNoGps:135' tmp40 = covP(1, 8)*sJ4_1 + covP(2, 8)*sJ4_2 + covP(3, 8)*sJ4_3 + covP(4, 8)*sJ4_4 + covP(6, 8)*sJ4_6 + covP(7, 8)*sJ4_7 + tmp39; 
  tmp40 = (((((stateJac[21] * covP[126] + stateJac[22] * covP[127]) + stateJac
              [23] * covP[128]) + stateJac[24] * covP[129]) + stateJac[25] *
            covP[131]) + stateJac[26] * covP[132]) + tmp39;

  // 'updateCovPNoGps:136' tmp41 = covP(5, 1)*1;
  tmp41 = covP[4];

  // 'updateCovPNoGps:137' tmp42 = covP(5, 2)*1;
  tmp42 = covP[22];

  // 'updateCovPNoGps:138' tmp43 = covP(5, 3)*1;
  tmp43 = covP[40];

  // 'updateCovPNoGps:139' tmp44 = covP(5, 4)*1;
  tmp44 = covP[58];

  // 'updateCovPNoGps:140' tmp45 = covP(5, 6)*1;
  tmp45 = covP[94];

  // 'updateCovPNoGps:141' tmp46 = covP(5, 7)*1;
  tmp46 = covP[112];

  // 'updateCovPNoGps:142' tmp47 = covP(5, 8)*1;
  tmp47 = covP[130];

  // 'updateCovPNoGps:143' tmp48 = 1*1;
  // 'updateCovPNoGps:144' tmp49 = 1*1;
  // 'updateCovPNoGps:145' tmp50 = 1*1;
  // 'updateCovPNoGps:146' tmp51 = 1*1;
  // 'updateCovPNoGps:147' tmp52 = 1*1;
  // 'updateCovPNoGps:148' tmp53 = 1*1;
  // 'updateCovPNoGps:149' tmp54 = 1*1;
  // 'updateCovPNoGps:150' tmp55 = 1*1;
  // 'updateCovPNoGps:151' tmp56 = 1*1;
  // 'updateCovPNoGps:152' tmp57 = 1*1;
  // 'updateCovPNoGps:153' tmp58 = covP(6, 1)*1;
  tmp58 = covP[5];

  // 'updateCovPNoGps:154' tmp59 = covP(6, 2)*1;
  tmp59 = covP[23];

  // 'updateCovPNoGps:155' tmp60 = covP(6, 3)*1;
  tmp60 = covP[41];

  // 'updateCovPNoGps:156' tmp61 = covP(6, 4)*1;
  tmp61 = covP[59];

  // 'updateCovPNoGps:157' tmp62 = covP(6, 7)*1;
  tmp62 = covP[113];

  // 'updateCovPNoGps:158' tmp63 = covP(6, 8)*1;
  tmp63 = covP[131];

  // 'updateCovPNoGps:159' tmp64 = 1*1;
  // 'updateCovPNoGps:160' tmp65 = 1*1;
  // 'updateCovPNoGps:161' tmp66 = 1*1;
  // 'updateCovPNoGps:162' tmp67 = 1*1;
  // 'updateCovPNoGps:163' tmp68 = 1*1;
  // 'updateCovPNoGps:164' tmp69 = 1*1;
  // 'updateCovPNoGps:165' tmp70 = 1*1;
  // 'updateCovPNoGps:166' tmp71 = 1*1;
  // 'updateCovPNoGps:167' tmp72 = 1*1;
  // 'updateCovPNoGps:168' tmp73 = 1*1;
  // 'updateCovPNoGps:169' tmp74 = covP(7, 1)*1;
  tmp74 = covP[6];

  // 'updateCovPNoGps:170' tmp75 = covP(7, 2)*1;
  tmp75 = covP[24];

  // 'updateCovPNoGps:171' tmp76 = covP(7, 3)*1;
  tmp76 = covP[42];

  // 'updateCovPNoGps:172' tmp77 = covP(7, 4)*1;
  tmp77 = covP[60];

  // 'updateCovPNoGps:173' tmp78 = covP(7, 6)*1;
  tmp78 = covP[96];

  // 'updateCovPNoGps:174' tmp79 = covP(7, 8)*1;
  tmp79 = covP[132];

  // 'updateCovPNoGps:175' tmp80 = 1*1;
  // 'updateCovPNoGps:176' tmp81 = 1*1;
  // 'updateCovPNoGps:177' tmp82 = 1*1;
  // 'updateCovPNoGps:178' tmp83 = 1*1;
  // 'updateCovPNoGps:179' tmp84 = 1*1;
  // 'updateCovPNoGps:180' tmp85 = 1*1;
  // 'updateCovPNoGps:181' tmp86 = 1*1;
  // 'updateCovPNoGps:182' tmp87 = 1*1;
  // 'updateCovPNoGps:183' tmp88 = 1*1;
  // 'updateCovPNoGps:184' tmp89 = 1*1;
  // 'updateCovPNoGps:185' tmp90 = covP(8, 1)*1;
  tmp90 = covP[7];

  // 'updateCovPNoGps:186' tmp91 = covP(8, 2)*1;
  tmp91 = covP[25];

  // 'updateCovPNoGps:187' tmp92 = covP(8, 3)*1;
  tmp92 = covP[43];

  // 'updateCovPNoGps:188' tmp93 = covP(8, 4)*1;
  tmp93 = covP[61];

  // 'updateCovPNoGps:189' tmp94 = covP(8, 6)*1;
  tmp94 = covP[97];

  // 'updateCovPNoGps:190' tmp95 = covP(8, 7)*1;
  tmp95 = covP[115];

  // 'updateCovPNoGps:191' tmp96 = 1*1;
  // 'updateCovPNoGps:192' tmp97 = 1*1;
  // 'updateCovPNoGps:193' tmp98 = 1*1;
  // 'updateCovPNoGps:194' tmp99 = 1*1;
  // 'updateCovPNoGps:195' tmp100 = 1*1;
  // 'updateCovPNoGps:196' tmp101 = 1*1;
  // 'updateCovPNoGps:197' tmp102 = 1*1;
  // 'updateCovPNoGps:198' tmp103 = 1*1;
  // 'updateCovPNoGps:199' tmp104 = 1*1;
  // 'updateCovPNoGps:200' tmp105 = 1*1;
  // 'updateCovPNoGps:201' tmp106 = covP(9, 1)*1;
  tmp106 = covP[8];

  // 'updateCovPNoGps:202' tmp107 = covP(9, 2)*1;
  tmp107 = covP[26];

  // 'updateCovPNoGps:203' tmp108 = covP(9, 3)*1;
  tmp108 = covP[44];

  // 'updateCovPNoGps:204' tmp109 = covP(9, 4)*1;
  tmp109 = covP[62];

  // 'updateCovPNoGps:205' tmp110 = covP(9, 6)*1;
  tmp110 = covP[98];

  // 'updateCovPNoGps:206' tmp111 = covP(9, 7)*1;
  tmp111 = covP[116];

  // 'updateCovPNoGps:207' tmp112 = covP(9, 8)*1;
  tmp112 = covP[134];

  // 'updateCovPNoGps:208' tmp113 = 1*1;
  // 'updateCovPNoGps:209' tmp114 = 1*1;
  // 'updateCovPNoGps:210' tmp115 = 1*1;
  // 'updateCovPNoGps:211' tmp116 = 1*1;
  // 'updateCovPNoGps:212' tmp117 = 1*1;
  // 'updateCovPNoGps:213' tmp118 = 1*1;
  // 'updateCovPNoGps:214' tmp119 = 1*1;
  // 'updateCovPNoGps:215' tmp120 = 1*1;
  // 'updateCovPNoGps:216' tmp121 = 1*1;
  // 'updateCovPNoGps:217' tmp122 = covP(10, 1)*1;
  tmp122 = covP[9];

  // 'updateCovPNoGps:218' tmp123 = covP(10, 2)*1;
  tmp123 = covP[27];

  // 'updateCovPNoGps:219' tmp124 = covP(10, 3)*1;
  tmp124 = covP[45];

  // 'updateCovPNoGps:220' tmp125 = covP(10, 4)*1;
  tmp125 = covP[63];

  // 'updateCovPNoGps:221' tmp126 = covP(10, 6)*1;
  tmp126 = covP[99];

  // 'updateCovPNoGps:222' tmp127 = covP(10, 7)*1;
  tmp127 = covP[117];

  // 'updateCovPNoGps:223' tmp128 = covP(10, 8)*1;
  tmp128 = covP[135];

  // 'updateCovPNoGps:224' tmp129 = 1*1;
  // 'updateCovPNoGps:225' tmp130 = 1*1;
  // 'updateCovPNoGps:226' tmp131 = 1*1;
  // 'updateCovPNoGps:227' tmp132 = 1*1;
  // 'updateCovPNoGps:228' tmp133 = 1*1;
  // 'updateCovPNoGps:229' tmp134 = 1*1;
  // 'updateCovPNoGps:230' tmp135 = 1*1;
  // 'updateCovPNoGps:231' tmp136 = 1*1;
  // 'updateCovPNoGps:232' tmp137 = covP(11, 1)*1;
  tmp137 = covP[10];

  // 'updateCovPNoGps:233' tmp138 = covP(11, 2)*1;
  tmp138 = covP[28];

  // 'updateCovPNoGps:234' tmp139 = covP(11, 3)*1;
  tmp139 = covP[46];

  // 'updateCovPNoGps:235' tmp140 = covP(11, 4)*1;
  tmp140 = covP[64];

  // 'updateCovPNoGps:236' tmp141 = covP(11, 6)*1;
  tmp141 = covP[100];

  // 'updateCovPNoGps:237' tmp142 = covP(11, 7)*1;
  tmp142 = covP[118];

  // 'updateCovPNoGps:238' tmp143 = covP(11, 8)*1;
  tmp143 = covP[136];

  // 'updateCovPNoGps:239' tmp144 = 1*1;
  // 'updateCovPNoGps:240' tmp145 = 1*1;
  // 'updateCovPNoGps:241' tmp146 = 1*1;
  // 'updateCovPNoGps:242' tmp147 = 1*1;
  // 'updateCovPNoGps:243' tmp148 = 1*1;
  // 'updateCovPNoGps:244' tmp149 = 1*1;
  // 'updateCovPNoGps:245' tmp150 = 1*1;
  // 'updateCovPNoGps:246' tmp151 = covP(12, 1)*1;
  tmp151 = covP[11];

  // 'updateCovPNoGps:247' tmp152 = covP(12, 2)*1;
  tmp152 = covP[29];

  // 'updateCovPNoGps:248' tmp153 = covP(12, 3)*1;
  tmp153 = covP[47];

  // 'updateCovPNoGps:249' tmp154 = covP(12, 4)*1;
  tmp154 = covP[65];

  // 'updateCovPNoGps:250' tmp155 = covP(12, 6)*1;
  tmp155 = covP[101];

  // 'updateCovPNoGps:251' tmp156 = covP(12, 7)*1;
  tmp156 = covP[119];

  // 'updateCovPNoGps:252' tmp157 = covP(12, 8)*1;
  tmp157 = covP[137];

  // 'updateCovPNoGps:253' tmp158 = 1*1;
  // 'updateCovPNoGps:254' tmp159 = 1*1;
  // 'updateCovPNoGps:255' tmp160 = 1*1;
  // 'updateCovPNoGps:256' tmp161 = 1*1;
  // 'updateCovPNoGps:257' tmp162 = 1*1;
  // 'updateCovPNoGps:258' tmp163 = 1*1;
  // 'updateCovPNoGps:259' tmp164 = covP(13, 1)*1;
  tmp164 = covP[12];

  // 'updateCovPNoGps:260' tmp165 = covP(13, 2)*1;
  tmp165 = covP[30];

  // 'updateCovPNoGps:261' tmp166 = covP(13, 3)*1;
  tmp166 = covP[48];

  // 'updateCovPNoGps:262' tmp167 = covP(13, 4)*1;
  tmp167 = covP[66];

  // 'updateCovPNoGps:263' tmp168 = covP(13, 6)*1;
  tmp168 = covP[102];

  // 'updateCovPNoGps:264' tmp169 = covP(13, 7)*1;
  tmp169 = covP[120];

  // 'updateCovPNoGps:265' tmp170 = covP(13, 8)*1;
  tmp170 = covP[138];

  // 'updateCovPNoGps:266' tmp171 = 1*1;
  // 'updateCovPNoGps:267' tmp172 = 1*1;
  // 'updateCovPNoGps:268' tmp173 = 1*1;
  // 'updateCovPNoGps:269' tmp174 = 1*1;
  // 'updateCovPNoGps:270' tmp175 = 1*1;
  // 'updateCovPNoGps:271' tmp176 = covP(14, 1)*1;
  tmp176 = covP[13];

  // 'updateCovPNoGps:272' tmp177 = covP(14, 2)*1;
  tmp177 = covP[31];

  // 'updateCovPNoGps:273' tmp178 = covP(14, 3)*1;
  tmp178 = covP[49];

  // 'updateCovPNoGps:274' tmp179 = covP(14, 4)*1;
  tmp179 = covP[67];

  // 'updateCovPNoGps:275' tmp180 = covP(14, 6)*1;
  tmp180 = covP[103];

  // 'updateCovPNoGps:276' tmp181 = covP(14, 7)*1;
  tmp181 = covP[121];

  // 'updateCovPNoGps:277' tmp182 = covP(14, 8)*1;
  tmp182 = covP[139];

  // 'updateCovPNoGps:278' tmp183 = 1*1;
  // 'updateCovPNoGps:279' tmp184 = 1*1;
  // 'updateCovPNoGps:280' tmp185 = 1*1;
  // 'updateCovPNoGps:281' tmp186 = 1*1;
  // 'updateCovPNoGps:282' tmp187 = covP(15, 1)*1;
  tmp187 = covP[14];

  // 'updateCovPNoGps:283' tmp188 = covP(15, 2)*1;
  tmp188 = covP[32];

  // 'updateCovPNoGps:284' tmp189 = covP(15, 3)*1;
  tmp189 = covP[50];

  // 'updateCovPNoGps:285' tmp190 = covP(15, 4)*1;
  tmp190 = covP[68];

  // 'updateCovPNoGps:286' tmp191 = covP(15, 6)*1;
  tmp191 = covP[104];

  // 'updateCovPNoGps:287' tmp192 = covP(15, 7)*1;
  tmp192 = covP[122];

  // 'updateCovPNoGps:288' tmp193 = covP(15, 8)*1;
  tmp193 = covP[140];

  // 'updateCovPNoGps:289' tmp194 = 1*1;
  // 'updateCovPNoGps:290' tmp195 = 1*1;
  // 'updateCovPNoGps:291' tmp196 = 1*1;
  // 'updateCovPNoGps:292' tmp197 = covP(16, 1)*1;
  tmp197 = covP[15];

  // 'updateCovPNoGps:293' tmp198 = covP(16, 2)*1;
  tmp198 = covP[33];

  // 'updateCovPNoGps:294' tmp199 = covP(16, 3)*1;
  tmp199 = covP[51];

  // 'updateCovPNoGps:295' tmp200 = covP(16, 4)*1;
  tmp200 = covP[69];

  // 'updateCovPNoGps:296' tmp201 = covP(16, 6)*1;
  tmp201 = covP[105];

  // 'updateCovPNoGps:297' tmp202 = covP(16, 7)*1;
  tmp202 = covP[123];

  // 'updateCovPNoGps:298' tmp203 = covP(16, 8)*1;
  tmp203 = covP[141];

  // 'updateCovPNoGps:299' tmp204 = 1*1;
  // 'updateCovPNoGps:300' tmp205 = 1*1;
  // 'updateCovPNoGps:301' tmp206 = covP(17, 1)*1;
  tmp206 = covP[16];

  // 'updateCovPNoGps:302' tmp207 = covP(17, 2)*1;
  tmp207 = covP[34];

  // 'updateCovPNoGps:303' tmp208 = covP(17, 3)*1;
  tmp208 = covP[52];

  // 'updateCovPNoGps:304' tmp209 = covP(17, 4)*1;
  tmp209 = covP[70];

  // 'updateCovPNoGps:305' tmp210 = covP(17, 6)*1;
  tmp210 = covP[106];

  // 'updateCovPNoGps:306' tmp211 = covP(17, 7)*1;
  tmp211 = covP[124];

  // 'updateCovPNoGps:307' tmp212 = covP(17, 8)*1;
  tmp212 = covP[142];

  // 'updateCovPNoGps:308' tmp213 = 1*1;
  // 'updateCovPNoGps:309' tmp214 = covP(18, 1)*1;
  tmp214 = covP[17];

  // 'updateCovPNoGps:310' tmp215 = covP(18, 2)*1;
  tmp215 = covP[35];

  // 'updateCovPNoGps:311' tmp216 = covP(18, 3)*1;
  tmp216 = covP[53];

  // 'updateCovPNoGps:312' tmp217 = covP(18, 4)*1;
  tmp217 = covP[71];

  // 'updateCovPNoGps:313' tmp218 = covP(18, 6)*1;
  tmp218 = covP[107];

  // 'updateCovPNoGps:314' tmp219 = covP(18, 7)*1;
  tmp219 = covP[125];

  // 'updateCovPNoGps:315' tmp220 = covP(18, 8)*1;
  tmp220 = covP[143];

  // 'updateCovPNoGps:316' covP(1, 1) = processNoiseQ(1, 1) + sJ1_1*tmp1 + sJ1_2*tmp2 + sJ1_3*tmp3 + sJ1_4*tmp4 + sJ1_6*tmp6 + sJ1_7*tmp8 + sJ1_8*tmp10; 
  covP[0] = ((((((stateJac[0] * tmp1 + processNoiseQ[0]) + stateJac[1] * tmp2) +
                stateJac[2] * tmp3) + stateJac[3] * tmp4) + stateJac[4] * tmp6)
             + stateJac[5] * tmp8) + stateJac[6] * tmp10;

  // 'updateCovPNoGps:317' covP(1, 2) = sJ2_1*tmp1 + sJ2_2*tmp2 + sJ2_3*tmp3 + sJ2_4*tmp4 + sJ2_6*tmp6 + sJ2_7*tmp8 + sJ2_8*tmp10; 
  covP[18] = (((((stateJac[7] * tmp1 + stateJac[8] * tmp2) + stateJac[9] * tmp3)
                + stateJac[10] * tmp4) + stateJac[11] * tmp6) + stateJac[12] *
              tmp8) + stateJac[13] * tmp10;

  // 'updateCovPNoGps:318' covP(1, 3) = sJ3_1*tmp1 + sJ3_2*tmp2 + sJ3_3*tmp3 + sJ3_4*tmp4 + sJ3_6*tmp6 + sJ3_7*tmp8 + sJ3_8*tmp10; 
  covP[36] = (((((stateJac[14] * tmp1 + stateJac[15] * tmp2) + stateJac[16] *
                 tmp3) + stateJac[17] * tmp4) + stateJac[18] * tmp6) + stateJac
              [19] * tmp8) + stateJac[20] * tmp10;

  // 'updateCovPNoGps:319' covP(1, 4) = sJ4_1*tmp1 + sJ4_2*tmp2 + sJ4_3*tmp3 + sJ4_4*tmp4 + sJ4_6*tmp6 + sJ4_7*tmp8 + sJ4_8*tmp10; 
  covP[54] = (((((stateJac[21] * tmp1 + stateJac[22] * tmp2) + stateJac[23] *
                 tmp3) + stateJac[24] * tmp4) + stateJac[25] * tmp6) + stateJac
              [26] * tmp8) + stateJac[27] * tmp10;

  // 'updateCovPNoGps:320' covP(1, 5) = 1*(covP(1, 5)*sJ1_1 + covP(2, 5)*sJ1_2 + covP(3, 5)*sJ1_3 + covP(4, 5)*sJ1_4 + covP(6, 5)*sJ1_6 + covP(7, 5)*sJ1_7 + covP(8, 5)*sJ1_8); 
  covP[72] = (((((stateJac[0] * covP[72] + stateJac[1] * covP[73]) + stateJac[2]
                 * covP[74]) + stateJac[3] * covP[75]) + stateJac[4] * covP[77])
              + stateJac[5] * covP[78]) + stateJac[6] * covP[79];

  // 'updateCovPNoGps:321' covP(1, 6) = 1*tmp6;
  covP[90] = tmp6;

  // 'updateCovPNoGps:322' covP(1, 7) = 1*tmp8;
  covP[108] = tmp8;

  // 'updateCovPNoGps:323' covP(1, 8) = 1*tmp10;
  covP[126] = tmp10;

  // 'updateCovPNoGps:324' covP(1, 9) = 1*(covP(1, 9)*sJ1_1 + covP(2, 9)*sJ1_2 + covP(3, 9)*sJ1_3 + covP(4, 9)*sJ1_4 + covP(6, 9)*sJ1_6 + covP(7, 9)*sJ1_7 + covP(8, 9)*sJ1_8); 
  covP[144] = (((((stateJac[0] * covP[144] + stateJac[1] * covP[145]) +
                  stateJac[2] * covP[146]) + stateJac[3] * covP[147]) +
                stateJac[4] * covP[149]) + stateJac[5] * covP[150]) + stateJac[6]
    * covP[151];

  // 'updateCovPNoGps:325' covP(1, 10) = 1*(covP(1, 10)*sJ1_1 + covP(2, 10)*sJ1_2 + covP(3, 10)*sJ1_3 + covP(4, 10)*sJ1_4 + covP(6, 10)*sJ1_6 + covP(7, 10)*sJ1_7 + covP(8, 10)*sJ1_8); 
  covP[162] = (((((stateJac[0] * covP[162] + stateJac[1] * covP[163]) +
                  stateJac[2] * covP[164]) + stateJac[3] * covP[165]) +
                stateJac[4] * covP[167]) + stateJac[5] * covP[168]) + stateJac[6]
    * covP[169];

  // 'updateCovPNoGps:326' covP(1, 11) = 1*(covP(1, 11)*sJ1_1 + covP(2, 11)*sJ1_2 + covP(3, 11)*sJ1_3 + covP(4, 11)*sJ1_4 + covP(6, 11)*sJ1_6 + covP(7, 11)*sJ1_7 + covP(8, 11)*sJ1_8); 
  covP[180] = (((((stateJac[0] * covP[180] + stateJac[1] * covP[181]) +
                  stateJac[2] * covP[182]) + stateJac[3] * covP[183]) +
                stateJac[4] * covP[185]) + stateJac[5] * covP[186]) + stateJac[6]
    * covP[187];

  // 'updateCovPNoGps:327' covP(1, 12) = 1*(covP(1, 12)*sJ1_1 + covP(2, 12)*sJ1_2 + covP(3, 12)*sJ1_3 + covP(4, 12)*sJ1_4 + covP(6, 12)*sJ1_6 + covP(7, 12)*sJ1_7 + covP(8, 12)*sJ1_8); 
  covP[198] = (((((stateJac[0] * covP[198] + stateJac[1] * covP[199]) +
                  stateJac[2] * covP[200]) + stateJac[3] * covP[201]) +
                stateJac[4] * covP[203]) + stateJac[5] * covP[204]) + stateJac[6]
    * covP[205];

  // 'updateCovPNoGps:328' covP(1, 13) = 1*(covP(1, 13)*sJ1_1 + covP(2, 13)*sJ1_2 + covP(3, 13)*sJ1_3 + covP(4, 13)*sJ1_4 + covP(6, 13)*sJ1_6 + covP(7, 13)*sJ1_7 + covP(8, 13)*sJ1_8); 
  covP[216] = (((((stateJac[0] * covP[216] + stateJac[1] * covP[217]) +
                  stateJac[2] * covP[218]) + stateJac[3] * covP[219]) +
                stateJac[4] * covP[221]) + stateJac[5] * covP[222]) + stateJac[6]
    * covP[223];

  // 'updateCovPNoGps:329' covP(1, 14) = 1*(covP(1, 14)*sJ1_1 + covP(2, 14)*sJ1_2 + covP(3, 14)*sJ1_3 + covP(4, 14)*sJ1_4 + covP(6, 14)*sJ1_6 + covP(7, 14)*sJ1_7 + covP(8, 14)*sJ1_8); 
  covP[234] = (((((stateJac[0] * covP[234] + stateJac[1] * covP[235]) +
                  stateJac[2] * covP[236]) + stateJac[3] * covP[237]) +
                stateJac[4] * covP[239]) + stateJac[5] * covP[240]) + stateJac[6]
    * covP[241];

  // 'updateCovPNoGps:330' covP(1, 15) = 1*(covP(1, 15)*sJ1_1 + covP(2, 15)*sJ1_2 + covP(3, 15)*sJ1_3 + covP(4, 15)*sJ1_4 + covP(6, 15)*sJ1_6 + covP(7, 15)*sJ1_7 + covP(8, 15)*sJ1_8); 
  covP[252] = (((((stateJac[0] * covP[252] + stateJac[1] * covP[253]) +
                  stateJac[2] * covP[254]) + stateJac[3] * covP[255]) +
                stateJac[4] * covP[257]) + stateJac[5] * covP[258]) + stateJac[6]
    * covP[259];

  // 'updateCovPNoGps:331' covP(1, 16) = 1*(covP(1, 16)*sJ1_1 + covP(2, 16)*sJ1_2 + covP(3, 16)*sJ1_3 + covP(4, 16)*sJ1_4 + covP(6, 16)*sJ1_6 + covP(7, 16)*sJ1_7 + covP(8, 16)*sJ1_8); 
  covP[270] = (((((stateJac[0] * covP[270] + stateJac[1] * covP[271]) +
                  stateJac[2] * covP[272]) + stateJac[3] * covP[273]) +
                stateJac[4] * covP[275]) + stateJac[5] * covP[276]) + stateJac[6]
    * covP[277];

  // 'updateCovPNoGps:332' covP(1, 17) = 1*(covP(1, 17)*sJ1_1 + covP(2, 17)*sJ1_2 + covP(3, 17)*sJ1_3 + covP(4, 17)*sJ1_4 + covP(6, 17)*sJ1_6 + covP(7, 17)*sJ1_7 + covP(8, 17)*sJ1_8); 
  covP[288] = (((((stateJac[0] * covP[288] + stateJac[1] * covP[289]) +
                  stateJac[2] * covP[290]) + stateJac[3] * covP[291]) +
                stateJac[4] * covP[293]) + stateJac[5] * covP[294]) + stateJac[6]
    * covP[295];

  // 'updateCovPNoGps:333' covP(1, 18) = 1*(covP(1, 18)*sJ1_1 + covP(2, 18)*sJ1_2 + covP(3, 18)*sJ1_3 + covP(4, 18)*sJ1_4 + covP(6, 18)*sJ1_6 + covP(7, 18)*sJ1_7 + covP(8, 18)*sJ1_8); 
  covP[306] = (((((stateJac[0] * covP[306] + stateJac[1] * covP[307]) +
                  stateJac[2] * covP[308]) + stateJac[3] * covP[309]) +
                stateJac[4] * covP[311]) + stateJac[5] * covP[312]) + stateJac[6]
    * covP[313];

  // 'updateCovPNoGps:334' covP(2, 1) = sJ1_1*tmp11 + sJ1_2*tmp12 + sJ1_3*tmp13 + sJ1_4*tmp14 + sJ1_6*tmp16 + sJ1_7*tmp18 + sJ1_8*tmp20; 
  covP[1] = (((((stateJac[0] * tmp11 + stateJac[1] * tmp12) + stateJac[2] *
                tmp13) + stateJac[3] * tmp14) + stateJac[4] * tmp16) + stateJac
             [5] * tmp18) + stateJac[6] * tmp20;

  // 'updateCovPNoGps:335' covP(2, 2) = processNoiseQ(2, 2) + sJ2_1*tmp11 + sJ2_2*tmp12 + sJ2_3*tmp13 + sJ2_4*tmp14 + sJ2_6*tmp16 + sJ2_7*tmp18 + sJ2_8*tmp20; 
  covP[19] = ((((((stateJac[7] * tmp11 + processNoiseQ[19]) + stateJac[8] *
                  tmp12) + stateJac[9] * tmp13) + stateJac[10] * tmp14) +
               stateJac[11] * tmp16) + stateJac[12] * tmp18) + stateJac[13] *
    tmp20;

  // 'updateCovPNoGps:336' covP(2, 3) = sJ3_1*tmp11 + sJ3_2*tmp12 + sJ3_3*tmp13 + sJ3_4*tmp14 + sJ3_6*tmp16 + sJ3_7*tmp18 + sJ3_8*tmp20; 
  covP[37] = (((((stateJac[14] * tmp11 + stateJac[15] * tmp12) + stateJac[16] *
                 tmp13) + stateJac[17] * tmp14) + stateJac[18] * tmp16) +
              stateJac[19] * tmp18) + stateJac[20] * tmp20;

  // 'updateCovPNoGps:337' covP(2, 4) = sJ4_1*tmp11 + sJ4_2*tmp12 + sJ4_3*tmp13 + sJ4_4*tmp14 + sJ4_6*tmp16 + sJ4_7*tmp18 + sJ4_8*tmp20; 
  covP[55] = (((((stateJac[21] * tmp11 + stateJac[22] * tmp12) + stateJac[23] *
                 tmp13) + stateJac[24] * tmp14) + stateJac[25] * tmp16) +
              stateJac[26] * tmp18) + stateJac[27] * tmp20;

  // 'updateCovPNoGps:338' covP(2, 5) = 1*(covP(1, 5)*sJ2_1 + covP(2, 5)*sJ2_2 + covP(3, 5)*sJ2_3 + covP(4, 5)*sJ2_4 + covP(6, 5)*sJ2_6 + covP(7, 5)*sJ2_7 + covP(8, 5)*sJ2_8); 
  covP[73] = (((((stateJac[7] * covP[72] + stateJac[8] * covP[73]) + stateJac[9]
                 * covP[74]) + stateJac[10] * covP[75]) + stateJac[11] * covP[77])
              + stateJac[12] * covP[78]) + stateJac[13] * covP[79];

  // 'updateCovPNoGps:339' covP(2, 6) = 1*tmp16;
  covP[91] = tmp16;

  // 'updateCovPNoGps:340' covP(2, 7) = 1*tmp18;
  covP[109] = tmp18;

  // 'updateCovPNoGps:341' covP(2, 8) = 1*tmp20;
  covP[127] = tmp20;

  // 'updateCovPNoGps:342' covP(2, 9) = 1*(covP(1, 9)*sJ2_1 + covP(2, 9)*sJ2_2 + covP(3, 9)*sJ2_3 + covP(4, 9)*sJ2_4 + covP(6, 9)*sJ2_6 + covP(7, 9)*sJ2_7 + covP(8, 9)*sJ2_8); 
  covP[145] = (((((stateJac[7] * covP[144] + stateJac[8] * covP[145]) +
                  stateJac[9] * covP[146]) + stateJac[10] * covP[147]) +
                stateJac[11] * covP[149]) + stateJac[12] * covP[150]) +
    stateJac[13] * covP[151];

  // 'updateCovPNoGps:343' covP(2, 10) = 1*(covP(1, 10)*sJ2_1 + covP(2, 10)*sJ2_2 + covP(3, 10)*sJ2_3 + covP(4, 10)*sJ2_4 + covP(6, 10)*sJ2_6 + covP(7, 10)*sJ2_7 + covP(8, 10)*sJ2_8); 
  covP[163] = (((((stateJac[7] * covP[162] + stateJac[8] * covP[163]) +
                  stateJac[9] * covP[164]) + stateJac[10] * covP[165]) +
                stateJac[11] * covP[167]) + stateJac[12] * covP[168]) +
    stateJac[13] * covP[169];

  // 'updateCovPNoGps:344' covP(2, 11) = 1*(covP(1, 11)*sJ2_1 + covP(2, 11)*sJ2_2 + covP(3, 11)*sJ2_3 + covP(4, 11)*sJ2_4 + covP(6, 11)*sJ2_6 + covP(7, 11)*sJ2_7 + covP(8, 11)*sJ2_8); 
  covP[181] = (((((stateJac[7] * covP[180] + stateJac[8] * covP[181]) +
                  stateJac[9] * covP[182]) + stateJac[10] * covP[183]) +
                stateJac[11] * covP[185]) + stateJac[12] * covP[186]) +
    stateJac[13] * covP[187];

  // 'updateCovPNoGps:345' covP(2, 12) = 1*(covP(1, 12)*sJ2_1 + covP(2, 12)*sJ2_2 + covP(3, 12)*sJ2_3 + covP(4, 12)*sJ2_4 + covP(6, 12)*sJ2_6 + covP(7, 12)*sJ2_7 + covP(8, 12)*sJ2_8); 
  covP[199] = (((((stateJac[7] * covP[198] + stateJac[8] * covP[199]) +
                  stateJac[9] * covP[200]) + stateJac[10] * covP[201]) +
                stateJac[11] * covP[203]) + stateJac[12] * covP[204]) +
    stateJac[13] * covP[205];

  // 'updateCovPNoGps:346' covP(2, 13) = 1*(covP(1, 13)*sJ2_1 + covP(2, 13)*sJ2_2 + covP(3, 13)*sJ2_3 + covP(4, 13)*sJ2_4 + covP(6, 13)*sJ2_6 + covP(7, 13)*sJ2_7 + covP(8, 13)*sJ2_8); 
  covP[217] = (((((stateJac[7] * covP[216] + stateJac[8] * covP[217]) +
                  stateJac[9] * covP[218]) + stateJac[10] * covP[219]) +
                stateJac[11] * covP[221]) + stateJac[12] * covP[222]) +
    stateJac[13] * covP[223];

  // 'updateCovPNoGps:347' covP(2, 14) = 1*(covP(1, 14)*sJ2_1 + covP(2, 14)*sJ2_2 + covP(3, 14)*sJ2_3 + covP(4, 14)*sJ2_4 + covP(6, 14)*sJ2_6 + covP(7, 14)*sJ2_7 + covP(8, 14)*sJ2_8); 
  covP[235] = (((((stateJac[7] * covP[234] + stateJac[8] * covP[235]) +
                  stateJac[9] * covP[236]) + stateJac[10] * covP[237]) +
                stateJac[11] * covP[239]) + stateJac[12] * covP[240]) +
    stateJac[13] * covP[241];

  // 'updateCovPNoGps:348' covP(2, 15) = 1*(covP(1, 15)*sJ2_1 + covP(2, 15)*sJ2_2 + covP(3, 15)*sJ2_3 + covP(4, 15)*sJ2_4 + covP(6, 15)*sJ2_6 + covP(7, 15)*sJ2_7 + covP(8, 15)*sJ2_8); 
  covP[253] = (((((stateJac[7] * covP[252] + stateJac[8] * covP[253]) +
                  stateJac[9] * covP[254]) + stateJac[10] * covP[255]) +
                stateJac[11] * covP[257]) + stateJac[12] * covP[258]) +
    stateJac[13] * covP[259];

  // 'updateCovPNoGps:349' covP(2, 16) = 1*(covP(1, 16)*sJ2_1 + covP(2, 16)*sJ2_2 + covP(3, 16)*sJ2_3 + covP(4, 16)*sJ2_4 + covP(6, 16)*sJ2_6 + covP(7, 16)*sJ2_7 + covP(8, 16)*sJ2_8); 
  covP[271] = (((((stateJac[7] * covP[270] + stateJac[8] * covP[271]) +
                  stateJac[9] * covP[272]) + stateJac[10] * covP[273]) +
                stateJac[11] * covP[275]) + stateJac[12] * covP[276]) +
    stateJac[13] * covP[277];

  // 'updateCovPNoGps:350' covP(2, 17) = 1*(covP(1, 17)*sJ2_1 + covP(2, 17)*sJ2_2 + covP(3, 17)*sJ2_3 + covP(4, 17)*sJ2_4 + covP(6, 17)*sJ2_6 + covP(7, 17)*sJ2_7 + covP(8, 17)*sJ2_8); 
  covP[289] = (((((stateJac[7] * covP[288] + stateJac[8] * covP[289]) +
                  stateJac[9] * covP[290]) + stateJac[10] * covP[291]) +
                stateJac[11] * covP[293]) + stateJac[12] * covP[294]) +
    stateJac[13] * covP[295];

  // 'updateCovPNoGps:351' covP(2, 18) = 1*(covP(1, 18)*sJ2_1 + covP(2, 18)*sJ2_2 + covP(3, 18)*sJ2_3 + covP(4, 18)*sJ2_4 + covP(6, 18)*sJ2_6 + covP(7, 18)*sJ2_7 + covP(8, 18)*sJ2_8); 
  covP[307] = (((((stateJac[7] * covP[306] + stateJac[8] * covP[307]) +
                  stateJac[9] * covP[308]) + stateJac[10] * covP[309]) +
                stateJac[11] * covP[311]) + stateJac[12] * covP[312]) +
    stateJac[13] * covP[313];

  // 'updateCovPNoGps:352' covP(3, 1) = sJ1_1*tmp21 + sJ1_2*tmp22 + sJ1_3*tmp23 + sJ1_4*tmp24 + sJ1_6*tmp26 + sJ1_7*tmp28 + sJ1_8*tmp30; 
  covP[2] = (((((stateJac[0] * tmp21 + stateJac[1] * tmp22) + stateJac[2] *
                tmp23) + stateJac[3] * tmp24) + stateJac[4] * tmp26) + stateJac
             [5] * tmp28) + stateJac[6] * tmp30;

  // 'updateCovPNoGps:353' covP(3, 2) = sJ2_1*tmp21 + sJ2_2*tmp22 + sJ2_3*tmp23 + sJ2_4*tmp24 + sJ2_6*tmp26 + sJ2_7*tmp28 + sJ2_8*tmp30; 
  covP[20] = (((((stateJac[7] * tmp21 + stateJac[8] * tmp22) + stateJac[9] *
                 tmp23) + stateJac[10] * tmp24) + stateJac[11] * tmp26) +
              stateJac[12] * tmp28) + stateJac[13] * tmp30;

  // 'updateCovPNoGps:354' covP(3, 3) = processNoiseQ(3, 3) + sJ3_1*tmp21 + sJ3_2*tmp22 + sJ3_3*tmp23 + sJ3_4*tmp24 + sJ3_6*tmp26 + sJ3_7*tmp28 + sJ3_8*tmp30; 
  covP[38] = ((((((stateJac[14] * tmp21 + processNoiseQ[38]) + stateJac[15] *
                  tmp22) + stateJac[16] * tmp23) + stateJac[17] * tmp24) +
               stateJac[18] * tmp26) + stateJac[19] * tmp28) + stateJac[20] *
    tmp30;

  // 'updateCovPNoGps:355' covP(3, 4) = sJ4_1*tmp21 + sJ4_2*tmp22 + sJ4_3*tmp23 + sJ4_4*tmp24 + sJ4_6*tmp26 + sJ4_7*tmp28 + sJ4_8*tmp30; 
  covP[56] = (((((stateJac[21] * tmp21 + stateJac[22] * tmp22) + stateJac[23] *
                 tmp23) + stateJac[24] * tmp24) + stateJac[25] * tmp26) +
              stateJac[26] * tmp28) + stateJac[27] * tmp30;

  // 'updateCovPNoGps:356' covP(3, 5) = 1*(covP(1, 5)*sJ3_1 + covP(2, 5)*sJ3_2 + covP(3, 5)*sJ3_3 + covP(4, 5)*sJ3_4 + covP(6, 5)*sJ3_6 + covP(7, 5)*sJ3_7 + covP(8, 5)*sJ3_8); 
  covP[74] = (((((stateJac[14] * covP[72] + stateJac[15] * covP[73]) + stateJac
                 [16] * covP[74]) + stateJac[17] * covP[75]) + stateJac[18] *
               covP[77]) + stateJac[19] * covP[78]) + stateJac[20] * covP[79];

  // 'updateCovPNoGps:357' covP(3, 6) = 1*tmp26;
  covP[92] = tmp26;

  // 'updateCovPNoGps:358' covP(3, 7) = 1*tmp28;
  covP[110] = tmp28;

  // 'updateCovPNoGps:359' covP(3, 8) = 1*tmp30;
  covP[128] = tmp30;

  // 'updateCovPNoGps:360' covP(3, 9) = 1*(covP(1, 9)*sJ3_1 + covP(2, 9)*sJ3_2 + covP(3, 9)*sJ3_3 + covP(4, 9)*sJ3_4 + covP(6, 9)*sJ3_6 + covP(7, 9)*sJ3_7 + covP(8, 9)*sJ3_8); 
  covP[146] = (((((stateJac[14] * covP[144] + stateJac[15] * covP[145]) +
                  stateJac[16] * covP[146]) + stateJac[17] * covP[147]) +
                stateJac[18] * covP[149]) + stateJac[19] * covP[150]) +
    stateJac[20] * covP[151];

  // 'updateCovPNoGps:361' covP(3, 10) = 1*(covP(1, 10)*sJ3_1 + covP(2, 10)*sJ3_2 + covP(3, 10)*sJ3_3 + covP(4, 10)*sJ3_4 + covP(6, 10)*sJ3_6 + covP(7, 10)*sJ3_7 + covP(8, 10)*sJ3_8); 
  covP[164] = (((((stateJac[14] * covP[162] + stateJac[15] * covP[163]) +
                  stateJac[16] * covP[164]) + stateJac[17] * covP[165]) +
                stateJac[18] * covP[167]) + stateJac[19] * covP[168]) +
    stateJac[20] * covP[169];

  // 'updateCovPNoGps:362' covP(3, 11) = 1*(covP(1, 11)*sJ3_1 + covP(2, 11)*sJ3_2 + covP(3, 11)*sJ3_3 + covP(4, 11)*sJ3_4 + covP(6, 11)*sJ3_6 + covP(7, 11)*sJ3_7 + covP(8, 11)*sJ3_8); 
  covP[182] = (((((stateJac[14] * covP[180] + stateJac[15] * covP[181]) +
                  stateJac[16] * covP[182]) + stateJac[17] * covP[183]) +
                stateJac[18] * covP[185]) + stateJac[19] * covP[186]) +
    stateJac[20] * covP[187];

  // 'updateCovPNoGps:363' covP(3, 12) = 1*(covP(1, 12)*sJ3_1 + covP(2, 12)*sJ3_2 + covP(3, 12)*sJ3_3 + covP(4, 12)*sJ3_4 + covP(6, 12)*sJ3_6 + covP(7, 12)*sJ3_7 + covP(8, 12)*sJ3_8); 
  covP[200] = (((((stateJac[14] * covP[198] + stateJac[15] * covP[199]) +
                  stateJac[16] * covP[200]) + stateJac[17] * covP[201]) +
                stateJac[18] * covP[203]) + stateJac[19] * covP[204]) +
    stateJac[20] * covP[205];

  // 'updateCovPNoGps:364' covP(3, 13) = 1*(covP(1, 13)*sJ3_1 + covP(2, 13)*sJ3_2 + covP(3, 13)*sJ3_3 + covP(4, 13)*sJ3_4 + covP(6, 13)*sJ3_6 + covP(7, 13)*sJ3_7 + covP(8, 13)*sJ3_8); 
  covP[218] = (((((stateJac[14] * covP[216] + stateJac[15] * covP[217]) +
                  stateJac[16] * covP[218]) + stateJac[17] * covP[219]) +
                stateJac[18] * covP[221]) + stateJac[19] * covP[222]) +
    stateJac[20] * covP[223];

  // 'updateCovPNoGps:365' covP(3, 14) = 1*(covP(1, 14)*sJ3_1 + covP(2, 14)*sJ3_2 + covP(3, 14)*sJ3_3 + covP(4, 14)*sJ3_4 + covP(6, 14)*sJ3_6 + covP(7, 14)*sJ3_7 + covP(8, 14)*sJ3_8); 
  covP[236] = (((((stateJac[14] * covP[234] + stateJac[15] * covP[235]) +
                  stateJac[16] * covP[236]) + stateJac[17] * covP[237]) +
                stateJac[18] * covP[239]) + stateJac[19] * covP[240]) +
    stateJac[20] * covP[241];

  // 'updateCovPNoGps:366' covP(3, 15) = 1*(covP(1, 15)*sJ3_1 + covP(2, 15)*sJ3_2 + covP(3, 15)*sJ3_3 + covP(4, 15)*sJ3_4 + covP(6, 15)*sJ3_6 + covP(7, 15)*sJ3_7 + covP(8, 15)*sJ3_8); 
  covP[254] = (((((stateJac[14] * covP[252] + stateJac[15] * covP[253]) +
                  stateJac[16] * covP[254]) + stateJac[17] * covP[255]) +
                stateJac[18] * covP[257]) + stateJac[19] * covP[258]) +
    stateJac[20] * covP[259];

  // 'updateCovPNoGps:367' covP(3, 16) = 1*(covP(1, 16)*sJ3_1 + covP(2, 16)*sJ3_2 + covP(3, 16)*sJ3_3 + covP(4, 16)*sJ3_4 + covP(6, 16)*sJ3_6 + covP(7, 16)*sJ3_7 + covP(8, 16)*sJ3_8); 
  covP[272] = (((((stateJac[14] * covP[270] + stateJac[15] * covP[271]) +
                  stateJac[16] * covP[272]) + stateJac[17] * covP[273]) +
                stateJac[18] * covP[275]) + stateJac[19] * covP[276]) +
    stateJac[20] * covP[277];

  // 'updateCovPNoGps:368' covP(3, 17) = 1*(covP(1, 17)*sJ3_1 + covP(2, 17)*sJ3_2 + covP(3, 17)*sJ3_3 + covP(4, 17)*sJ3_4 + covP(6, 17)*sJ3_6 + covP(7, 17)*sJ3_7 + covP(8, 17)*sJ3_8); 
  covP[290] = (((((stateJac[14] * covP[288] + stateJac[15] * covP[289]) +
                  stateJac[16] * covP[290]) + stateJac[17] * covP[291]) +
                stateJac[18] * covP[293]) + stateJac[19] * covP[294]) +
    stateJac[20] * covP[295];

  // 'updateCovPNoGps:369' covP(3, 18) = 1*(covP(1, 18)*sJ3_1 + covP(2, 18)*sJ3_2 + covP(3, 18)*sJ3_3 + covP(4, 18)*sJ3_4 + covP(6, 18)*sJ3_6 + covP(7, 18)*sJ3_7 + covP(8, 18)*sJ3_8); 
  covP[308] = (((((stateJac[14] * covP[306] + stateJac[15] * covP[307]) +
                  stateJac[16] * covP[308]) + stateJac[17] * covP[309]) +
                stateJac[18] * covP[311]) + stateJac[19] * covP[312]) +
    stateJac[20] * covP[313];

  // 'updateCovPNoGps:370' covP(4, 1) = sJ1_1*tmp31 + sJ1_2*tmp32 + sJ1_3*tmp33 + sJ1_4*tmp34 + sJ1_6*tmp36 + sJ1_7*tmp38 + sJ1_8*tmp40; 
  covP[3] = (((((stateJac[0] * tmp31 + stateJac[1] * tmp32) + stateJac[2] *
                tmp33) + stateJac[3] * tmp34) + stateJac[4] * tmp36) + stateJac
             [5] * tmp38) + stateJac[6] * tmp40;

  // 'updateCovPNoGps:371' covP(4, 2) = sJ2_1*tmp31 + sJ2_2*tmp32 + sJ2_3*tmp33 + sJ2_4*tmp34 + sJ2_6*tmp36 + sJ2_7*tmp38 + sJ2_8*tmp40; 
  covP[21] = (((((stateJac[7] * tmp31 + stateJac[8] * tmp32) + stateJac[9] *
                 tmp33) + stateJac[10] * tmp34) + stateJac[11] * tmp36) +
              stateJac[12] * tmp38) + stateJac[13] * tmp40;

  // 'updateCovPNoGps:372' covP(4, 3) = sJ3_1*tmp31 + sJ3_2*tmp32 + sJ3_3*tmp33 + sJ3_4*tmp34 + sJ3_6*tmp36 + sJ3_7*tmp38 + sJ3_8*tmp40; 
  covP[39] = (((((stateJac[14] * tmp31 + stateJac[15] * tmp32) + stateJac[16] *
                 tmp33) + stateJac[17] * tmp34) + stateJac[18] * tmp36) +
              stateJac[19] * tmp38) + stateJac[20] * tmp40;

  // 'updateCovPNoGps:373' covP(4, 4) = processNoiseQ(4, 4) + sJ4_1*tmp31 + sJ4_2*tmp32 + sJ4_3*tmp33 + sJ4_4*tmp34 + sJ4_6*tmp36 + sJ4_7*tmp38 + sJ4_8*tmp40; 
  covP[57] = ((((((stateJac[21] * tmp31 + processNoiseQ[57]) + stateJac[22] *
                  tmp32) + stateJac[23] * tmp33) + stateJac[24] * tmp34) +
               stateJac[25] * tmp36) + stateJac[26] * tmp38) + stateJac[27] *
    tmp40;

  // 'updateCovPNoGps:374' covP(4, 5) = 1*(covP(1, 5)*sJ4_1 + covP(2, 5)*sJ4_2 + covP(3, 5)*sJ4_3 + covP(4, 5)*sJ4_4 + covP(6, 5)*sJ4_6 + covP(7, 5)*sJ4_7 + covP(8, 5)*sJ4_8); 
  covP[75] = (((((stateJac[21] * covP[72] + stateJac[22] * covP[73]) + stateJac
                 [23] * covP[74]) + stateJac[24] * covP[75]) + stateJac[25] *
               covP[77]) + stateJac[26] * covP[78]) + stateJac[27] * covP[79];

  // 'updateCovPNoGps:375' covP(4, 6) = 1*tmp36;
  covP[93] = tmp36;

  // 'updateCovPNoGps:376' covP(4, 7) = 1*tmp38;
  covP[111] = tmp38;

  // 'updateCovPNoGps:377' covP(4, 8) = 1*tmp40;
  covP[129] = tmp40;

  // 'updateCovPNoGps:378' covP(4, 9) = 1*(covP(1, 9)*sJ4_1 + covP(2, 9)*sJ4_2 + covP(3, 9)*sJ4_3 + covP(4, 9)*sJ4_4 + covP(6, 9)*sJ4_6 + covP(7, 9)*sJ4_7 + covP(8, 9)*sJ4_8); 
  covP[147] = (((((stateJac[21] * covP[144] + stateJac[22] * covP[145]) +
                  stateJac[23] * covP[146]) + stateJac[24] * covP[147]) +
                stateJac[25] * covP[149]) + stateJac[26] * covP[150]) +
    stateJac[27] * covP[151];

  // 'updateCovPNoGps:379' covP(4, 10) = 1*(covP(1, 10)*sJ4_1 + covP(2, 10)*sJ4_2 + covP(3, 10)*sJ4_3 + covP(4, 10)*sJ4_4 + covP(6, 10)*sJ4_6 + covP(7, 10)*sJ4_7 + covP(8, 10)*sJ4_8); 
  covP[165] = (((((stateJac[21] * covP[162] + stateJac[22] * covP[163]) +
                  stateJac[23] * covP[164]) + stateJac[24] * covP[165]) +
                stateJac[25] * covP[167]) + stateJac[26] * covP[168]) +
    stateJac[27] * covP[169];

  // 'updateCovPNoGps:380' covP(4, 11) = 1*(covP(1, 11)*sJ4_1 + covP(2, 11)*sJ4_2 + covP(3, 11)*sJ4_3 + covP(4, 11)*sJ4_4 + covP(6, 11)*sJ4_6 + covP(7, 11)*sJ4_7 + covP(8, 11)*sJ4_8); 
  covP[183] = (((((stateJac[21] * covP[180] + stateJac[22] * covP[181]) +
                  stateJac[23] * covP[182]) + stateJac[24] * covP[183]) +
                stateJac[25] * covP[185]) + stateJac[26] * covP[186]) +
    stateJac[27] * covP[187];

  // 'updateCovPNoGps:381' covP(4, 12) = 1*(covP(1, 12)*sJ4_1 + covP(2, 12)*sJ4_2 + covP(3, 12)*sJ4_3 + covP(4, 12)*sJ4_4 + covP(6, 12)*sJ4_6 + covP(7, 12)*sJ4_7 + covP(8, 12)*sJ4_8); 
  covP[201] = (((((stateJac[21] * covP[198] + stateJac[22] * covP[199]) +
                  stateJac[23] * covP[200]) + stateJac[24] * covP[201]) +
                stateJac[25] * covP[203]) + stateJac[26] * covP[204]) +
    stateJac[27] * covP[205];

  // 'updateCovPNoGps:382' covP(4, 13) = 1*(covP(1, 13)*sJ4_1 + covP(2, 13)*sJ4_2 + covP(3, 13)*sJ4_3 + covP(4, 13)*sJ4_4 + covP(6, 13)*sJ4_6 + covP(7, 13)*sJ4_7 + covP(8, 13)*sJ4_8); 
  covP[219] = (((((stateJac[21] * covP[216] + stateJac[22] * covP[217]) +
                  stateJac[23] * covP[218]) + stateJac[24] * covP[219]) +
                stateJac[25] * covP[221]) + stateJac[26] * covP[222]) +
    stateJac[27] * covP[223];

  // 'updateCovPNoGps:383' covP(4, 14) = 1*(covP(1, 14)*sJ4_1 + covP(2, 14)*sJ4_2 + covP(3, 14)*sJ4_3 + covP(4, 14)*sJ4_4 + covP(6, 14)*sJ4_6 + covP(7, 14)*sJ4_7 + covP(8, 14)*sJ4_8); 
  covP[237] = (((((stateJac[21] * covP[234] + stateJac[22] * covP[235]) +
                  stateJac[23] * covP[236]) + stateJac[24] * covP[237]) +
                stateJac[25] * covP[239]) + stateJac[26] * covP[240]) +
    stateJac[27] * covP[241];

  // 'updateCovPNoGps:384' covP(4, 15) = 1*(covP(1, 15)*sJ4_1 + covP(2, 15)*sJ4_2 + covP(3, 15)*sJ4_3 + covP(4, 15)*sJ4_4 + covP(6, 15)*sJ4_6 + covP(7, 15)*sJ4_7 + covP(8, 15)*sJ4_8); 
  covP[255] = (((((stateJac[21] * covP[252] + stateJac[22] * covP[253]) +
                  stateJac[23] * covP[254]) + stateJac[24] * covP[255]) +
                stateJac[25] * covP[257]) + stateJac[26] * covP[258]) +
    stateJac[27] * covP[259];

  // 'updateCovPNoGps:385' covP(4, 16) = 1*(covP(1, 16)*sJ4_1 + covP(2, 16)*sJ4_2 + covP(3, 16)*sJ4_3 + covP(4, 16)*sJ4_4 + covP(6, 16)*sJ4_6 + covP(7, 16)*sJ4_7 + covP(8, 16)*sJ4_8); 
  covP[273] = (((((stateJac[21] * covP[270] + stateJac[22] * covP[271]) +
                  stateJac[23] * covP[272]) + stateJac[24] * covP[273]) +
                stateJac[25] * covP[275]) + stateJac[26] * covP[276]) +
    stateJac[27] * covP[277];

  // 'updateCovPNoGps:386' covP(4, 17) = 1*(covP(1, 17)*sJ4_1 + covP(2, 17)*sJ4_2 + covP(3, 17)*sJ4_3 + covP(4, 17)*sJ4_4 + covP(6, 17)*sJ4_6 + covP(7, 17)*sJ4_7 + covP(8, 17)*sJ4_8); 
  covP[291] = (((((stateJac[21] * covP[288] + stateJac[22] * covP[289]) +
                  stateJac[23] * covP[290]) + stateJac[24] * covP[291]) +
                stateJac[25] * covP[293]) + stateJac[26] * covP[294]) +
    stateJac[27] * covP[295];

  // 'updateCovPNoGps:387' covP(4, 18) = 1*(covP(1, 18)*sJ4_1 + covP(2, 18)*sJ4_2 + covP(3, 18)*sJ4_3 + covP(4, 18)*sJ4_4 + covP(6, 18)*sJ4_6 + covP(7, 18)*sJ4_7 + covP(8, 18)*sJ4_8); 
  covP[309] = (((((stateJac[21] * covP[306] + stateJac[22] * covP[307]) +
                  stateJac[23] * covP[308]) + stateJac[24] * covP[309]) +
                stateJac[25] * covP[311]) + stateJac[26] * covP[312]) +
    stateJac[27] * covP[313];

  // 'updateCovPNoGps:388' covP(5, 1) = sJ1_1*tmp41 + sJ1_2*tmp42 + sJ1_3*tmp43 + sJ1_4*tmp44 + sJ1_6*tmp45 + sJ1_7*tmp46 + sJ1_8*tmp47; 
  covP[4] = (((((stateJac[0] * tmp41 + stateJac[1] * tmp42) + stateJac[2] *
                tmp43) + stateJac[3] * tmp44) + stateJac[4] * tmp45) + stateJac
             [5] * tmp46) + stateJac[6] * tmp47;

  // 'updateCovPNoGps:389' covP(5, 2) = sJ2_1*tmp41 + sJ2_2*tmp42 + sJ2_3*tmp43 + sJ2_4*tmp44 + sJ2_6*tmp45 + sJ2_7*tmp46 + sJ2_8*tmp47; 
  covP[22] = (((((stateJac[7] * tmp41 + stateJac[8] * tmp42) + stateJac[9] *
                 tmp43) + stateJac[10] * tmp44) + stateJac[11] * tmp45) +
              stateJac[12] * tmp46) + stateJac[13] * tmp47;

  // 'updateCovPNoGps:390' covP(5, 3) = sJ3_1*tmp41 + sJ3_2*tmp42 + sJ3_3*tmp43 + sJ3_4*tmp44 + sJ3_6*tmp45 + sJ3_7*tmp46 + sJ3_8*tmp47; 
  covP[40] = (((((stateJac[14] * tmp41 + stateJac[15] * tmp42) + stateJac[16] *
                 tmp43) + stateJac[17] * tmp44) + stateJac[18] * tmp45) +
              stateJac[19] * tmp46) + stateJac[20] * tmp47;

  // 'updateCovPNoGps:391' covP(5, 4) = sJ4_1*tmp41 + sJ4_2*tmp42 + sJ4_3*tmp43 + sJ4_4*tmp44 + sJ4_6*tmp45 + sJ4_7*tmp46 + sJ4_8*tmp47; 
  covP[58] = (((((stateJac[21] * tmp41 + stateJac[22] * tmp42) + stateJac[23] *
                 tmp43) + stateJac[24] * tmp44) + stateJac[25] * tmp45) +
              stateJac[26] * tmp46) + stateJac[27] * tmp47;

  // 'updateCovPNoGps:392' covP(5, 5) = covP(5, 5)*1^2 + processNoiseQ(5, 5);
  covP[76] += processNoiseQ[76];

  // 'updateCovPNoGps:393' covP(5, 6) = 1*tmp45;
  covP[94] = tmp45;

  // 'updateCovPNoGps:394' covP(5, 7) = 1*tmp46;
  covP[112] = tmp46;

  // 'updateCovPNoGps:395' covP(5, 8) = 1*tmp47;
  covP[130] = tmp47;

  // 'updateCovPNoGps:396' covP(5, 9) = covP(5, 9)*tmp48;
  // 'updateCovPNoGps:397' covP(5, 10) = covP(5, 10)*tmp49;
  // 'updateCovPNoGps:398' covP(5, 11) = covP(5, 11)*tmp50;
  // 'updateCovPNoGps:399' covP(5, 12) = covP(5, 12)*tmp51;
  // 'updateCovPNoGps:400' covP(5, 13) = covP(5, 13)*tmp52;
  // 'updateCovPNoGps:401' covP(5, 14) = covP(5, 14)*tmp53;
  // 'updateCovPNoGps:402' covP(5, 15) = covP(5, 15)*tmp54;
  // 'updateCovPNoGps:403' covP(5, 16) = covP(5, 16)*tmp55;
  // 'updateCovPNoGps:404' covP(5, 17) = covP(5, 17)*tmp56;
  // 'updateCovPNoGps:405' covP(5, 18) = covP(5, 18)*tmp57;
  // 'updateCovPNoGps:406' covP(6, 1) = sJ1_1*tmp58 + sJ1_2*tmp59 + sJ1_3*tmp60 + sJ1_4*tmp61 + sJ1_7*tmp62 + sJ1_8*tmp63 + 1*tmp5; 
  covP[5] = (((((stateJac[0] * tmp58 + stateJac[1] * tmp59) + stateJac[2] *
                tmp60) + stateJac[3] * tmp61) + stateJac[5] * tmp62) + stateJac
             [6] * tmp63) + tmp5;

  // 'updateCovPNoGps:407' covP(6, 2) = sJ2_1*tmp58 + sJ2_2*tmp59 + sJ2_3*tmp60 + sJ2_4*tmp61 + sJ2_7*tmp62 + sJ2_8*tmp63 + 1*tmp15; 
  covP[23] = (((((stateJac[7] * tmp58 + stateJac[8] * tmp59) + stateJac[9] *
                 tmp60) + stateJac[10] * tmp61) + stateJac[12] * tmp62) +
              stateJac[13] * tmp63) + tmp15;

  // 'updateCovPNoGps:408' covP(6, 3) = sJ3_1*tmp58 + sJ3_2*tmp59 + sJ3_3*tmp60 + sJ3_4*tmp61 + sJ3_7*tmp62 + sJ3_8*tmp63 + 1*tmp25; 
  covP[41] = (((((stateJac[14] * tmp58 + stateJac[15] * tmp59) + stateJac[16] *
                 tmp60) + stateJac[17] * tmp61) + stateJac[19] * tmp62) +
              stateJac[20] * tmp63) + tmp25;

  // 'updateCovPNoGps:409' covP(6, 4) = sJ4_1*tmp58 + sJ4_2*tmp59 + sJ4_3*tmp60 + sJ4_4*tmp61 + sJ4_7*tmp62 + sJ4_8*tmp63 + 1*tmp35; 
  covP[59] = (((((stateJac[21] * tmp58 + stateJac[22] * tmp59) + stateJac[23] *
                 tmp60) + stateJac[24] * tmp61) + stateJac[26] * tmp62) +
              stateJac[27] * tmp63) + tmp35;

  // 'updateCovPNoGps:410' covP(6, 5) = covP(6, 5)*1*1;
  // 'updateCovPNoGps:411' covP(6, 6) = covP(6, 6)*1^2 + processNoiseQ(6, 6);
  covP[95] += processNoiseQ[95];

  // 'updateCovPNoGps:412' covP(6, 7) = 1*tmp62;
  covP[113] = tmp62;

  // 'updateCovPNoGps:413' covP(6, 8) = 1*tmp63;
  covP[131] = tmp63;

  // 'updateCovPNoGps:414' covP(6, 9) = covP(6, 9)*tmp64;
  // 'updateCovPNoGps:415' covP(6, 10) = covP(6, 10)*tmp65;
  // 'updateCovPNoGps:416' covP(6, 11) = covP(6, 11)*tmp66;
  // 'updateCovPNoGps:417' covP(6, 12) = covP(6, 12)*tmp67;
  // 'updateCovPNoGps:418' covP(6, 13) = covP(6, 13)*tmp68;
  // 'updateCovPNoGps:419' covP(6, 14) = covP(6, 14)*tmp69;
  // 'updateCovPNoGps:420' covP(6, 15) = covP(6, 15)*tmp70;
  // 'updateCovPNoGps:421' covP(6, 16) = covP(6, 16)*tmp71;
  // 'updateCovPNoGps:422' covP(6, 17) = covP(6, 17)*tmp72;
  // 'updateCovPNoGps:423' covP(6, 18) = covP(6, 18)*tmp73;
  // 'updateCovPNoGps:424' covP(7, 1) = sJ1_1*tmp74 + sJ1_2*tmp75 + sJ1_3*tmp76 + sJ1_4*tmp77 + sJ1_6*tmp78 + sJ1_8*tmp79 + 1*tmp7; 
  covP[6] = (((((stateJac[0] * tmp74 + stateJac[1] * tmp75) + stateJac[2] *
                tmp76) + stateJac[3] * tmp77) + stateJac[4] * tmp78) + stateJac
             [6] * tmp79) + tmp7;

  // 'updateCovPNoGps:425' covP(7, 2) = sJ2_1*tmp74 + sJ2_2*tmp75 + sJ2_3*tmp76 + sJ2_4*tmp77 + sJ2_6*tmp78 + sJ2_8*tmp79 + 1*tmp17; 
  covP[24] = (((((stateJac[7] * tmp74 + stateJac[8] * tmp75) + stateJac[9] *
                 tmp76) + stateJac[10] * tmp77) + stateJac[11] * tmp78) +
              stateJac[13] * tmp79) + tmp17;

  // 'updateCovPNoGps:426' covP(7, 3) = sJ3_1*tmp74 + sJ3_2*tmp75 + sJ3_3*tmp76 + sJ3_4*tmp77 + sJ3_6*tmp78 + sJ3_8*tmp79 + 1*tmp27; 
  covP[42] = (((((stateJac[14] * tmp74 + stateJac[15] * tmp75) + stateJac[16] *
                 tmp76) + stateJac[17] * tmp77) + stateJac[18] * tmp78) +
              stateJac[20] * tmp79) + tmp27;

  // 'updateCovPNoGps:427' covP(7, 4) = sJ4_1*tmp74 + sJ4_2*tmp75 + sJ4_3*tmp76 + sJ4_4*tmp77 + sJ4_6*tmp78 + sJ4_8*tmp79 + 1*tmp37; 
  covP[60] = (((((stateJac[21] * tmp74 + stateJac[22] * tmp75) + stateJac[23] *
                 tmp76) + stateJac[24] * tmp77) + stateJac[25] * tmp78) +
              stateJac[27] * tmp79) + tmp37;

  // 'updateCovPNoGps:428' covP(7, 5) = covP(7, 5)*1*1;
  // 'updateCovPNoGps:429' covP(7, 6) = 1*tmp78;
  covP[96] = tmp78;

  // 'updateCovPNoGps:430' covP(7, 7) = covP(7, 7)*1^2 + processNoiseQ(7, 7);
  covP[114] += processNoiseQ[114];

  // 'updateCovPNoGps:431' covP(7, 8) = 1*tmp79;
  covP[132] = tmp79;

  // 'updateCovPNoGps:432' covP(7, 9) = covP(7, 9)*tmp80;
  // 'updateCovPNoGps:433' covP(7, 10) = covP(7, 10)*tmp81;
  // 'updateCovPNoGps:434' covP(7, 11) = covP(7, 11)*tmp82;
  // 'updateCovPNoGps:435' covP(7, 12) = covP(7, 12)*tmp83;
  // 'updateCovPNoGps:436' covP(7, 13) = covP(7, 13)*tmp84;
  // 'updateCovPNoGps:437' covP(7, 14) = covP(7, 14)*tmp85;
  // 'updateCovPNoGps:438' covP(7, 15) = covP(7, 15)*tmp86;
  // 'updateCovPNoGps:439' covP(7, 16) = covP(7, 16)*tmp87;
  // 'updateCovPNoGps:440' covP(7, 17) = covP(7, 17)*tmp88;
  // 'updateCovPNoGps:441' covP(7, 18) = covP(7, 18)*tmp89;
  // 'updateCovPNoGps:442' covP(8, 1) = sJ1_1*tmp90 + sJ1_2*tmp91 + sJ1_3*tmp92 + sJ1_4*tmp93 + sJ1_6*tmp94 + sJ1_7*tmp95 + 1*tmp9; 
  covP[7] = (((((stateJac[0] * tmp90 + stateJac[1] * tmp91) + stateJac[2] *
                tmp92) + stateJac[3] * tmp93) + stateJac[4] * tmp94) + stateJac
             [5] * tmp95) + tmp9;

  // 'updateCovPNoGps:443' covP(8, 2) = sJ2_1*tmp90 + sJ2_2*tmp91 + sJ2_3*tmp92 + sJ2_4*tmp93 + sJ2_6*tmp94 + sJ2_7*tmp95 + 1*tmp19; 
  covP[25] = (((((stateJac[7] * tmp90 + stateJac[8] * tmp91) + stateJac[9] *
                 tmp92) + stateJac[10] * tmp93) + stateJac[11] * tmp94) +
              stateJac[12] * tmp95) + tmp19;

  // 'updateCovPNoGps:444' covP(8, 3) = sJ3_1*tmp90 + sJ3_2*tmp91 + sJ3_3*tmp92 + sJ3_4*tmp93 + sJ3_6*tmp94 + sJ3_7*tmp95 + 1*tmp29; 
  covP[43] = (((((stateJac[14] * tmp90 + stateJac[15] * tmp91) + stateJac[16] *
                 tmp92) + stateJac[17] * tmp93) + stateJac[18] * tmp94) +
              stateJac[19] * tmp95) + tmp29;

  // 'updateCovPNoGps:445' covP(8, 4) = sJ4_1*tmp90 + sJ4_2*tmp91 + sJ4_3*tmp92 + sJ4_4*tmp93 + sJ4_6*tmp94 + sJ4_7*tmp95 + 1*tmp39; 
  covP[61] = (((((stateJac[21] * tmp90 + stateJac[22] * tmp91) + stateJac[23] *
                 tmp92) + stateJac[24] * tmp93) + stateJac[25] * tmp94) +
              stateJac[26] * tmp95) + tmp39;

  // 'updateCovPNoGps:446' covP(8, 5) = covP(8, 5)*1*1;
  // 'updateCovPNoGps:447' covP(8, 6) = 1*tmp94;
  covP[97] = tmp94;

  // 'updateCovPNoGps:448' covP(8, 7) = 1*tmp95;
  covP[115] = tmp95;

  // 'updateCovPNoGps:449' covP(8, 8) = covP(8, 8)*1^2 + processNoiseQ(8, 8);
  covP[133] += processNoiseQ[133];

  // 'updateCovPNoGps:450' covP(8, 9) = covP(8, 9)*tmp96;
  // 'updateCovPNoGps:451' covP(8, 10) = covP(8, 10)*tmp97;
  // 'updateCovPNoGps:452' covP(8, 11) = covP(8, 11)*tmp98;
  // 'updateCovPNoGps:453' covP(8, 12) = covP(8, 12)*tmp99;
  // 'updateCovPNoGps:454' covP(8, 13) = covP(8, 13)*tmp100;
  // 'updateCovPNoGps:455' covP(8, 14) = covP(8, 14)*tmp101;
  // 'updateCovPNoGps:456' covP(8, 15) = covP(8, 15)*tmp102;
  // 'updateCovPNoGps:457' covP(8, 16) = covP(8, 16)*tmp103;
  // 'updateCovPNoGps:458' covP(8, 17) = covP(8, 17)*tmp104;
  // 'updateCovPNoGps:459' covP(8, 18) = covP(8, 18)*tmp105;
  // 'updateCovPNoGps:460' covP(9, 1) = sJ1_1*tmp106 + sJ1_2*tmp107 + sJ1_3*tmp108 + sJ1_4*tmp109 + sJ1_6*tmp110 + sJ1_7*tmp111 + sJ1_8*tmp112; 
  covP[8] = (((((stateJac[0] * tmp106 + stateJac[1] * tmp107) + stateJac[2] *
                tmp108) + stateJac[3] * tmp109) + stateJac[4] * tmp110) +
             stateJac[5] * tmp111) + stateJac[6] * tmp112;

  // 'updateCovPNoGps:461' covP(9, 2) = sJ2_1*tmp106 + sJ2_2*tmp107 + sJ2_3*tmp108 + sJ2_4*tmp109 + sJ2_6*tmp110 + sJ2_7*tmp111 + sJ2_8*tmp112; 
  covP[26] = (((((stateJac[7] * tmp106 + stateJac[8] * tmp107) + stateJac[9] *
                 tmp108) + stateJac[10] * tmp109) + stateJac[11] * tmp110) +
              stateJac[12] * tmp111) + stateJac[13] * tmp112;

  // 'updateCovPNoGps:462' covP(9, 3) = sJ3_1*tmp106 + sJ3_2*tmp107 + sJ3_3*tmp108 + sJ3_4*tmp109 + sJ3_6*tmp110 + sJ3_7*tmp111 + sJ3_8*tmp112; 
  covP[44] = (((((stateJac[14] * tmp106 + stateJac[15] * tmp107) + stateJac[16] *
                 tmp108) + stateJac[17] * tmp109) + stateJac[18] * tmp110) +
              stateJac[19] * tmp111) + stateJac[20] * tmp112;

  // 'updateCovPNoGps:463' covP(9, 4) = sJ4_1*tmp106 + sJ4_2*tmp107 + sJ4_3*tmp108 + sJ4_4*tmp109 + sJ4_6*tmp110 + sJ4_7*tmp111 + sJ4_8*tmp112; 
  covP[62] = (((((stateJac[21] * tmp106 + stateJac[22] * tmp107) + stateJac[23] *
                 tmp108) + stateJac[24] * tmp109) + stateJac[25] * tmp110) +
              stateJac[26] * tmp111) + stateJac[27] * tmp112;

  // 'updateCovPNoGps:464' covP(9, 5) = covP(9, 5)*tmp48;
  // 'updateCovPNoGps:465' covP(9, 6) = covP(9, 6)*tmp64;
  // 'updateCovPNoGps:466' covP(9, 7) = covP(9, 7)*tmp80;
  // 'updateCovPNoGps:467' covP(9, 8) = covP(9, 8)*tmp96;
  // 'updateCovPNoGps:468' covP(9, 9) = covP(9, 9)*1^2 + processNoiseQ(9, 9);
  covP[152] += processNoiseQ[152];

  // 'updateCovPNoGps:469' covP(9, 10) = covP(9, 10)*tmp113;
  // 'updateCovPNoGps:470' covP(9, 11) = covP(9, 11)*tmp114;
  // 'updateCovPNoGps:471' covP(9, 12) = covP(9, 12)*tmp115;
  // 'updateCovPNoGps:472' covP(9, 13) = covP(9, 13)*tmp116;
  // 'updateCovPNoGps:473' covP(9, 14) = covP(9, 14)*tmp117;
  // 'updateCovPNoGps:474' covP(9, 15) = covP(9, 15)*tmp118;
  // 'updateCovPNoGps:475' covP(9, 16) = covP(9, 16)*tmp119;
  // 'updateCovPNoGps:476' covP(9, 17) = covP(9, 17)*tmp120;
  // 'updateCovPNoGps:477' covP(9, 18) = covP(9, 18)*tmp121;
  // 'updateCovPNoGps:478' covP(10, 1) = sJ1_1*tmp122 + sJ1_2*tmp123 + sJ1_3*tmp124 + sJ1_4*tmp125 + sJ1_6*tmp126 + sJ1_7*tmp127 + sJ1_8*tmp128; 
  covP[9] = (((((stateJac[0] * tmp122 + stateJac[1] * tmp123) + stateJac[2] *
                tmp124) + stateJac[3] * tmp125) + stateJac[4] * tmp126) +
             stateJac[5] * tmp127) + stateJac[6] * tmp128;

  // 'updateCovPNoGps:479' covP(10, 2) = sJ2_1*tmp122 + sJ2_2*tmp123 + sJ2_3*tmp124 + sJ2_4*tmp125 + sJ2_6*tmp126 + sJ2_7*tmp127 + sJ2_8*tmp128; 
  covP[27] = (((((stateJac[7] * tmp122 + stateJac[8] * tmp123) + stateJac[9] *
                 tmp124) + stateJac[10] * tmp125) + stateJac[11] * tmp126) +
              stateJac[12] * tmp127) + stateJac[13] * tmp128;

  // 'updateCovPNoGps:480' covP(10, 3) = sJ3_1*tmp122 + sJ3_2*tmp123 + sJ3_3*tmp124 + sJ3_4*tmp125 + sJ3_6*tmp126 + sJ3_7*tmp127 + sJ3_8*tmp128; 
  covP[45] = (((((stateJac[14] * tmp122 + stateJac[15] * tmp123) + stateJac[16] *
                 tmp124) + stateJac[17] * tmp125) + stateJac[18] * tmp126) +
              stateJac[19] * tmp127) + stateJac[20] * tmp128;

  // 'updateCovPNoGps:481' covP(10, 4) = sJ4_1*tmp122 + sJ4_2*tmp123 + sJ4_3*tmp124 + sJ4_4*tmp125 + sJ4_6*tmp126 + sJ4_7*tmp127 + sJ4_8*tmp128; 
  covP[63] = (((((stateJac[21] * tmp122 + stateJac[22] * tmp123) + stateJac[23] *
                 tmp124) + stateJac[24] * tmp125) + stateJac[25] * tmp126) +
              stateJac[26] * tmp127) + stateJac[27] * tmp128;

  // 'updateCovPNoGps:482' covP(10, 5) = covP(10, 5)*tmp49;
  // 'updateCovPNoGps:483' covP(10, 6) = covP(10, 6)*tmp65;
  // 'updateCovPNoGps:484' covP(10, 7) = covP(10, 7)*tmp81;
  // 'updateCovPNoGps:485' covP(10, 8) = covP(10, 8)*tmp97;
  // 'updateCovPNoGps:486' covP(10, 9) = covP(10, 9)*tmp113;
  // 'updateCovPNoGps:487' covP(10, 10) = covP(10, 10)*1^2 + processNoiseQ(10, 10); 
  covP[171] += processNoiseQ[171];

  // 'updateCovPNoGps:488' covP(10, 11) = covP(10, 11)*tmp129;
  // 'updateCovPNoGps:489' covP(10, 12) = covP(10, 12)*tmp130;
  // 'updateCovPNoGps:490' covP(10, 13) = covP(10, 13)*tmp131;
  // 'updateCovPNoGps:491' covP(10, 14) = covP(10, 14)*tmp132;
  // 'updateCovPNoGps:492' covP(10, 15) = covP(10, 15)*tmp133;
  // 'updateCovPNoGps:493' covP(10, 16) = covP(10, 16)*tmp134;
  // 'updateCovPNoGps:494' covP(10, 17) = covP(10, 17)*tmp135;
  // 'updateCovPNoGps:495' covP(10, 18) = covP(10, 18)*tmp136;
  // 'updateCovPNoGps:496' covP(11, 1) = sJ1_1*tmp137 + sJ1_2*tmp138 + sJ1_3*tmp139 + sJ1_4*tmp140 + sJ1_6*tmp141 + sJ1_7*tmp142 + sJ1_8*tmp143; 
  covP[10] = (((((stateJac[0] * tmp137 + stateJac[1] * tmp138) + stateJac[2] *
                 tmp139) + stateJac[3] * tmp140) + stateJac[4] * tmp141) +
              stateJac[5] * tmp142) + stateJac[6] * tmp143;

  // 'updateCovPNoGps:497' covP(11, 2) = sJ2_1*tmp137 + sJ2_2*tmp138 + sJ2_3*tmp139 + sJ2_4*tmp140 + sJ2_6*tmp141 + sJ2_7*tmp142 + sJ2_8*tmp143; 
  covP[28] = (((((stateJac[7] * tmp137 + stateJac[8] * tmp138) + stateJac[9] *
                 tmp139) + stateJac[10] * tmp140) + stateJac[11] * tmp141) +
              stateJac[12] * tmp142) + stateJac[13] * tmp143;

  // 'updateCovPNoGps:498' covP(11, 3) = sJ3_1*tmp137 + sJ3_2*tmp138 + sJ3_3*tmp139 + sJ3_4*tmp140 + sJ3_6*tmp141 + sJ3_7*tmp142 + sJ3_8*tmp143; 
  covP[46] = (((((stateJac[14] * tmp137 + stateJac[15] * tmp138) + stateJac[16] *
                 tmp139) + stateJac[17] * tmp140) + stateJac[18] * tmp141) +
              stateJac[19] * tmp142) + stateJac[20] * tmp143;

  // 'updateCovPNoGps:499' covP(11, 4) = sJ4_1*tmp137 + sJ4_2*tmp138 + sJ4_3*tmp139 + sJ4_4*tmp140 + sJ4_6*tmp141 + sJ4_7*tmp142 + sJ4_8*tmp143; 
  covP[64] = (((((stateJac[21] * tmp137 + stateJac[22] * tmp138) + stateJac[23] *
                 tmp139) + stateJac[24] * tmp140) + stateJac[25] * tmp141) +
              stateJac[26] * tmp142) + stateJac[27] * tmp143;

  // 'updateCovPNoGps:500' covP(11, 5) = covP(11, 5)*tmp50;
  // 'updateCovPNoGps:501' covP(11, 6) = covP(11, 6)*tmp66;
  // 'updateCovPNoGps:502' covP(11, 7) = covP(11, 7)*tmp82;
  // 'updateCovPNoGps:503' covP(11, 8) = covP(11, 8)*tmp98;
  // 'updateCovPNoGps:504' covP(11, 9) = covP(11, 9)*tmp114;
  // 'updateCovPNoGps:505' covP(11, 10) = covP(11, 10)*tmp129;
  // 'updateCovPNoGps:506' covP(11, 11) = covP(11, 11)*1^2 + processNoiseQ(11, 11); 
  covP[190] += processNoiseQ[190];

  // 'updateCovPNoGps:507' covP(11, 12) = covP(11, 12)*tmp144;
  // 'updateCovPNoGps:508' covP(11, 13) = covP(11, 13)*tmp145;
  // 'updateCovPNoGps:509' covP(11, 14) = covP(11, 14)*tmp146;
  // 'updateCovPNoGps:510' covP(11, 15) = covP(11, 15)*tmp147;
  // 'updateCovPNoGps:511' covP(11, 16) = covP(11, 16)*tmp148;
  // 'updateCovPNoGps:512' covP(11, 17) = covP(11, 17)*tmp149;
  // 'updateCovPNoGps:513' covP(11, 18) = covP(11, 18)*tmp150;
  // 'updateCovPNoGps:514' covP(12, 1) = sJ1_1*tmp151 + sJ1_2*tmp152 + sJ1_3*tmp153 + sJ1_4*tmp154 + sJ1_6*tmp155 + sJ1_7*tmp156 + sJ1_8*tmp157; 
  covP[11] = (((((stateJac[0] * tmp151 + stateJac[1] * tmp152) + stateJac[2] *
                 tmp153) + stateJac[3] * tmp154) + stateJac[4] * tmp155) +
              stateJac[5] * tmp156) + stateJac[6] * tmp157;

  // 'updateCovPNoGps:515' covP(12, 2) = sJ2_1*tmp151 + sJ2_2*tmp152 + sJ2_3*tmp153 + sJ2_4*tmp154 + sJ2_6*tmp155 + sJ2_7*tmp156 + sJ2_8*tmp157; 
  covP[29] = (((((stateJac[7] * tmp151 + stateJac[8] * tmp152) + stateJac[9] *
                 tmp153) + stateJac[10] * tmp154) + stateJac[11] * tmp155) +
              stateJac[12] * tmp156) + stateJac[13] * tmp157;

  // 'updateCovPNoGps:516' covP(12, 3) = sJ3_1*tmp151 + sJ3_2*tmp152 + sJ3_3*tmp153 + sJ3_4*tmp154 + sJ3_6*tmp155 + sJ3_7*tmp156 + sJ3_8*tmp157; 
  covP[47] = (((((stateJac[14] * tmp151 + stateJac[15] * tmp152) + stateJac[16] *
                 tmp153) + stateJac[17] * tmp154) + stateJac[18] * tmp155) +
              stateJac[19] * tmp156) + stateJac[20] * tmp157;

  // 'updateCovPNoGps:517' covP(12, 4) = sJ4_1*tmp151 + sJ4_2*tmp152 + sJ4_3*tmp153 + sJ4_4*tmp154 + sJ4_6*tmp155 + sJ4_7*tmp156 + sJ4_8*tmp157; 
  covP[65] = (((((stateJac[21] * tmp151 + stateJac[22] * tmp152) + stateJac[23] *
                 tmp153) + stateJac[24] * tmp154) + stateJac[25] * tmp155) +
              stateJac[26] * tmp156) + stateJac[27] * tmp157;

  // 'updateCovPNoGps:518' covP(12, 5) = covP(12, 5)*tmp51;
  // 'updateCovPNoGps:519' covP(12, 6) = covP(12, 6)*tmp67;
  // 'updateCovPNoGps:520' covP(12, 7) = covP(12, 7)*tmp83;
  // 'updateCovPNoGps:521' covP(12, 8) = covP(12, 8)*tmp99;
  // 'updateCovPNoGps:522' covP(12, 9) = covP(12, 9)*tmp115;
  // 'updateCovPNoGps:523' covP(12, 10) = covP(12, 10)*tmp130;
  // 'updateCovPNoGps:524' covP(12, 11) = covP(12, 11)*tmp144;
  // 'updateCovPNoGps:525' covP(12, 12) = covP(12, 12)*1^2 + processNoiseQ(12, 12); 
  covP[209] += processNoiseQ[209];

  // 'updateCovPNoGps:526' covP(12, 13) = covP(12, 13)*tmp158;
  // 'updateCovPNoGps:527' covP(12, 14) = covP(12, 14)*tmp159;
  // 'updateCovPNoGps:528' covP(12, 15) = covP(12, 15)*tmp160;
  // 'updateCovPNoGps:529' covP(12, 16) = covP(12, 16)*tmp161;
  // 'updateCovPNoGps:530' covP(12, 17) = covP(12, 17)*tmp162;
  // 'updateCovPNoGps:531' covP(12, 18) = covP(12, 18)*tmp163;
  // 'updateCovPNoGps:532' covP(13, 1) = sJ1_1*tmp164 + sJ1_2*tmp165 + sJ1_3*tmp166 + sJ1_4*tmp167 + sJ1_6*tmp168 + sJ1_7*tmp169 + sJ1_8*tmp170; 
  covP[12] = (((((stateJac[0] * tmp164 + stateJac[1] * tmp165) + stateJac[2] *
                 tmp166) + stateJac[3] * tmp167) + stateJac[4] * tmp168) +
              stateJac[5] * tmp169) + stateJac[6] * tmp170;

  // 'updateCovPNoGps:533' covP(13, 2) = sJ2_1*tmp164 + sJ2_2*tmp165 + sJ2_3*tmp166 + sJ2_4*tmp167 + sJ2_6*tmp168 + sJ2_7*tmp169 + sJ2_8*tmp170; 
  covP[30] = (((((stateJac[7] * tmp164 + stateJac[8] * tmp165) + stateJac[9] *
                 tmp166) + stateJac[10] * tmp167) + stateJac[11] * tmp168) +
              stateJac[12] * tmp169) + stateJac[13] * tmp170;

  // 'updateCovPNoGps:534' covP(13, 3) = sJ3_1*tmp164 + sJ3_2*tmp165 + sJ3_3*tmp166 + sJ3_4*tmp167 + sJ3_6*tmp168 + sJ3_7*tmp169 + sJ3_8*tmp170; 
  covP[48] = (((((stateJac[14] * tmp164 + stateJac[15] * tmp165) + stateJac[16] *
                 tmp166) + stateJac[17] * tmp167) + stateJac[18] * tmp168) +
              stateJac[19] * tmp169) + stateJac[20] * tmp170;

  // 'updateCovPNoGps:535' covP(13, 4) = sJ4_1*tmp164 + sJ4_2*tmp165 + sJ4_3*tmp166 + sJ4_4*tmp167 + sJ4_6*tmp168 + sJ4_7*tmp169 + sJ4_8*tmp170; 
  covP[66] = (((((stateJac[21] * tmp164 + stateJac[22] * tmp165) + stateJac[23] *
                 tmp166) + stateJac[24] * tmp167) + stateJac[25] * tmp168) +
              stateJac[26] * tmp169) + stateJac[27] * tmp170;

  // 'updateCovPNoGps:536' covP(13, 5) = covP(13, 5)*tmp52;
  // 'updateCovPNoGps:537' covP(13, 6) = covP(13, 6)*tmp68;
  // 'updateCovPNoGps:538' covP(13, 7) = covP(13, 7)*tmp84;
  // 'updateCovPNoGps:539' covP(13, 8) = covP(13, 8)*tmp100;
  // 'updateCovPNoGps:540' covP(13, 9) = covP(13, 9)*tmp116;
  // 'updateCovPNoGps:541' covP(13, 10) = covP(13, 10)*tmp131;
  // 'updateCovPNoGps:542' covP(13, 11) = covP(13, 11)*tmp145;
  // 'updateCovPNoGps:543' covP(13, 12) = covP(13, 12)*tmp158;
  // 'updateCovPNoGps:544' covP(13, 13) = covP(13, 13)*1^2 + processNoiseQ(13, 13); 
  covP[228] += processNoiseQ[228];

  // 'updateCovPNoGps:545' covP(13, 14) = covP(13, 14)*tmp171;
  // 'updateCovPNoGps:546' covP(13, 15) = covP(13, 15)*tmp172;
  // 'updateCovPNoGps:547' covP(13, 16) = covP(13, 16)*tmp173;
  // 'updateCovPNoGps:548' covP(13, 17) = covP(13, 17)*tmp174;
  // 'updateCovPNoGps:549' covP(13, 18) = covP(13, 18)*tmp175;
  // 'updateCovPNoGps:550' covP(14, 1) = sJ1_1*tmp176 + sJ1_2*tmp177 + sJ1_3*tmp178 + sJ1_4*tmp179 + sJ1_6*tmp180 + sJ1_7*tmp181 + sJ1_8*tmp182; 
  covP[13] = (((((stateJac[0] * tmp176 + stateJac[1] * tmp177) + stateJac[2] *
                 tmp178) + stateJac[3] * tmp179) + stateJac[4] * tmp180) +
              stateJac[5] * tmp181) + stateJac[6] * tmp182;

  // 'updateCovPNoGps:551' covP(14, 2) = sJ2_1*tmp176 + sJ2_2*tmp177 + sJ2_3*tmp178 + sJ2_4*tmp179 + sJ2_6*tmp180 + sJ2_7*tmp181 + sJ2_8*tmp182; 
  covP[31] = (((((stateJac[7] * tmp176 + stateJac[8] * tmp177) + stateJac[9] *
                 tmp178) + stateJac[10] * tmp179) + stateJac[11] * tmp180) +
              stateJac[12] * tmp181) + stateJac[13] * tmp182;

  // 'updateCovPNoGps:552' covP(14, 3) = sJ3_1*tmp176 + sJ3_2*tmp177 + sJ3_3*tmp178 + sJ3_4*tmp179 + sJ3_6*tmp180 + sJ3_7*tmp181 + sJ3_8*tmp182; 
  covP[49] = (((((stateJac[14] * tmp176 + stateJac[15] * tmp177) + stateJac[16] *
                 tmp178) + stateJac[17] * tmp179) + stateJac[18] * tmp180) +
              stateJac[19] * tmp181) + stateJac[20] * tmp182;

  // 'updateCovPNoGps:553' covP(14, 4) = sJ4_1*tmp176 + sJ4_2*tmp177 + sJ4_3*tmp178 + sJ4_4*tmp179 + sJ4_6*tmp180 + sJ4_7*tmp181 + sJ4_8*tmp182; 
  covP[67] = (((((stateJac[21] * tmp176 + stateJac[22] * tmp177) + stateJac[23] *
                 tmp178) + stateJac[24] * tmp179) + stateJac[25] * tmp180) +
              stateJac[26] * tmp181) + stateJac[27] * tmp182;

  // 'updateCovPNoGps:554' covP(14, 5) = covP(14, 5)*tmp53;
  // 'updateCovPNoGps:555' covP(14, 6) = covP(14, 6)*tmp69;
  // 'updateCovPNoGps:556' covP(14, 7) = covP(14, 7)*tmp85;
  // 'updateCovPNoGps:557' covP(14, 8) = covP(14, 8)*tmp101;
  // 'updateCovPNoGps:558' covP(14, 9) = covP(14, 9)*tmp117;
  // 'updateCovPNoGps:559' covP(14, 10) = covP(14, 10)*tmp132;
  // 'updateCovPNoGps:560' covP(14, 11) = covP(14, 11)*tmp146;
  // 'updateCovPNoGps:561' covP(14, 12) = covP(14, 12)*tmp159;
  // 'updateCovPNoGps:562' covP(14, 13) = covP(14, 13)*tmp171;
  // 'updateCovPNoGps:563' covP(14, 14) = covP(14, 14)*1^2 + processNoiseQ(14, 14); 
  covP[247] += processNoiseQ[247];

  // 'updateCovPNoGps:564' covP(14, 15) = covP(14, 15)*tmp183;
  // 'updateCovPNoGps:565' covP(14, 16) = covP(14, 16)*tmp184;
  // 'updateCovPNoGps:566' covP(14, 17) = covP(14, 17)*tmp185;
  // 'updateCovPNoGps:567' covP(14, 18) = covP(14, 18)*tmp186;
  // 'updateCovPNoGps:568' covP(15, 1) = sJ1_1*tmp187 + sJ1_2*tmp188 + sJ1_3*tmp189 + sJ1_4*tmp190 + sJ1_6*tmp191 + sJ1_7*tmp192 + sJ1_8*tmp193; 
  covP[14] = (((((stateJac[0] * tmp187 + stateJac[1] * tmp188) + stateJac[2] *
                 tmp189) + stateJac[3] * tmp190) + stateJac[4] * tmp191) +
              stateJac[5] * tmp192) + stateJac[6] * tmp193;

  // 'updateCovPNoGps:569' covP(15, 2) = sJ2_1*tmp187 + sJ2_2*tmp188 + sJ2_3*tmp189 + sJ2_4*tmp190 + sJ2_6*tmp191 + sJ2_7*tmp192 + sJ2_8*tmp193; 
  covP[32] = (((((stateJac[7] * tmp187 + stateJac[8] * tmp188) + stateJac[9] *
                 tmp189) + stateJac[10] * tmp190) + stateJac[11] * tmp191) +
              stateJac[12] * tmp192) + stateJac[13] * tmp193;

  // 'updateCovPNoGps:570' covP(15, 3) = sJ3_1*tmp187 + sJ3_2*tmp188 + sJ3_3*tmp189 + sJ3_4*tmp190 + sJ3_6*tmp191 + sJ3_7*tmp192 + sJ3_8*tmp193; 
  covP[50] = (((((stateJac[14] * tmp187 + stateJac[15] * tmp188) + stateJac[16] *
                 tmp189) + stateJac[17] * tmp190) + stateJac[18] * tmp191) +
              stateJac[19] * tmp192) + stateJac[20] * tmp193;

  // 'updateCovPNoGps:571' covP(15, 4) = sJ4_1*tmp187 + sJ4_2*tmp188 + sJ4_3*tmp189 + sJ4_4*tmp190 + sJ4_6*tmp191 + sJ4_7*tmp192 + sJ4_8*tmp193; 
  covP[68] = (((((stateJac[21] * tmp187 + stateJac[22] * tmp188) + stateJac[23] *
                 tmp189) + stateJac[24] * tmp190) + stateJac[25] * tmp191) +
              stateJac[26] * tmp192) + stateJac[27] * tmp193;

  // 'updateCovPNoGps:572' covP(15, 5) = covP(15, 5)*tmp54;
  // 'updateCovPNoGps:573' covP(15, 6) = covP(15, 6)*tmp70;
  // 'updateCovPNoGps:574' covP(15, 7) = covP(15, 7)*tmp86;
  // 'updateCovPNoGps:575' covP(15, 8) = covP(15, 8)*tmp102;
  // 'updateCovPNoGps:576' covP(15, 9) = covP(15, 9)*tmp118;
  // 'updateCovPNoGps:577' covP(15, 10) = covP(15, 10)*tmp133;
  // 'updateCovPNoGps:578' covP(15, 11) = covP(15, 11)*tmp147;
  // 'updateCovPNoGps:579' covP(15, 12) = covP(15, 12)*tmp160;
  // 'updateCovPNoGps:580' covP(15, 13) = covP(15, 13)*tmp172;
  // 'updateCovPNoGps:581' covP(15, 14) = covP(15, 14)*tmp183;
  // 'updateCovPNoGps:582' covP(15, 15) = covP(15, 15)*1^2 + processNoiseQ(15, 15); 
  covP[266] += processNoiseQ[266];

  // 'updateCovPNoGps:583' covP(15, 16) = covP(15, 16)*tmp194;
  // 'updateCovPNoGps:584' covP(15, 17) = covP(15, 17)*tmp195;
  // 'updateCovPNoGps:585' covP(15, 18) = covP(15, 18)*tmp196;
  // 'updateCovPNoGps:586' covP(16, 1) = sJ1_1*tmp197 + sJ1_2*tmp198 + sJ1_3*tmp199 + sJ1_4*tmp200 + sJ1_6*tmp201 + sJ1_7*tmp202 + sJ1_8*tmp203; 
  covP[15] = (((((stateJac[0] * tmp197 + stateJac[1] * tmp198) + stateJac[2] *
                 tmp199) + stateJac[3] * tmp200) + stateJac[4] * tmp201) +
              stateJac[5] * tmp202) + stateJac[6] * tmp203;

  // 'updateCovPNoGps:587' covP(16, 2) = sJ2_1*tmp197 + sJ2_2*tmp198 + sJ2_3*tmp199 + sJ2_4*tmp200 + sJ2_6*tmp201 + sJ2_7*tmp202 + sJ2_8*tmp203; 
  covP[33] = (((((stateJac[7] * tmp197 + stateJac[8] * tmp198) + stateJac[9] *
                 tmp199) + stateJac[10] * tmp200) + stateJac[11] * tmp201) +
              stateJac[12] * tmp202) + stateJac[13] * tmp203;

  // 'updateCovPNoGps:588' covP(16, 3) = sJ3_1*tmp197 + sJ3_2*tmp198 + sJ3_3*tmp199 + sJ3_4*tmp200 + sJ3_6*tmp201 + sJ3_7*tmp202 + sJ3_8*tmp203; 
  covP[51] = (((((stateJac[14] * tmp197 + stateJac[15] * tmp198) + stateJac[16] *
                 tmp199) + stateJac[17] * tmp200) + stateJac[18] * tmp201) +
              stateJac[19] * tmp202) + stateJac[20] * tmp203;

  // 'updateCovPNoGps:589' covP(16, 4) = sJ4_1*tmp197 + sJ4_2*tmp198 + sJ4_3*tmp199 + sJ4_4*tmp200 + sJ4_6*tmp201 + sJ4_7*tmp202 + sJ4_8*tmp203; 
  covP[69] = (((((stateJac[21] * tmp197 + stateJac[22] * tmp198) + stateJac[23] *
                 tmp199) + stateJac[24] * tmp200) + stateJac[25] * tmp201) +
              stateJac[26] * tmp202) + stateJac[27] * tmp203;

  // 'updateCovPNoGps:590' covP(16, 5) = covP(16, 5)*tmp55;
  // 'updateCovPNoGps:591' covP(16, 6) = covP(16, 6)*tmp71;
  // 'updateCovPNoGps:592' covP(16, 7) = covP(16, 7)*tmp87;
  // 'updateCovPNoGps:593' covP(16, 8) = covP(16, 8)*tmp103;
  // 'updateCovPNoGps:594' covP(16, 9) = covP(16, 9)*tmp119;
  // 'updateCovPNoGps:595' covP(16, 10) = covP(16, 10)*tmp134;
  // 'updateCovPNoGps:596' covP(16, 11) = covP(16, 11)*tmp148;
  // 'updateCovPNoGps:597' covP(16, 12) = covP(16, 12)*tmp161;
  // 'updateCovPNoGps:598' covP(16, 13) = covP(16, 13)*tmp173;
  // 'updateCovPNoGps:599' covP(16, 14) = covP(16, 14)*tmp184;
  // 'updateCovPNoGps:600' covP(16, 15) = covP(16, 15)*tmp194;
  // 'updateCovPNoGps:601' covP(16, 16) = covP(16, 16)*1^2 + processNoiseQ(16, 16); 
  covP[285] += processNoiseQ[285];

  // 'updateCovPNoGps:602' covP(16, 17) = covP(16, 17)*tmp204;
  // 'updateCovPNoGps:603' covP(16, 18) = covP(16, 18)*tmp205;
  // 'updateCovPNoGps:604' covP(17, 1) = sJ1_1*tmp206 + sJ1_2*tmp207 + sJ1_3*tmp208 + sJ1_4*tmp209 + sJ1_6*tmp210 + sJ1_7*tmp211 + sJ1_8*tmp212; 
  covP[16] = (((((stateJac[0] * tmp206 + stateJac[1] * tmp207) + stateJac[2] *
                 tmp208) + stateJac[3] * tmp209) + stateJac[4] * tmp210) +
              stateJac[5] * tmp211) + stateJac[6] * tmp212;

  // 'updateCovPNoGps:605' covP(17, 2) = sJ2_1*tmp206 + sJ2_2*tmp207 + sJ2_3*tmp208 + sJ2_4*tmp209 + sJ2_6*tmp210 + sJ2_7*tmp211 + sJ2_8*tmp212; 
  covP[34] = (((((stateJac[7] * tmp206 + stateJac[8] * tmp207) + stateJac[9] *
                 tmp208) + stateJac[10] * tmp209) + stateJac[11] * tmp210) +
              stateJac[12] * tmp211) + stateJac[13] * tmp212;

  // 'updateCovPNoGps:606' covP(17, 3) = sJ3_1*tmp206 + sJ3_2*tmp207 + sJ3_3*tmp208 + sJ3_4*tmp209 + sJ3_6*tmp210 + sJ3_7*tmp211 + sJ3_8*tmp212; 
  covP[52] = (((((stateJac[14] * tmp206 + stateJac[15] * tmp207) + stateJac[16] *
                 tmp208) + stateJac[17] * tmp209) + stateJac[18] * tmp210) +
              stateJac[19] * tmp211) + stateJac[20] * tmp212;

  // 'updateCovPNoGps:607' covP(17, 4) = sJ4_1*tmp206 + sJ4_2*tmp207 + sJ4_3*tmp208 + sJ4_4*tmp209 + sJ4_6*tmp210 + sJ4_7*tmp211 + sJ4_8*tmp212; 
  covP[70] = (((((stateJac[21] * tmp206 + stateJac[22] * tmp207) + stateJac[23] *
                 tmp208) + stateJac[24] * tmp209) + stateJac[25] * tmp210) +
              stateJac[26] * tmp211) + stateJac[27] * tmp212;

  // 'updateCovPNoGps:608' covP(17, 5) = covP(17, 5)*tmp56;
  // 'updateCovPNoGps:609' covP(17, 6) = covP(17, 6)*tmp72;
  // 'updateCovPNoGps:610' covP(17, 7) = covP(17, 7)*tmp88;
  // 'updateCovPNoGps:611' covP(17, 8) = covP(17, 8)*tmp104;
  // 'updateCovPNoGps:612' covP(17, 9) = covP(17, 9)*tmp120;
  // 'updateCovPNoGps:613' covP(17, 10) = covP(17, 10)*tmp135;
  // 'updateCovPNoGps:614' covP(17, 11) = covP(17, 11)*tmp149;
  // 'updateCovPNoGps:615' covP(17, 12) = covP(17, 12)*tmp162;
  // 'updateCovPNoGps:616' covP(17, 13) = covP(17, 13)*tmp174;
  // 'updateCovPNoGps:617' covP(17, 14) = covP(17, 14)*tmp185;
  // 'updateCovPNoGps:618' covP(17, 15) = covP(17, 15)*tmp195;
  // 'updateCovPNoGps:619' covP(17, 16) = covP(17, 16)*tmp204;
  // 'updateCovPNoGps:620' covP(17, 17) = covP(17, 17)*1^2 + processNoiseQ(17, 17); 
  covP[304] += processNoiseQ[304];

  // 'updateCovPNoGps:621' covP(17, 18) = covP(17, 18)*tmp213;
  // 'updateCovPNoGps:622' covP(18, 1) = sJ1_1*tmp214 + sJ1_2*tmp215 + sJ1_3*tmp216 + sJ1_4*tmp217 + sJ1_6*tmp218 + sJ1_7*tmp219 + sJ1_8*tmp220; 
  covP[17] = (((((stateJac[0] * tmp214 + stateJac[1] * tmp215) + stateJac[2] *
                 tmp216) + stateJac[3] * tmp217) + stateJac[4] * tmp218) +
              stateJac[5] * tmp219) + stateJac[6] * tmp220;

  // 'updateCovPNoGps:623' covP(18, 2) = sJ2_1*tmp214 + sJ2_2*tmp215 + sJ2_3*tmp216 + sJ2_4*tmp217 + sJ2_6*tmp218 + sJ2_7*tmp219 + sJ2_8*tmp220; 
  covP[35] = (((((stateJac[7] * tmp214 + stateJac[8] * tmp215) + stateJac[9] *
                 tmp216) + stateJac[10] * tmp217) + stateJac[11] * tmp218) +
              stateJac[12] * tmp219) + stateJac[13] * tmp220;

  // 'updateCovPNoGps:624' covP(18, 3) = sJ3_1*tmp214 + sJ3_2*tmp215 + sJ3_3*tmp216 + sJ3_4*tmp217 + sJ3_6*tmp218 + sJ3_7*tmp219 + sJ3_8*tmp220; 
  covP[53] = (((((stateJac[14] * tmp214 + stateJac[15] * tmp215) + stateJac[16] *
                 tmp216) + stateJac[17] * tmp217) + stateJac[18] * tmp218) +
              stateJac[19] * tmp219) + stateJac[20] * tmp220;

  // 'updateCovPNoGps:625' covP(18, 4) = sJ4_1*tmp214 + sJ4_2*tmp215 + sJ4_3*tmp216 + sJ4_4*tmp217 + sJ4_6*tmp218 + sJ4_7*tmp219 + sJ4_8*tmp220; 
  covP[71] = (((((stateJac[21] * tmp214 + stateJac[22] * tmp215) + stateJac[23] *
                 tmp216) + stateJac[24] * tmp217) + stateJac[25] * tmp218) +
              stateJac[26] * tmp219) + stateJac[27] * tmp220;

  // 'updateCovPNoGps:626' covP(18, 5) = covP(18, 5)*tmp57;
  // 'updateCovPNoGps:627' covP(18, 6) = covP(18, 6)*tmp73;
  // 'updateCovPNoGps:628' covP(18, 7) = covP(18, 7)*tmp89;
  // 'updateCovPNoGps:629' covP(18, 8) = covP(18, 8)*tmp105;
  // 'updateCovPNoGps:630' covP(18, 9) = covP(18, 9)*tmp121;
  // 'updateCovPNoGps:631' covP(18, 10) = covP(18, 10)*tmp136;
  // 'updateCovPNoGps:632' covP(18, 11) = covP(18, 11)*tmp150;
  // 'updateCovPNoGps:633' covP(18, 12) = covP(18, 12)*tmp163;
  // 'updateCovPNoGps:634' covP(18, 13) = covP(18, 13)*tmp175;
  // 'updateCovPNoGps:635' covP(18, 14) = covP(18, 14)*tmp186;
  // 'updateCovPNoGps:636' covP(18, 15) = covP(18, 15)*tmp196;
  // 'updateCovPNoGps:637' covP(18, 16) = covP(18, 16)*tmp205;
  // 'updateCovPNoGps:638' covP(18, 17) = covP(18, 17)*tmp213;
  // 'updateCovPNoGps:639' covP(18, 18) = covP(18, 18)*1^2 + processNoiseQ(18, 18); 
  covP[323] += processNoiseQ[323];
}

//
// File trailer for generated code.
//
// [EOF]
//
