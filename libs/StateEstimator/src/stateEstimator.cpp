//
// File: stateEstimator.cpp
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM Cortex-M
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. ROM efficiency
// Validation result: Not run
//
#include "stateEstimator.h"
#include "stateEstimator_types.h"
#include "rtwtypes.h"
#include <cstring>
#include <cmath>
#include "norm_7MzYkgry.h"
#include "computeStateJac_YJuzBPAK.h"
#include "updateCovP_UFEbdQNU.h"
#include "mrdiv_0UppHIxu.h"
#include "mrdiv_8PldARhB.h"
#include "updateCovPNoGps_vX01OT1j.h"
#include "mrdiv_2CyhSJN4.h"
#include "computeMagMeasJacNoGps_mcxOPfaA.h"
#include "stateEstimator_private.h"

// Named constants for Chart: '<Root>/estimatorStateMachine'
const uint8_T stateEstima_IN_RUN_GPS_NOT_INIT{ 4U };

const uint8_T stateEstimator_IN_INITIALIZE{ 1U };

const uint8_T stateEstimator_IN_RUN{ 2U };

const uint8_T stateEstimator_IN_RUN_GPS_LOST{ 3U };

const uint8_T stateEstimator_IN_RUN_INIT_GPS{ 5U };

// Function for Chart: '<Root>/estimatorStateMachine'
void stateEstimator::stateEstimator_INITIALIZE(enumStateEstimateMode *mode,
  boolean_T *isMagValid, boolean_T *isGpsValid, real32_T *baroAltOut_m,
  boolean_T *isBaroValid, real32_T bodyAccelsOut_mps2[3], const busMagData
  *rtu_magData, const busGpsData *rtu_gpsData, const busBaroData *rtu_baroData,
  const busStateEstSmParams *rtu_stateEstSmParams)
{
  *mode = enumStateEstimateMode::INITIALIZE;
  stateEstimator_DW.resetStates = true;

  // During 'INITIALIZE': '<S5>:1'
  // '<S5>:44:1' sf_internal_predicateOutput = isAttInitialized && isPosInitialized && isBaroInitialized; 
  if (static_cast<boolean_T>(static_cast<boolean_T>
       (stateEstimator_DW.isAttInitialized & stateEstimator_DW.isPosInitialized)
       & stateEstimator_DW.isBaroInitialized)) {
    // Transition: '<S5>:44'
    stateEstimator_DW.durationCounter_1 = 0;
    stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_RUN;

    // Entry 'RUN': '<S5>:43'
    // FULL EKF WITH GPS IS RUNNING
    // '<S5>:43:4' resetStates = false;
    stateEstimator_DW.resetStates = false;

    // '<S5>:43:5' mode = enumStateEstimateMode.RUN;
    *mode = enumStateEstimateMode::RUN;

    // Chart: '<Root>/estimatorStateMachine'
    // '<S5>:43:6' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
    // '<S5>:43:7' bodyRatesOut_radps = filtBodyRatesIn_radps;
    // '<S5>:43:8' normMagVecOut_nd = normMagVecIn_nd;
    // '<S5>:43:9' isMagValid = isMagDataValid;
    *isMagValid = rtu_magData->isMagDataValid;

    // Product: '<S4>/Product'
    // '<S5>:43:10' latLonAltOut = latLonAltIn;
    bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];

    // SignalConversion generated from: '<S5>/ SFunction '
    stateEstimator_DW.bodyRatesOut_radps[0] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];

    // Product: '<S9>/Divide'
    stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];

    // Chart: '<Root>/estimatorStateMachine'
    stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];

    // Product: '<S4>/Product'
    bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];

    // SignalConversion generated from: '<S5>/ SFunction '
    stateEstimator_DW.bodyRatesOut_radps[1] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];

    // Product: '<S9>/Divide'
    stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];

    // Chart: '<Root>/estimatorStateMachine'
    stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];

    // Product: '<S4>/Product'
    bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];

    // SignalConversion generated from: '<S5>/ SFunction '
    stateEstimator_DW.bodyRatesOut_radps[2] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];

    // Product: '<S9>/Divide'
    stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];

    // Chart: '<Root>/estimatorStateMachine'
    stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

    // '<S5>:43:11' isGpsValid = isGpsDataValid;
    *isGpsValid = rtu_gpsData->isGpsDataValid;

    // '<S5>:43:12' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
    *baroAltOut_m = stateEstimator_DW.Divide1 -
      stateEstimator_DW.baroInitAltMean;

    // Chart: '<Root>/estimatorStateMachine'
    // '<S5>:43:13' isBaroValid = isBaroDataValid;
    *isBaroValid = rtu_baroData->isBaroDataValid;

    //

    // Chart: '<Root>/estimatorStateMachine' incorporates:
    //   Product: '<S4>/Product'
    //   Product: '<S9>/Divide'
    //   SignalConversion generated from: '<S5>/ SFunction '

    // '<S5>:63:1' sf_internal_predicateOutput = isAttInitialized && ~isGpsInitialized &&  isBaroInitialized; 
  } else if (static_cast<boolean_T>(static_cast<boolean_T>(static_cast<boolean_T>
               (rtu_gpsData->isGpsInitialized ^ 1) &
               stateEstimator_DW.isAttInitialized) &
              stateEstimator_DW.isBaroInitialized)) {
    // Transition: '<S5>:63'
    stateEstimator_DW.durationCounter_1_p = 0;
    stateEstimator_DW.is_c3_stateEstimator = stateEstima_IN_RUN_GPS_NOT_INIT;

    // Entry 'RUN_GPS_NOT_INIT': '<S5>:62'
    // EKF STARTED RUNNING WITHOUT GPS
    // '<S5>:62:4' resetStates = false;
    stateEstimator_DW.resetStates = false;

    // '<S5>:62:5' isPosInitialized = false;
    stateEstimator_DW.isPosInitialized = false;

    // '<S5>:62:6' gpsValidCount = 0;
    stateEstimator_DW.gpsValidCount = 0U;

    // '<S5>:62:7' mode = enumStateEstimateMode.RUN_GPS_NOT_INIT;
    *mode = enumStateEstimateMode::RUN_GPS_NOT_INIT;

    // '<S5>:62:8' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
    // '<S5>:62:9' bodyRatesOut_radps = filtBodyRatesIn_radps;
    // '<S5>:62:10' normMagVecOut_nd = normMagVecIn_nd;
    // '<S5>:62:11' isMagValid = isMagDataValid;
    *isMagValid = rtu_magData->isMagDataValid;

    // '<S5>:62:12' latLonAltOut = latLonAltIn;
    bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
    stateEstimator_DW.bodyRatesOut_radps[0] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
    stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
    stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
    bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
    stateEstimator_DW.bodyRatesOut_radps[1] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
    stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
    stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
    bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
    stateEstimator_DW.bodyRatesOut_radps[2] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
    stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
    stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

    // '<S5>:62:13' isGpsValid = isGpsDataValid;
    *isGpsValid = rtu_gpsData->isGpsDataValid;

    // '<S5>:62:14' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
    *baroAltOut_m = stateEstimator_DW.Divide1 -
      stateEstimator_DW.baroInitAltMean;

    // '<S5>:62:15' isBaroValid = isBaroDataValid;
    *isBaroValid = rtu_baroData->isBaroDataValid;

    //
  } else {
    real32_T imuDelta_idx_0;
    real32_T imuDelta_idx_1;
    real32_T imuDelta_idx_2;
    real32_T imuDelta_idx_3;
    real32_T imuDelta_idx_4;
    real32_T imuDelta_idx_5;

    // Take first n(user specified) data to find the initial attitude, accel and gyro biases and NED 
    // origin lat, lon and alt
    //
    // Compute running mean and bias of IMU data
    // '<S5>:1:48' if ( imuIdx < max(stateEstSmParams.imuInitCount, 1) )
    if (stateEstimator_DW.imuIdx < std::fmax(rtu_stateEstSmParams->imuInitCount,
         1.0F)) {
      // '<S5>:1:49' imuIdx = imuIdx + 1;
      stateEstimator_DW.imuIdx++;

      // '<S5>:1:50' imuDelta = [filtBodyAccelsIn_mps2; filtBodyRatesIn_radps] -  ... 
      // '<S5>:1:51'         imuMean;
      imuDelta_idx_0 = stateEstimator_DW.Product[0] - stateEstimator_DW.imuMean
        [0];
      imuDelta_idx_3 = stateEstimator_DW.TmpSignalConversionAtSFunctionI[0] -
        stateEstimator_DW.imuMean[3];
      imuDelta_idx_1 = stateEstimator_DW.Product[1] - stateEstimator_DW.imuMean
        [1];
      imuDelta_idx_4 = stateEstimator_DW.TmpSignalConversionAtSFunctionI[1] -
        stateEstimator_DW.imuMean[4];
      imuDelta_idx_2 = stateEstimator_DW.Product[2] - stateEstimator_DW.imuMean
        [2];
      imuDelta_idx_5 = stateEstimator_DW.TmpSignalConversionAtSFunctionI[2] -
        stateEstimator_DW.imuMean[5];

      // '<S5>:1:52' imuMean = imuMean + imuDelta / imuIdx;
      stateEstimator_DW.imuMean[0] += imuDelta_idx_0 / stateEstimator_DW.imuIdx;
      stateEstimator_DW.imuMean[1] += imuDelta_idx_1 / stateEstimator_DW.imuIdx;
      stateEstimator_DW.imuMean[2] += imuDelta_idx_2 / stateEstimator_DW.imuIdx;
      stateEstimator_DW.imuMean[3] += imuDelta_idx_3 / stateEstimator_DW.imuIdx;
      stateEstimator_DW.imuMean[4] += imuDelta_idx_4 / stateEstimator_DW.imuIdx;
      stateEstimator_DW.imuMean[5] += imuDelta_idx_5 / stateEstimator_DW.imuIdx;

      // '<S5>:1:53' imuM2 = imuM2 + imuDelta .* ( [filtBodyAccelsIn_mps2; filtBodyRatesIn_radps] - ... 
      // '<S5>:1:54'         imuMean);
      stateEstimator_DW.imuM2[0] += (stateEstimator_DW.Product[0] -
        stateEstimator_DW.imuMean[0]) * imuDelta_idx_0;
      stateEstimator_DW.imuM2[1] += (stateEstimator_DW.Product[1] -
        stateEstimator_DW.imuMean[1]) * imuDelta_idx_1;
      stateEstimator_DW.imuM2[2] += (stateEstimator_DW.Product[2] -
        stateEstimator_DW.imuMean[2]) * imuDelta_idx_2;
      stateEstimator_DW.imuM2[3] +=
        (stateEstimator_DW.TmpSignalConversionAtSFunctionI[0] -
         stateEstimator_DW.imuMean[3]) * imuDelta_idx_3;
      stateEstimator_DW.imuM2[4] +=
        (stateEstimator_DW.TmpSignalConversionAtSFunctionI[1] -
         stateEstimator_DW.imuMean[4]) * imuDelta_idx_4;
      stateEstimator_DW.imuM2[5] +=
        (stateEstimator_DW.TmpSignalConversionAtSFunctionI[2] -
         stateEstimator_DW.imuMean[5]) * imuDelta_idx_5;
    }

    //
    // '<S5>:1:57' if (imuIdx >= stateEstSmParams.imuInitCount)
    if (stateEstimator_DW.imuIdx >= rtu_stateEstSmParams->imuInitCount) {
      // '<S5>:1:58' if (stateEstSmParams.imuInitCount > 1)
      if (rtu_stateEstSmParams->imuInitCount > 1.0F) {
        // compute accel and gyro biases
        // '<S5>:1:60' accelBias_mps2 = sqrt( imuM2(1 : 3) / (stateEstSmParams.imuInitCount - 1) ); 
        // '<S5>:1:61' gyroBias_radps = sqrt( imuM2(4 : 6) / (stateEstSmParams.imuInitCount - 1) ); 
        stateEstimator_DW.accelBias_mps2[0] = std::sqrt(stateEstimator_DW.imuM2
          [0] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
        stateEstimator_DW.gyroBias_radps[0] = std::sqrt(stateEstimator_DW.imuM2
          [3] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
        stateEstimator_DW.accelBias_mps2[1] = std::sqrt(stateEstimator_DW.imuM2
          [1] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
        stateEstimator_DW.gyroBias_radps[1] = std::sqrt(stateEstimator_DW.imuM2
          [4] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
        stateEstimator_DW.accelBias_mps2[2] = std::sqrt(stateEstimator_DW.imuM2
          [2] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
        stateEstimator_DW.gyroBias_radps[2] = std::sqrt(stateEstimator_DW.imuM2
          [5] / (rtu_stateEstSmParams->imuInitCount - 1.0F));
      } else {
        // '<S5>:1:62' else
        // '<S5>:1:63' accelBias_mps2 = [0; 0; 0];
        // '<S5>:1:64' gyroBias_radps = [0; 0; 0];
        stateEstimator_DW.accelBias_mps2[0] = 0.0F;
        stateEstimator_DW.gyroBias_radps[0] = 0.0F;
        stateEstimator_DW.accelBias_mps2[1] = 0.0F;
        stateEstimator_DW.gyroBias_radps[1] = 0.0F;
        stateEstimator_DW.accelBias_mps2[2] = 0.0F;
        stateEstimator_DW.gyroBias_radps[2] = 0.0F;
      }
    }

    //
    // Compute running mean and bias of MAG data
    // '<S5>:1:69' if (isMagDataValid)
    if (rtu_magData->isMagDataValid) {
      // '<S5>:1:70' if( magIdx < max(stateEstSmParams.magInitCount, 1) )
      if (stateEstimator_DW.magIdx < std::fmax
          (rtu_stateEstSmParams->magInitCount, 1.0F)) {
        // '<S5>:1:71' magIdx = magIdx + 1;
        stateEstimator_DW.magIdx++;

        // '<S5>:1:72' magDelta = (normMagVecIn_nd - magMean);
        // '<S5>:1:73' magMean = magMean + magDelta / magIdx;
        // '<S5>:1:74' magM2 = magM2 + magDelta .* (normMagVecIn_nd - magMean);
        imuDelta_idx_0 = stateEstimator_DW.Divide[0] -
          stateEstimator_DW.magMean[0];
        stateEstimator_DW.magMean[0] += imuDelta_idx_0 /
          stateEstimator_DW.magIdx;
        stateEstimator_DW.magM2[0] += (stateEstimator_DW.Divide[0] -
          stateEstimator_DW.magMean[0]) * imuDelta_idx_0;
        imuDelta_idx_0 = stateEstimator_DW.Divide[1] -
          stateEstimator_DW.magMean[1];
        stateEstimator_DW.magMean[1] += imuDelta_idx_0 /
          stateEstimator_DW.magIdx;
        stateEstimator_DW.magM2[1] += (stateEstimator_DW.Divide[1] -
          stateEstimator_DW.magMean[1]) * imuDelta_idx_0;
        imuDelta_idx_0 = stateEstimator_DW.Divide[2] -
          stateEstimator_DW.magMean[2];
        stateEstimator_DW.magMean[2] += imuDelta_idx_0 /
          stateEstimator_DW.magIdx;
        stateEstimator_DW.magM2[2] += (stateEstimator_DW.Divide[2] -
          stateEstimator_DW.magMean[2]) * imuDelta_idx_0;
      }

      //
      // '<S5>:1:77' if(magIdx >= stateEstSmParams.magInitCount)
      if (stateEstimator_DW.magIdx >= rtu_stateEstSmParams->magInitCount) {
        real32_T nedMagVecNorm_nd_tmp;
        real32_T sHalfPsi;
        real32_T scale;
        real32_T yaw_rad_tmp;

        // '<S5>:1:78' if(stateEstSmParams.magInitCount > 1)
        if (rtu_stateEstSmParams->magInitCount > 1.0F) {
          // compute mag bias
          // '<S5>:1:80' magBias_nd = sqrt( magM2 / (stateEstSmParams.magInitCount - 1) ); 
          stateEstimator_DW.magBias_nd[0] = std::sqrt(stateEstimator_DW.magM2[0]
            / (rtu_stateEstSmParams->magInitCount - 1.0F));
          stateEstimator_DW.magBias_nd[1] = std::sqrt(stateEstimator_DW.magM2[1]
            / (rtu_stateEstSmParams->magInitCount - 1.0F));
          stateEstimator_DW.magBias_nd[2] = std::sqrt(stateEstimator_DW.magM2[2]
            / (rtu_stateEstSmParams->magInitCount - 1.0F));
        } else {
          // '<S5>:1:81' else
          // '<S5>:1:82' magBias_nd = [0; 0; 0];
          stateEstimator_DW.magBias_nd[0] = 0.0F;
          stateEstimator_DW.magBias_nd[1] = 0.0F;
          stateEstimator_DW.magBias_nd[2] = 0.0F;
        }

        // '<S5>:1:84' [initialQuat, nedMagVecNorm_nd] = computeInitialAttitude(imuMean(1:3), magMean, stateEstSmParams.initMagDec_rad); 
        // COMPUTEINITIALATTITUDE Computes the IMU attitude using IMU and mag data 
        //
        // Inputs:
        // bodyAccels_mps2:    3x1 array of sum of accel_[x, y, z] readings
        // magVecNorm_nd:      3x1 unit vector obtained from magnetometer strapped to 
        // the body
        // magDec_rad:         Magnetic Declination
        //
        // Outputs:
        // quat:               4x1 array of attitude represented as
        // quaterion with first entry being the real entry
        // nedMagVecNorm_nd:   Mag unit vector in NED frame
        // 'computeInitialAttitude:14' roll_rad = atan2( -bodyAccels_mps2(2), -bodyAccels_mps2(3) ); 
        imuDelta_idx_4 = std::atan2(-stateEstimator_DW.imuMean[1],
          -stateEstimator_DW.imuMean[2]);

        // 'computeInitialAttitude:15' pitch_rad = atan2( -bodyAccels_mps2(1), norm(bodyAccels_mps2(2:3)) ); 
        scale = 1.29246971E-26F;
        imuDelta_idx_2 = std::abs(stateEstimator_DW.imuMean[1]);
        if (imuDelta_idx_2 > 1.29246971E-26F) {
          imuDelta_idx_0 = 1.0F;
          scale = imuDelta_idx_2;
        } else {
          imuDelta_idx_5 = imuDelta_idx_2 / 1.29246971E-26F;
          imuDelta_idx_0 = imuDelta_idx_5 * imuDelta_idx_5;
        }

        imuDelta_idx_2 = std::abs(stateEstimator_DW.imuMean[2]);
        if (imuDelta_idx_2 > scale) {
          imuDelta_idx_5 = scale / imuDelta_idx_2;
          imuDelta_idx_0 = imuDelta_idx_0 * imuDelta_idx_5 * imuDelta_idx_5 +
            1.0F;
          scale = imuDelta_idx_2;
        } else {
          imuDelta_idx_5 = imuDelta_idx_2 / scale;
          imuDelta_idx_0 += imuDelta_idx_5 * imuDelta_idx_5;
        }

        imuDelta_idx_2 = std::atan2(-stateEstimator_DW.imuMean[0], scale * std::
          sqrt(imuDelta_idx_0));

        //  Compute corrected magnetometer readings
        // 'computeInitialAttitude:18' cPhi = cos(roll_rad);
        imuDelta_idx_0 = std::cos(imuDelta_idx_4);

        // 'computeInitialAttitude:19' sPhi = sin(roll_rad);
        scale = std::sin(imuDelta_idx_4);

        // 'computeInitialAttitude:20' sTheta = sin(pitch_rad);
        imuDelta_idx_5 = std::sin(imuDelta_idx_2);

        // 'computeInitialAttitude:21' cTheta = cos(pitch_rad);
        // 'computeInitialAttitude:22' yaw_rad = atan2(-magVecNorm_nd(2) * cPhi + magVecNorm_nd(3) * sPhi, ... 
        // 'computeInitialAttitude:23'                     magVecNorm_nd(1) * cTheta + magVecNorm_nd(2) * sPhi * sTheta + ... 
        // 'computeInitialAttitude:24'                     magVecNorm_nd(3) * cPhi * sTheta) + magDec_rad; 
        imuDelta_idx_3 = std::atan2(-stateEstimator_DW.magMean[1] *
          imuDelta_idx_0 + stateEstimator_DW.magMean[2] * scale,
          (stateEstimator_DW.magMean[1] * scale * imuDelta_idx_5 +
           stateEstimator_DW.magMean[0] * std::cos(imuDelta_idx_2)) +
          stateEstimator_DW.magMean[2] * imuDelta_idx_0 * imuDelta_idx_5) +
          rtu_stateEstSmParams->initMagDec_rad;

        //  Compute half-angles
        // 'computeInitialAttitude:27' cHalfPhi= cos(roll_rad/2);
        imuDelta_idx_0 = std::cos(imuDelta_idx_4 / 2.0F);

        // 'computeInitialAttitude:28' sHalfPhi = sin(roll_rad/2);
        imuDelta_idx_4 = std::sin(imuDelta_idx_4 / 2.0F);

        // 'computeInitialAttitude:32' cHalfTheta = cos(pitch_rad/2);
        scale = std::cos(imuDelta_idx_2 / 2.0F);

        // 'computeInitialAttitude:33' sHalfTheta = sin(pitch_rad/2);
        imuDelta_idx_2 = std::sin(imuDelta_idx_2 / 2.0F);

        // 'computeInitialAttitude:35' cHalfPsi = cos(yaw_rad/2);
        imuDelta_idx_5 = std::cos(imuDelta_idx_3 / 2.0F);

        // 'computeInitialAttitude:36' sHalfPsi = sin(yaw_rad/2);
        sHalfPsi = std::sin(imuDelta_idx_3 / 2.0F);

        //  Compute quaternion components
        // 'computeInitialAttitude:39' q0 = cHalfPsi * cHalfTheta * cHalfPhi + sHalfPsi * sHalfTheta * sHalfPhi; 
        imuDelta_idx_1 = imuDelta_idx_5 * scale;
        yaw_rad_tmp = sHalfPsi * imuDelta_idx_2;
        imuDelta_idx_3 = imuDelta_idx_1 * imuDelta_idx_0 + yaw_rad_tmp *
          imuDelta_idx_4;

        // 'computeInitialAttitude:40' q1 = cHalfPsi * cHalfTheta * sHalfPhi - sHalfPsi * sHalfTheta * cHalfPhi; 
        imuDelta_idx_1 = imuDelta_idx_1 * imuDelta_idx_4 - yaw_rad_tmp *
          imuDelta_idx_0;

        // 'computeInitialAttitude:41' q2 = cHalfPsi * sHalfTheta * cHalfPhi + sHalfPsi * cHalfTheta * sHalfPhi; 
        scale *= sHalfPsi;
        imuDelta_idx_5 *= imuDelta_idx_2;
        imuDelta_idx_2 = imuDelta_idx_5 * imuDelta_idx_0 + scale *
          imuDelta_idx_4;

        // 'computeInitialAttitude:42' q3 = sHalfPsi * cHalfTheta * cHalfPhi - cHalfPsi * sHalfTheta * sHalfPhi; 
        imuDelta_idx_0 = scale * imuDelta_idx_0 - imuDelta_idx_5 *
          imuDelta_idx_4;

        // 'computeInitialAttitude:43' quat = [q0; q1; q2; q3];
        //  Direction Cosine Matrix (DCM) from body cooridinates to NED coordinates 
        //  expressed using quaternions.
        // 'computeInitialAttitude:47' C_b2ned=[1-2*(q2^2+q3^2), 2*(q1*q2-q3*q0), 2*(q1*q3+q2*q0); 
        // 'computeInitialAttitude:48'          2*(q1*q2+q3*q0), 1-2*(q1^2+q3^2), 2*(q2*q3-q1*q0); 
        // 'computeInitialAttitude:49'          2*(q1*q3-q2*q0), 2*(q2*q3+q1*q0), 1-2*(q1^2+q2^2)]; 
        // 'computeInitialAttitude:50' nedMagVecNorm_nd = C_b2ned * magVecNorm_nd; 
        stateEstimator_DW.initialQuat[0] = imuDelta_idx_3;
        stateEstimator_DW.initialQuat[1] = imuDelta_idx_1;
        stateEstimator_DW.initialQuat[2] = imuDelta_idx_2;
        stateEstimator_DW.initialQuat[3] = imuDelta_idx_0;
        imuDelta_idx_4 = imuDelta_idx_0 * imuDelta_idx_0;
        imuDelta_idx_5 = imuDelta_idx_2 * imuDelta_idx_2;
        stateEstimator_DW.nedMagVecNorm_nd[0] = (1.0F - (imuDelta_idx_5 +
          imuDelta_idx_4) * 2.0F) * stateEstimator_DW.magMean[0];
        scale = imuDelta_idx_1 * imuDelta_idx_2;
        sHalfPsi = imuDelta_idx_0 * imuDelta_idx_3;
        stateEstimator_DW.nedMagVecNorm_nd[0] += (scale - sHalfPsi) * 2.0F *
          stateEstimator_DW.magMean[1];
        yaw_rad_tmp = imuDelta_idx_1 * imuDelta_idx_0;
        nedMagVecNorm_nd_tmp = imuDelta_idx_2 * imuDelta_idx_3;
        stateEstimator_DW.nedMagVecNorm_nd[0] += (yaw_rad_tmp +
          nedMagVecNorm_nd_tmp) * 2.0F * stateEstimator_DW.magMean[2];
        stateEstimator_DW.nedMagVecNorm_nd[1] = (scale + sHalfPsi) * 2.0F *
          stateEstimator_DW.magMean[0];
        scale = imuDelta_idx_1 * imuDelta_idx_1;
        stateEstimator_DW.nedMagVecNorm_nd[1] += (1.0F - (scale + imuDelta_idx_4)
          * 2.0F) * stateEstimator_DW.magMean[1];
        imuDelta_idx_4 = imuDelta_idx_2 * imuDelta_idx_0;
        sHalfPsi = imuDelta_idx_1 * imuDelta_idx_3;
        stateEstimator_DW.nedMagVecNorm_nd[1] += (imuDelta_idx_4 - sHalfPsi) *
          2.0F * stateEstimator_DW.magMean[2];
        stateEstimator_DW.nedMagVecNorm_nd[2] = (yaw_rad_tmp -
          nedMagVecNorm_nd_tmp) * 2.0F * stateEstimator_DW.magMean[0];
        stateEstimator_DW.nedMagVecNorm_nd[2] += (imuDelta_idx_4 + sHalfPsi) *
          2.0F * stateEstimator_DW.magMean[1];
        stateEstimator_DW.nedMagVecNorm_nd[2] += (1.0F - (scale + imuDelta_idx_5)
          * 2.0F) * stateEstimator_DW.magMean[2];

        // '<S5>:1:85' isAttInitialized = true;
        stateEstimator_DW.isAttInitialized = true;
      }
    }

    //
    // Compute running mean of GPS data for NED origin Lat, Lon and Alt
    // '<S5>:1:90' if (isGpsDataValid && isGpsInitialized)
    if (static_cast<boolean_T>(rtu_gpsData->isGpsDataValid &
         rtu_gpsData->isGpsInitialized)) {
      // '<S5>:1:91' if( gpsIdx < max(stateEstSmParams.gpsInitCount, 1) )
      if (stateEstimator_DW.gpsIdx < std::fmax
          (rtu_stateEstSmParams->gpsInitCount, 1.0F)) {
        // '<S5>:1:92' gpsIdx = gpsIdx + 1;
        stateEstimator_DW.gpsIdx++;

        // '<S5>:1:93' llhDelta = (latLonAltIn - refLatLonAlt);
        // '<S5>:1:94' refLatLonAlt = refLatLonAlt + llhDelta/gpsIdx;
        stateEstimator_DW.refLatLonAlt[0] += (rtu_gpsData->latLonAlt[0] -
          stateEstimator_DW.refLatLonAlt[0]) / stateEstimator_DW.gpsIdx;
        stateEstimator_DW.refLatLonAlt[1] += (rtu_gpsData->latLonAlt[1] -
          stateEstimator_DW.refLatLonAlt[1]) / stateEstimator_DW.gpsIdx;
        stateEstimator_DW.refLatLonAlt[2] += (rtu_gpsData->latLonAlt[2] -
          stateEstimator_DW.refLatLonAlt[2]) / stateEstimator_DW.gpsIdx;
      }

      //
      // '<S5>:1:97' if(gpsIdx >= stateEstSmParams.gpsInitCount)
      if (stateEstimator_DW.gpsIdx >= rtu_stateEstSmParams->gpsInitCount) {
        // '<S5>:1:98' isPosInitialized = true;
        stateEstimator_DW.isPosInitialized = true;
      }
    }

    //
    // Compute running mean and bias of baro data
    // '<S5>:1:103' if (isBaroDataValid)
    if (rtu_baroData->isBaroDataValid) {
      // '<S5>:1:104' if( baroIdx < max(stateEstSmParams.baroInitCount, 1) )
      if (stateEstimator_DW.baroIdx < std::fmax
          (rtu_stateEstSmParams->baroInitCount, 1.0F)) {
        // '<S5>:1:105' baroIdx = baroIdx + 1;
        stateEstimator_DW.baroIdx++;

        // '<S5>:1:106' baroInitAltDelta = (baroPressAlt_m - baroInitAltMean);
        imuDelta_idx_0 = stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean;

        // '<S5>:1:107' baroInitAltMean = baroInitAltMean + baroInitAltDelta / baroIdx; 
        stateEstimator_DW.baroInitAltMean += imuDelta_idx_0 /
          stateEstimator_DW.baroIdx;

        // '<S5>:1:108' baroInitAltM2 = baroInitAltM2 + baroInitAltDelta .* (baroPressAlt_m - baroInitAltMean); 
        stateEstimator_DW.baroInitAltM2 += (stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean) * imuDelta_idx_0;
      }

      //
      // '<S5>:1:111' if(baroIdx >= stateEstSmParams.baroInitCount)
      if (stateEstimator_DW.baroIdx >= rtu_stateEstSmParams->baroInitCount) {
        // '<S5>:1:112' if(stateEstSmParams.baroInitCount > 1)
        if (rtu_stateEstSmParams->baroInitCount > 1.0F) {
          // compute baro bias
          // '<S5>:1:114' baroBias_m = sqrt( baroInitAltM2 / (stateEstSmParams.baroInitCount - 1) ); 
          stateEstimator_DW.baroBias_m = stateEstimator_DW.baroInitAltM2 /
            (rtu_stateEstSmParams->baroInitCount - 1.0F);
          stateEstimator_DW.baroBias_m = std::sqrt(stateEstimator_DW.baroBias_m);
        } else {
          // '<S5>:1:115' else
          // '<S5>:1:116' baroBias_m = 0;
          stateEstimator_DW.baroBias_m = 0.0F;
        }

        // '<S5>:1:118' isBaroInitialized = true;
        stateEstimator_DW.isBaroInitialized = true;
      }
    }

    // '<S5>:1:121' if(isGpsInitialized)
    if (rtu_gpsData->isGpsInitialized) {
      // Status of the initialization
      // '<S5>:1:123' stateEstInitPct = (imuIdx + magIdx + gpsIdx + baroIdx) / (stateEstSmParams.imuInitCount + ... 
      // '<S5>:1:124'         stateEstSmParams.magInitCount + stateEstSmParams.gpsInitCount +  ... 
      // '<S5>:1:125'         stateEstSmParams.baroInitCount) * 100;
      stateEstimator_DW.stateEstInitPct = (((stateEstimator_DW.imuIdx +
        stateEstimator_DW.magIdx) + static_cast<real32_T>
        (stateEstimator_DW.gpsIdx)) + stateEstimator_DW.baroIdx) /
        (((rtu_stateEstSmParams->imuInitCount +
           rtu_stateEstSmParams->magInitCount) +
          rtu_stateEstSmParams->gpsInitCount) +
         rtu_stateEstSmParams->baroInitCount) * 100.0F;
    } else {
      // '<S5>:1:126' else
      // Status of the initialization
      // '<S5>:1:128' stateEstInitPct = (imuIdx + magIdx + baroIdx) / (stateEstSmParams.imuInitCount + ... 
      // '<S5>:1:129'         stateEstSmParams.magInitCount + stateEstSmParams.baroInitCount) * 100; 
      stateEstimator_DW.stateEstInitPct = ((stateEstimator_DW.imuIdx +
        stateEstimator_DW.magIdx) + stateEstimator_DW.baroIdx) /
        ((rtu_stateEstSmParams->imuInitCount +
          rtu_stateEstSmParams->magInitCount) +
         rtu_stateEstSmParams->baroInitCount) * 100.0F;
    }

    //
    // '<S5>:1:133' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
    // '<S5>:1:134' bodyRatesOut_radps = filtBodyRatesIn_radps;
    // '<S5>:1:135' normMagVecOut_nd = normMagVecIn_nd;
    // '<S5>:1:136' isMagValid = isMagDataValid;
    *isMagValid = rtu_magData->isMagDataValid;

    // '<S5>:1:137' latLonAltOut = latLonAltIn;
    bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
    stateEstimator_DW.bodyRatesOut_radps[0] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
    stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
    stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
    bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
    stateEstimator_DW.bodyRatesOut_radps[1] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
    stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
    stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
    bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
    stateEstimator_DW.bodyRatesOut_radps[2] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
    stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
    stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

    // '<S5>:1:138' isGpsValid = isGpsDataValid;
    *isGpsValid = rtu_gpsData->isGpsDataValid;

    // '<S5>:1:139' baroAltOut_m = 0;
    *baroAltOut_m = 0.0F;

    // '<S5>:1:140' isBaroValid = isBaroDataValid;
    *isBaroValid = rtu_baroData->isBaroDataValid;

    //
    // '<S5>:1:143' initialStates = [initialQuat; 0; 0; 0; 0; 0; 0; gyroBias_radps; accelBias_mps2; nedMagVecNorm_nd; magBias_nd; baroBias_m]; 
    stateEstimator_DW.initialStates[0] = stateEstimator_DW.initialQuat[0];
    stateEstimator_DW.initialStates[1] = stateEstimator_DW.initialQuat[1];
    stateEstimator_DW.initialStates[2] = stateEstimator_DW.initialQuat[2];
    stateEstimator_DW.initialStates[3] = stateEstimator_DW.initialQuat[3];
    stateEstimator_DW.initialStates[4] = 0.0F;
    stateEstimator_DW.initialStates[5] = 0.0F;
    stateEstimator_DW.initialStates[6] = 0.0F;
    stateEstimator_DW.initialStates[7] = 0.0F;
    stateEstimator_DW.initialStates[8] = 0.0F;
    stateEstimator_DW.initialStates[9] = 0.0F;
    stateEstimator_DW.initialStates[10] = stateEstimator_DW.gyroBias_radps[0];
    stateEstimator_DW.initialStates[13] = stateEstimator_DW.accelBias_mps2[0];
    stateEstimator_DW.initialStates[16] = stateEstimator_DW.nedMagVecNorm_nd[0];
    stateEstimator_DW.initialStates[19] = stateEstimator_DW.magBias_nd[0];
    stateEstimator_DW.initialStates[11] = stateEstimator_DW.gyroBias_radps[1];
    stateEstimator_DW.initialStates[14] = stateEstimator_DW.accelBias_mps2[1];
    stateEstimator_DW.initialStates[17] = stateEstimator_DW.nedMagVecNorm_nd[1];
    stateEstimator_DW.initialStates[20] = stateEstimator_DW.magBias_nd[1];
    stateEstimator_DW.initialStates[12] = stateEstimator_DW.gyroBias_radps[2];
    stateEstimator_DW.initialStates[15] = stateEstimator_DW.accelBias_mps2[2];
    stateEstimator_DW.initialStates[18] = stateEstimator_DW.nedMagVecNorm_nd[2];
    stateEstimator_DW.initialStates[21] = stateEstimator_DW.magBias_nd[2];
    stateEstimator_DW.initialStates[22] = stateEstimator_DW.baroBias_m;

    // Compute Body To NED DCM
    // '<S5>:1:145' initialDcmBodyToNed = quatToDcm_function(initialStates(1:4)); 
    // Quaternions
    // 'quatToDcm_function:3' q0 = quat(1);
    // 'quatToDcm_function:4' q1 = quat(2);
    // 'quatToDcm_function:5' q2 = quat(3);
    // 'quatToDcm_function:6' q3 = quat(4);
    //  Direction Cosine Matrix (DCM) from body cooridinates to NED coordinates
    //  expressed using quaternions.
    // 'quatToDcm_function:10' dcmBodyToNed = [1-2*(q2^2+q3^2), 2*(q1*q2-q3*q0), 2*(q1*q3+q2*q0); 
    // 'quatToDcm_function:11'     2*(q1*q2+q3*q0), 1-2*(q1^2+q3^2), 2*(q2*q3-q1*q0); 
    // 'quatToDcm_function:12'     2*(q1*q3-q2*q0), 2*(q2*q3+q1*q0), 1-2*(q1^2+q2^2)]; 
    imuDelta_idx_0 = stateEstimator_DW.initialStates[3] *
      stateEstimator_DW.initialStates[3];
    imuDelta_idx_3 = stateEstimator_DW.initialStates[2] *
      stateEstimator_DW.initialStates[2];
    stateEstimator_DW.initialDcmBodyToNed[0] = 1.0F - (imuDelta_idx_3 +
      imuDelta_idx_0) * 2.0F;
    imuDelta_idx_1 = stateEstimator_DW.initialStates[1] *
      stateEstimator_DW.initialStates[2];
    imuDelta_idx_4 = stateEstimator_DW.initialStates[0] *
      stateEstimator_DW.initialStates[3];
    stateEstimator_DW.initialDcmBodyToNed[3] = (imuDelta_idx_1 - imuDelta_idx_4)
      * 2.0F;
    imuDelta_idx_2 = stateEstimator_DW.initialStates[1] *
      stateEstimator_DW.initialStates[3];
    imuDelta_idx_5 = stateEstimator_DW.initialStates[0] *
      stateEstimator_DW.initialStates[2];
    stateEstimator_DW.initialDcmBodyToNed[6] = (imuDelta_idx_2 + imuDelta_idx_5)
      * 2.0F;
    stateEstimator_DW.initialDcmBodyToNed[1] = (imuDelta_idx_1 + imuDelta_idx_4)
      * 2.0F;
    imuDelta_idx_1 = stateEstimator_DW.initialStates[1] *
      stateEstimator_DW.initialStates[1];
    stateEstimator_DW.initialDcmBodyToNed[4] = 1.0F - (imuDelta_idx_1 +
      imuDelta_idx_0) * 2.0F;
    imuDelta_idx_0 = stateEstimator_DW.initialStates[2] *
      stateEstimator_DW.initialStates[3];
    imuDelta_idx_4 = stateEstimator_DW.initialStates[0] *
      stateEstimator_DW.initialStates[1];
    stateEstimator_DW.initialDcmBodyToNed[7] = (imuDelta_idx_0 - imuDelta_idx_4)
      * 2.0F;
    stateEstimator_DW.initialDcmBodyToNed[2] = (imuDelta_idx_2 - imuDelta_idx_5)
      * 2.0F;
    stateEstimator_DW.initialDcmBodyToNed[5] = (imuDelta_idx_0 + imuDelta_idx_4)
      * 2.0F;
    stateEstimator_DW.initialDcmBodyToNed[8] = 1.0F - (imuDelta_idx_1 +
      imuDelta_idx_3) * 2.0F;
  }
}

// Function for Chart: '<Root>/estimatorStateMachine'
void stateEstimator::state_enter_atomic_RUN_INIT_GPS(enumStateEstimateMode *mode,
  boolean_T *isMagValid, boolean_T *isGpsValid, real32_T *baroAltOut_m,
  boolean_T *isBaroValid, real32_T bodyAccelsOut_mps2[3], const busMagData
  *rtu_magData, const busGpsData *rtu_gpsData, const busBaroData *rtu_baroData)
{
  // Entry 'RUN_INIT_GPS': '<S5>:64'
  // EKF RUNNING WITHOUT GPS BUT WE HAVE HEALTHY GPS SIGNAL
  // START INITIALIZING GPS
  // '<S5>:64:5' gpsValidCount = 0;
  stateEstimator_DW.gpsValidCount = 0U;

  //  Index to keep track of how many gps readings we have summed
  //  so far
  // '<S5>:64:8' gpsIdx = 0;
  stateEstimator_DW.gpsIdx = 0.0;

  // '<S5>:64:9' refLatLonAlt = [0; 0; 0];
  // '<S5>:64:10' isPosInitialized = false;
  stateEstimator_DW.isPosInitialized = false;

  // '<S5>:64:11' mode = enumStateEstimateMode.RUN_INIT_GPS;
  *mode = enumStateEstimateMode::RUN_INIT_GPS;

  // Chart: '<Root>/estimatorStateMachine'
  // '<S5>:64:12' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
  // '<S5>:64:13' bodyRatesOut_radps = filtBodyRatesIn_radps;
  // '<S5>:64:14' normMagVecOut_nd = normMagVecIn_nd;
  // '<S5>:64:15' isMagValid = isMagDataValid;
  *isMagValid = rtu_magData->isMagDataValid;

  // '<S5>:64:16' latLonAltOut = latLonAltIn;
  stateEstimator_DW.refLatLonAlt[0] = 0.0;

  // Product: '<S4>/Product'
  bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];

  // SignalConversion generated from: '<S5>/ SFunction '
  stateEstimator_DW.bodyRatesOut_radps[0] =
    stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];

  // Product: '<S9>/Divide'
  stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];

  // Chart: '<Root>/estimatorStateMachine'
  stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
  stateEstimator_DW.refLatLonAlt[1] = 0.0;

  // Product: '<S4>/Product'
  bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];

  // SignalConversion generated from: '<S5>/ SFunction '
  stateEstimator_DW.bodyRatesOut_radps[1] =
    stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];

  // Product: '<S9>/Divide'
  stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];

  // Chart: '<Root>/estimatorStateMachine'
  stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
  stateEstimator_DW.refLatLonAlt[2] = 0.0;

  // Product: '<S4>/Product'
  bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];

  // SignalConversion generated from: '<S5>/ SFunction '
  stateEstimator_DW.bodyRatesOut_radps[2] =
    stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];

  // Product: '<S9>/Divide'
  stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];

  // Chart: '<Root>/estimatorStateMachine'
  stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

  // '<S5>:64:17' isGpsValid = isGpsDataValid;
  *isGpsValid = rtu_gpsData->isGpsDataValid;

  // '<S5>:64:18' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
  *baroAltOut_m = stateEstimator_DW.Divide1 - stateEstimator_DW.baroInitAltMean;

  // Chart: '<Root>/estimatorStateMachine'
  // '<S5>:64:19' isBaroValid = isBaroDataValid;
  *isBaroValid = rtu_baroData->isBaroDataValid;

  //
}

// System initialize for referenced model: 'stateEstimator'
void stateEstimator::init(busStateEstimatorDebug *rty_stateEstimatorDebug)
{
  static const int8_T tmp[23]{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0 };

  // InitializeConditions for Delay: '<S1>/Delay'
  stateEstimator_DW.icLoad = true;

  // InitializeConditions for Delay: '<S1>/Delay2'
  stateEstimator_DW.icLoad_j = true;
  for (int32_T i{0}; i < 23; i++) {
    // InitializeConditions for UnitDelay: '<Root>/Unit Delay'
    stateEstimator_DW.UnitDelay_DSTATE[i] = rtCP_UnitDelay_InitialCondition[i];

    // SystemInitialize for Chart: '<Root>/estimatorStateMachine' incorporates:
    //   UnitDelay: '<Root>/Unit Delay'

    stateEstimator_DW.initialStates[i] = tmp[i];
  }

  // SystemInitialize for Chart: '<Root>/estimatorStateMachine'
  stateEstimator_DW.initialDcmBodyToNed[0] = 1.0F;
  stateEstimator_DW.initialDcmBodyToNed[1] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[2] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[3] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[4] = 1.0F;
  stateEstimator_DW.initialDcmBodyToNed[5] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[6] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[7] = 0.0F;
  stateEstimator_DW.initialDcmBodyToNed[8] = 1.0F;

  // SystemInitialize for BusCreator: '<Root>/Bus Creator' incorporates:
  //   Chart: '<Root>/estimatorStateMachine'

  rty_stateEstimatorDebug->stateEstInitPct = 0.0F;
  rty_stateEstimatorDebug->smMode = enumStateEstimateMode::NONE;
}

// Output and update for referenced model: 'stateEstimator'
void stateEstimator::step(const busImuData *rtu_imuData, const busMagData
  *rtu_magData, const busGpsData *rtu_gpsData, const busBaroData *rtu_baroData,
  const busLidarData *rtu_lidarData, const busImuNtchFiltParams
  *rtu_imuNotchFiltParams, const busAccelParams *rtu_accelParams, const
  busMagParams *rtu_magParams, const busLidarParams *rtu_lidarParams, const
  busStateEstSmParams *rtu_stateEstSmParams, const real32_T rtu_processNoiseQ
  [529], const real32_T rtu_measNoiseR[121], const real32_T rtu_initCovP[529],
  const real32_T rtu_processNoiseNoGpsQ[324], const real32_T
  rtu_measNoiseNoGpsR[64], const real32_T rtu_initCovNoGpsP[324], const real32_T
  *rtu_gEarth_mps2, real32_T rty_states[23], real32_T rty_eulAng_rad[3],
  real32_T rty_dcmNedToBody[9], real32_T rty_dcmNedToFep[9], real32_T
  rty_bodyAccels_mps2[3], busStateEstimatorDebug *rty_stateEstimatorDebug)
{
  static const int8_T d_measJac[138]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const int8_T e[138]{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const int8_T b_measJac[23]{ 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T c_measJac[23]{ 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const int8_T b_measJac_0[18]{ 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1 };

  static const int8_T c[18]{ 0, 1, 2, 3, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22 };

  static const int8_T c_measJac_0[18]{ 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0 };

  real_T rtb_nedPos_m_idx_0;
  real_T rtb_nedPos_m_idx_1;
  real_T rtb_nedPos_m_idx_2;
  int32_T b_tmp1_tmp;
  int32_T covP_tmp;
  int32_T d_tmp1_tmp;
  int32_T f_y_tmp_0;
  int32_T i;
  int32_T i_0;
  int32_T i_1;
  real32_T K_1[529];
  real32_T covP[529];
  real32_T K_2[324];
  real32_T covP_0[324];
  real32_T b_K[138];
  real32_T d_tmp1[138];
  real32_T K[69];
  real32_T b_tmp1[69];
  real32_T measJac[69];
  real32_T tmp_1[58];
  real32_T K_0[54];
  real32_T c_tmp1_0[54];
  real32_T measJac_0[54];
  real32_T b_K_tmp[36];
  real32_T stateJac[28];
  real32_T c_tmp1[23];
  real32_T rtb_UnitDelay[23];
  real32_T d_tmp1_0[18];
  real32_T rtb_states_a[18];
  real32_T tmp[12];
  real32_T measJac_1[9];
  real32_T tmp_0[4];
  real32_T bodyAccelsOut_mps2[3];
  real32_T baroAltOut_m;
  real32_T rtb_Product1;
  real32_T rtb_Product2;
  real32_T rtb_UnitDelay_e;
  real32_T rtb_XAxis;
  real32_T rtb_XAxis1;
  real32_T rtb_XAxis2;
  real32_T rtb_dcmBodyToNed_k_idx_0;
  real32_T rtb_dcmBodyToNed_k_idx_1;
  real32_T rtb_dcmBodyToNed_k_idx_2;
  real32_T rtb_dcmBodyToNed_k_idx_3;
  real32_T rtb_dcmBodyToNed_k_idx_4;
  real32_T rtb_dcmBodyToNed_k_idx_5;
  real32_T rtb_dcmBodyToNed_k_idx_6;
  real32_T rtb_dcmBodyToNed_k_idx_7;
  real32_T rtb_dcmBodyToNed_k_idx_8;
  real32_T tmp1;
  real32_T tmp11;
  real32_T tmp17;
  real32_T tmp2;
  real32_T tmp6;
  real32_T tmp8;
  real32_T tmp9;
  int8_T f_y_tmp[23];
  int8_T g_y_tmp[18];
  const int8_T *b_K_tmp_0;
  boolean_T isBaroValid;
  boolean_T isGpsValid;
  boolean_T isMagValid;
  boolean_T rtb_AND1;
  boolean_T rtb_Compare;
  enumStateEstimateMode mode;

  // UnitDelay: '<Root>/Unit Delay'
  std::memcpy(&rtb_UnitDelay[0], &stateEstimator_DW.UnitDelay_DSTATE[0], 23U *
              sizeof(real32_T));

  // DiscreteTransferFcn: '<S23>/X Axis'
  stateEstimator_DW.XAxis_tmp = (rtu_imuData->bodyAccels_mps2[0] -
    stateEstimator_DW.XAxis_states[0] *
    rtu_imuNotchFiltParams->accelNtchFilt.xDen[1]) -
    stateEstimator_DW.XAxis_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.xDen[2];
  rtb_Product2 = (rtu_imuNotchFiltParams->accelNtchFilt.xNum[0] *
                  stateEstimator_DW.XAxis_tmp + stateEstimator_DW.XAxis_states[0]
                  * rtu_imuNotchFiltParams->accelNtchFilt.xNum[1]) +
    stateEstimator_DW.XAxis_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.xNum[2];

  // DiscreteTransferFcn: '<S23>/X Axis1'
  stateEstimator_DW.XAxis1_tmp = (rtu_imuData->bodyAccels_mps2[1] -
    stateEstimator_DW.XAxis1_states[0] *
    rtu_imuNotchFiltParams->accelNtchFilt.yDen[1]) -
    stateEstimator_DW.XAxis1_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.yDen[2];
  rtb_Product1 = (rtu_imuNotchFiltParams->accelNtchFilt.yNum[0] *
                  stateEstimator_DW.XAxis1_tmp +
                  stateEstimator_DW.XAxis1_states[0] *
                  rtu_imuNotchFiltParams->accelNtchFilt.yNum[1]) +
    stateEstimator_DW.XAxis1_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.yNum[2];

  // DiscreteTransferFcn: '<S23>/X Axis2'
  stateEstimator_DW.XAxis2_tmp = (rtu_imuData->bodyAccels_mps2[2] -
    stateEstimator_DW.XAxis2_states[0] *
    rtu_imuNotchFiltParams->accelNtchFilt.zDen[1]) -
    stateEstimator_DW.XAxis2_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.zDen[2];
  rtb_UnitDelay_e = (rtu_imuNotchFiltParams->accelNtchFilt.zNum[0] *
                     stateEstimator_DW.XAxis2_tmp +
                     stateEstimator_DW.XAxis2_states[0] *
                     rtu_imuNotchFiltParams->accelNtchFilt.zNum[1]) +
    stateEstimator_DW.XAxis2_states[1] *
    rtu_imuNotchFiltParams->accelNtchFilt.zNum[2];

  // Product: '<S4>/Product' incorporates:
  //   Product: '<S4>/Matrix Multiply'
  //   SignalConversion generated from: '<S4>/Matrix Multiply'
  //   Sum: '<S4>/Sum'

  stateEstimator_DW.Product[0] = (((rtu_accelParams->scaleAlignMat_nd[0] *
    rtb_Product2 + rtu_accelParams->scaleAlignMat_nd[3] * rtb_Product1) +
    rtu_accelParams->scaleAlignMat_nd[6] * rtb_UnitDelay_e) +
    rtu_accelParams->offset_nd[0]) * *rtu_gEarth_mps2;
  stateEstimator_DW.Product[1] = (((rtu_accelParams->scaleAlignMat_nd[1] *
    rtb_Product2 + rtu_accelParams->scaleAlignMat_nd[4] * rtb_Product1) +
    rtu_accelParams->scaleAlignMat_nd[7] * rtb_UnitDelay_e) +
    rtu_accelParams->offset_nd[1]) * *rtu_gEarth_mps2;
  stateEstimator_DW.Product[2] = (((rtu_accelParams->scaleAlignMat_nd[2] *
    rtb_Product2 + rtu_accelParams->scaleAlignMat_nd[5] * rtb_Product1) +
    rtu_accelParams->scaleAlignMat_nd[8] * rtb_UnitDelay_e) +
    rtu_accelParams->offset_nd[2]) * *rtu_gEarth_mps2;

  // DiscreteTransferFcn: '<S24>/X Axis'
  stateEstimator_DW.XAxis_tmp_k = (rtu_imuData->bodyRates_radps[0] -
    stateEstimator_DW.XAxis_states_b[0] *
    rtu_imuNotchFiltParams->gyroNtchFilt.xDen[1]) -
    stateEstimator_DW.XAxis_states_b[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.xDen[2];

  // DiscreteTransferFcn: '<S24>/X Axis1'
  stateEstimator_DW.XAxis1_tmp_m = (rtu_imuData->bodyRates_radps[1] -
    stateEstimator_DW.XAxis1_states_m[0] *
    rtu_imuNotchFiltParams->gyroNtchFilt.yDen[1]) -
    stateEstimator_DW.XAxis1_states_m[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.yDen[2];

  // DiscreteTransferFcn: '<S24>/X Axis2'
  stateEstimator_DW.XAxis2_tmp_g = (rtu_imuData->bodyRates_radps[2] -
    stateEstimator_DW.XAxis2_states_h[0] *
    rtu_imuNotchFiltParams->gyroNtchFilt.zDen[1]) -
    stateEstimator_DW.XAxis2_states_h[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.zDen[2];

  // Sum: '<S9>/Sum'
  rtb_Product2 = rtu_magData->bodyMagVector_uT[0] - rtu_magParams->offset_uT[0];
  rtb_dcmBodyToNed_k_idx_3 = rtu_magData->bodyMagVector_uT[1] -
    rtu_magParams->offset_uT[1];
  rtb_dcmBodyToNed_k_idx_6 = rtu_magData->bodyMagVector_uT[2] -
    rtu_magParams->offset_uT[2];

  // Product: '<S9>/Matrix Multiply'
  baroAltOut_m = (rtu_magParams->scaleAlignMat_nd[0] * rtb_Product2 +
                  rtu_magParams->scaleAlignMat_nd[3] * rtb_dcmBodyToNed_k_idx_3)
    + rtu_magParams->scaleAlignMat_nd[6] * rtb_dcmBodyToNed_k_idx_6;
  rtb_dcmBodyToNed_k_idx_0 = (rtu_magParams->scaleAlignMat_nd[1] * rtb_Product2
    + rtu_magParams->scaleAlignMat_nd[4] * rtb_dcmBodyToNed_k_idx_3) +
    rtu_magParams->scaleAlignMat_nd[7] * rtb_dcmBodyToNed_k_idx_6;
  rtb_dcmBodyToNed_k_idx_3 = (rtu_magParams->scaleAlignMat_nd[2] * rtb_Product2
    + rtu_magParams->scaleAlignMat_nd[5] * rtb_dcmBodyToNed_k_idx_3) +
    rtu_magParams->scaleAlignMat_nd[8] * rtb_dcmBodyToNed_k_idx_6;

  // MinMax: '<S9>/Max' incorporates:
  //   Constant: '<S9>/Constant2'
  //   Product: '<S9>/Matrix Multiply1'
  //   Sqrt: '<S9>/Sqrt'

  rtb_Product2 = std::fmax(std::sqrt((baroAltOut_m * baroAltOut_m +
    rtb_dcmBodyToNed_k_idx_0 * rtb_dcmBodyToNed_k_idx_0) +
    rtb_dcmBodyToNed_k_idx_3 * rtb_dcmBodyToNed_k_idx_3), 1.0E-7F);

  // Product: '<S9>/Divide'
  stateEstimator_DW.Divide[0] = baroAltOut_m / rtb_Product2;
  stateEstimator_DW.Divide[1] = rtb_dcmBodyToNed_k_idx_0 / rtb_Product2;
  stateEstimator_DW.Divide[2] = rtb_dcmBodyToNed_k_idx_3 / rtb_Product2;

  // Product: '<S10>/Divide1' incorporates:
  //   Constant: '<S10>/Constant'
  //   Constant: '<S10>/Constant2'
  //   Constant: '<S10>/Constant3'
  //   Math: '<S10>/Power'
  //   Product: '<S10>/Divide'
  //   Sum: '<S10>/Sum'

  stateEstimator_DW.Divide1 = (1.0F - std::pow(rtu_baroData->pressure_pa /
    101325.0F, 0.190294951F)) * 44330.0F;

  // SignalConversion generated from: '<S5>/ SFunction ' incorporates:
  //   Chart: '<Root>/estimatorStateMachine'
  //   DiscreteTransferFcn: '<S24>/X Axis'
  //   DiscreteTransferFcn: '<S24>/X Axis1'
  //   DiscreteTransferFcn: '<S24>/X Axis2'

  stateEstimator_DW.TmpSignalConversionAtSFunctionI[0] =
    (rtu_imuNotchFiltParams->gyroNtchFilt.xNum[0] *
     stateEstimator_DW.XAxis_tmp_k + stateEstimator_DW.XAxis_states_b[0] *
     rtu_imuNotchFiltParams->gyroNtchFilt.xNum[1]) +
    stateEstimator_DW.XAxis_states_b[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.xNum[2];
  stateEstimator_DW.TmpSignalConversionAtSFunctionI[1] =
    (rtu_imuNotchFiltParams->gyroNtchFilt.yNum[0] *
     stateEstimator_DW.XAxis1_tmp_m + stateEstimator_DW.XAxis1_states_m[0] *
     rtu_imuNotchFiltParams->gyroNtchFilt.yNum[1]) +
    stateEstimator_DW.XAxis1_states_m[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.yNum[2];
  stateEstimator_DW.TmpSignalConversionAtSFunctionI[2] =
    (rtu_imuNotchFiltParams->gyroNtchFilt.zNum[0] *
     stateEstimator_DW.XAxis2_tmp_g + stateEstimator_DW.XAxis2_states_h[0] *
     rtu_imuNotchFiltParams->gyroNtchFilt.zNum[1]) +
    stateEstimator_DW.XAxis2_states_h[1] *
    rtu_imuNotchFiltParams->gyroNtchFilt.zNum[2];

  // Chart: '<Root>/estimatorStateMachine' incorporates:
  //   Product: '<S4>/Product'
  //   Product: '<S9>/Divide'
  //   SignalConversion generated from: '<S5>/ SFunction '

  // Gateway: estimatorStateMachine
  // During: estimatorStateMachine
  if (stateEstimator_DW.is_active_c3_stateEstimator == 0U) {
    // Entry: estimatorStateMachine
    stateEstimator_DW.is_active_c3_stateEstimator = 1U;

    // Entry Internal: estimatorStateMachine
    // Transition: '<S5>:2'
    stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_INITIALIZE;

    // Entry 'INITIALIZE': '<S5>:1'
    // Variables to set on entry
    // '<S5>:1:4' mode = enumStateEstimateMode.INITIALIZE;
    mode = enumStateEstimateMode::INITIALIZE;

    // '<S5>:1:5' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
    // '<S5>:1:6' bodyRatesOut_radps = filtBodyRatesIn_radps;
    // '<S5>:1:7' normMagVecOut_nd = normMagVecIn_nd;
    // '<S5>:1:8' isMagValid = isMagDataValid;
    isMagValid = rtu_magData->isMagDataValid;

    // '<S5>:1:9' latLonAltOut = latLonAltIn;
    bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
    stateEstimator_DW.bodyRatesOut_radps[0] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
    stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
    stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
    bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
    stateEstimator_DW.bodyRatesOut_radps[1] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
    stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
    stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
    bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
    stateEstimator_DW.bodyRatesOut_radps[2] =
      stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
    stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
    stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

    // '<S5>:1:10' isGpsValid = isGpsDataValid;
    isGpsValid = rtu_gpsData->isGpsDataValid;

    // '<S5>:1:11' baroAltOut_m = 0;
    baroAltOut_m = 0.0F;

    // '<S5>:1:12' isBaroValid = isBaroDataValid;
    isBaroValid = rtu_baroData->isBaroDataValid;

    // '<S5>:1:13' resetStates = true;
    stateEstimator_DW.resetStates = true;

    //  Index to keep track of how many imu readings we have summed
    //  so far
    // '<S5>:1:16' imuIdx = 0;
    stateEstimator_DW.imuIdx = 0.0F;

    // '<S5>:1:17' imuMean = [0; 0; 0; 0; 0; 0];
    stateEstimator_DW.imuMean[0] = 0.0F;
    stateEstimator_DW.imuMean[1] = 0.0F;
    stateEstimator_DW.imuMean[2] = 0.0F;
    stateEstimator_DW.imuMean[3] = 0.0F;
    stateEstimator_DW.imuMean[4] = 0.0F;
    stateEstimator_DW.imuMean[5] = 0.0F;

    // '<S5>:1:18' accelBias_mps2 = [0; 0; 0];
    // '<S5>:1:19' gyroBias_radps = [0; 0; 0];
    stateEstimator_DW.accelBias_mps2[0] = 0.0F;
    stateEstimator_DW.gyroBias_radps[0] = 0.0F;
    stateEstimator_DW.accelBias_mps2[1] = 0.0F;
    stateEstimator_DW.gyroBias_radps[1] = 0.0F;
    stateEstimator_DW.accelBias_mps2[2] = 0.0F;
    stateEstimator_DW.gyroBias_radps[2] = 0.0F;

    // '<S5>:1:20' imuM2 = [0; 0; 0; 0; 0; 0];
    stateEstimator_DW.imuM2[0] = 0.0F;
    stateEstimator_DW.imuM2[1] = 0.0F;
    stateEstimator_DW.imuM2[2] = 0.0F;
    stateEstimator_DW.imuM2[3] = 0.0F;
    stateEstimator_DW.imuM2[4] = 0.0F;
    stateEstimator_DW.imuM2[5] = 0.0F;

    // '<S5>:1:21' initialQuat = [1; 0; 0; 0];
    stateEstimator_DW.initialQuat[0] = 1.0F;
    stateEstimator_DW.initialQuat[1] = 0.0F;
    stateEstimator_DW.initialQuat[2] = 0.0F;
    stateEstimator_DW.initialQuat[3] = 0.0F;

    // '<S5>:1:22' isAttInitialized = false;
    stateEstimator_DW.isAttInitialized = false;

    //  Index to keep track of how many mag readings we have summed
    //  so far
    // '<S5>:1:25' magIdx = 0;
    stateEstimator_DW.magIdx = 0.0F;

    // '<S5>:1:26' magMean = [0; 0; 0];
    // '<S5>:1:27' magM2 = [0; 0; 0];
    // '<S5>:1:28' magBias_nd = [0; 0; 0];
    // '<S5>:1:29' nedMagVecNorm_nd = [0; 0; 0];
    //  Index to keep track of how many gps readings we have summed
    //  so far
    // '<S5>:1:32' gpsIdx = 0;
    stateEstimator_DW.gpsIdx = 0.0;

    // '<S5>:1:33' refLatLonAlt = [0; 0; 0];
    stateEstimator_DW.magMean[0] = 0.0F;
    stateEstimator_DW.magM2[0] = 0.0F;
    stateEstimator_DW.magBias_nd[0] = 0.0F;
    stateEstimator_DW.nedMagVecNorm_nd[0] = 0.0F;
    stateEstimator_DW.refLatLonAlt[0] = 0.0;
    stateEstimator_DW.magMean[1] = 0.0F;
    stateEstimator_DW.magM2[1] = 0.0F;
    stateEstimator_DW.magBias_nd[1] = 0.0F;
    stateEstimator_DW.nedMagVecNorm_nd[1] = 0.0F;
    stateEstimator_DW.refLatLonAlt[1] = 0.0;
    stateEstimator_DW.magMean[2] = 0.0F;
    stateEstimator_DW.magM2[2] = 0.0F;
    stateEstimator_DW.magBias_nd[2] = 0.0F;
    stateEstimator_DW.nedMagVecNorm_nd[2] = 0.0F;
    stateEstimator_DW.refLatLonAlt[2] = 0.0;

    // '<S5>:1:34' isPosInitialized = false;
    stateEstimator_DW.isPosInitialized = false;

    //  Index to keep track of how many baro readings we have summed
    //  so far
    // '<S5>:1:37' baroIdx = 0;
    stateEstimator_DW.baroIdx = 0.0F;

    // '<S5>:1:38' baroInitAltMean = 0;
    stateEstimator_DW.baroInitAltMean = 0.0F;

    // '<S5>:1:39' baroInitAltM2 = 0;
    stateEstimator_DW.baroInitAltM2 = 0.0F;

    // '<S5>:1:40' baroBias_m = 0;
    stateEstimator_DW.baroBias_m = 0.0F;

    // '<S5>:1:41' isBaroInitialized = false;
    stateEstimator_DW.isBaroInitialized = false;

    // End of entry stage
  } else {
    switch (stateEstimator_DW.is_c3_stateEstimator) {
     case stateEstimator_IN_INITIALIZE:
      stateEstimator_INITIALIZE(&mode, &isMagValid, &isGpsValid, &baroAltOut_m,
        &isBaroValid, bodyAccelsOut_mps2, rtu_magData, rtu_gpsData, rtu_baroData,
        rtu_stateEstSmParams);
      break;

     case stateEstimator_IN_RUN:
      stateEstimator_DW.resetStates = false;
      mode = enumStateEstimateMode::RUN;

      // During 'RUN': '<S5>:43'
      // '<S5>:68:1' sf_internal_predicateOutput = duration(~isGpsDataValid) >=  ... 
      // '<S5>:68:2' stateEstSmParams.gpsLossCheckDuration_s;
      if (rtu_gpsData->isGpsDataValid) {
        stateEstimator_DW.durationCounter_1 = 0;
      }

      if (static_cast<real_T>(stateEstimator_DW.durationCounter_1) >=
          rtu_stateEstSmParams->gpsLossCheckDuration_s * 250.0F) {
        // Transition: '<S5>:68'
        stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_RUN_GPS_LOST;

        // Entry 'RUN_GPS_LOST': '<S5>:67'
        // GPS WAS LOST
        // '<S5>:67:4' gpsValidCount = 0;
        stateEstimator_DW.gpsValidCount = 0U;

        // '<S5>:67:5' mode = enumStateEstimateMode.RUN_GPS_LOST;
        mode = enumStateEstimateMode::RUN_GPS_LOST;

        // '<S5>:67:6' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
        // '<S5>:67:7' bodyRatesOut_radps = filtBodyRatesIn_radps;
        // '<S5>:67:8' normMagVecOut_nd = normMagVecIn_nd;
        // '<S5>:67:9' isMagValid = isMagDataValid;
        isMagValid = rtu_magData->isMagDataValid;

        // '<S5>:67:10' latLonAltOut = latLonAltIn;
        bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
        stateEstimator_DW.bodyRatesOut_radps[0] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
        stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
        stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
        bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
        stateEstimator_DW.bodyRatesOut_radps[1] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
        stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
        stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
        bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
        stateEstimator_DW.bodyRatesOut_radps[2] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
        stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
        stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

        // '<S5>:67:11' isGpsValid = isGpsDataValid;
        isGpsValid = rtu_gpsData->isGpsDataValid;

        // '<S5>:67:12' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
        baroAltOut_m = stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean;

        // '<S5>:67:13' isBaroValid = isBaroDataValid;
        isBaroValid = rtu_baroData->isBaroDataValid;

        //
      } else {
        // FULL EKF WITH GPS IS RUNNING
        // '<S5>:43:16' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
        // '<S5>:43:17' bodyRatesOut_radps = filtBodyRatesIn_radps;
        // '<S5>:43:18' normMagVecOut_nd = normMagVecIn_nd;
        // '<S5>:43:19' isMagValid = isMagDataValid;
        isMagValid = rtu_magData->isMagDataValid;

        // '<S5>:43:20' latLonAltOut = latLonAltIn;
        bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
        stateEstimator_DW.bodyRatesOut_radps[0] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
        stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
        stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
        bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
        stateEstimator_DW.bodyRatesOut_radps[1] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
        stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
        stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
        bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
        stateEstimator_DW.bodyRatesOut_radps[2] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
        stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
        stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

        // '<S5>:43:21' isGpsValid = isGpsDataValid;
        isGpsValid = rtu_gpsData->isGpsDataValid;

        // '<S5>:43:22' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
        baroAltOut_m = stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean;

        // '<S5>:43:23' isBaroValid = isBaroDataValid;
        isBaroValid = rtu_baroData->isBaroDataValid;
      }
      break;

     case stateEstimator_IN_RUN_GPS_LOST:
      mode = enumStateEstimateMode::RUN_GPS_LOST;

      // During 'RUN_GPS_LOST': '<S5>:67'
      // '<S5>:71:1' sf_internal_predicateOutput = gpsValidCount >= stateEstSmParams.desValidGpsCount; 
      if (stateEstimator_DW.gpsValidCount >=
          rtu_stateEstSmParams->desValidGpsCount) {
        // Transition: '<S5>:71'
        stateEstimator_DW.durationCounter_1_c = 0;
        stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_RUN_INIT_GPS;
        state_enter_atomic_RUN_INIT_GPS(&mode, &isMagValid, &isGpsValid,
          &baroAltOut_m, &isBaroValid, bodyAccelsOut_mps2, rtu_magData,
          rtu_gpsData, rtu_baroData);
      } else {
        // GPS WAS LOST
        // '<S5>:67:16' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
        // '<S5>:67:17' bodyRatesOut_radps = filtBodyRatesIn_radps;
        // '<S5>:67:18' normMagVecOut_nd = normMagVecIn_nd;
        // '<S5>:67:19' isMagValid = isMagDataValid;
        isMagValid = rtu_magData->isMagDataValid;

        // '<S5>:67:20' latLonAltOut = latLonAltIn;
        bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
        stateEstimator_DW.bodyRatesOut_radps[0] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
        stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
        stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
        bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
        stateEstimator_DW.bodyRatesOut_radps[1] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
        stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
        stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
        bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
        stateEstimator_DW.bodyRatesOut_radps[2] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
        stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
        stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

        // '<S5>:67:21' isGpsValid = isGpsDataValid;
        isGpsValid = rtu_gpsData->isGpsDataValid;

        // '<S5>:67:22' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
        baroAltOut_m = stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean;

        // '<S5>:67:23' isBaroValid = isBaroDataValid;
        isBaroValid = rtu_baroData->isBaroDataValid;

        //
        // '<S5>:67:25' if(isGpsValid)
        if (rtu_gpsData->isGpsDataValid) {
          // '<S5>:67:26' gpsValidCount = gpsValidCount + 1;
          stateEstimator_DW.gpsValidCount = static_cast<uint16_T>
            (stateEstimator_DW.gpsValidCount + 1U);
        }
      }
      break;

     case stateEstima_IN_RUN_GPS_NOT_INIT:
      stateEstimator_DW.resetStates = false;
      mode = enumStateEstimateMode::RUN_GPS_NOT_INIT;

      // During 'RUN_GPS_NOT_INIT': '<S5>:62'
      // '<S5>:65:1' sf_internal_predicateOutput = gpsValidCount >= stateEstSmParams.desValidGpsCount; 
      if (stateEstimator_DW.gpsValidCount >=
          rtu_stateEstSmParams->desValidGpsCount) {
        // Transition: '<S5>:65'
        stateEstimator_DW.durationCounter_1_c = 0;
        stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_RUN_INIT_GPS;
        state_enter_atomic_RUN_INIT_GPS(&mode, &isMagValid, &isGpsValid,
          &baroAltOut_m, &isBaroValid, bodyAccelsOut_mps2, rtu_magData,
          rtu_gpsData, rtu_baroData);
      } else {
        // '<S5>:73:1' sf_internal_predicateOutput = duration(~isGpsDataValid) >=  ... 
        // '<S5>:73:2' stateEstSmParams.gpsLossCheckDuration_s;
        if (rtu_gpsData->isGpsDataValid) {
          stateEstimator_DW.durationCounter_1_p = 0;
        }

        if (static_cast<real_T>(stateEstimator_DW.durationCounter_1_p) >=
            rtu_stateEstSmParams->gpsLossCheckDuration_s * 250.0F) {
          // Transition: '<S5>:73'
          stateEstimator_DW.is_c3_stateEstimator =
            stateEstimator_IN_RUN_GPS_LOST;

          // Entry 'RUN_GPS_LOST': '<S5>:67'
          // GPS WAS LOST
          // '<S5>:67:4' gpsValidCount = 0;
          stateEstimator_DW.gpsValidCount = 0U;

          // '<S5>:67:5' mode = enumStateEstimateMode.RUN_GPS_LOST;
          mode = enumStateEstimateMode::RUN_GPS_LOST;

          // '<S5>:67:6' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
          // '<S5>:67:7' bodyRatesOut_radps = filtBodyRatesIn_radps;
          // '<S5>:67:8' normMagVecOut_nd = normMagVecIn_nd;
          // '<S5>:67:9' isMagValid = isMagDataValid;
          isMagValid = rtu_magData->isMagDataValid;

          // '<S5>:67:10' latLonAltOut = latLonAltIn;
          bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
          stateEstimator_DW.bodyRatesOut_radps[0] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
          stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
          stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
          bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
          stateEstimator_DW.bodyRatesOut_radps[1] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
          stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
          stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
          bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
          stateEstimator_DW.bodyRatesOut_radps[2] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
          stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
          stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

          // '<S5>:67:11' isGpsValid = isGpsDataValid;
          isGpsValid = rtu_gpsData->isGpsDataValid;

          // '<S5>:67:12' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
          baroAltOut_m = stateEstimator_DW.Divide1 -
            stateEstimator_DW.baroInitAltMean;

          // '<S5>:67:13' isBaroValid = isBaroDataValid;
          isBaroValid = rtu_baroData->isBaroDataValid;

          //
        } else {
          // EKF STARTED RUNNING WITHOUT GPS
          // '<S5>:62:18' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
          // '<S5>:62:19' bodyRatesOut_radps = filtBodyRatesIn_radps;
          // '<S5>:62:20' normMagVecOut_nd = normMagVecIn_nd;
          // '<S5>:62:21' isMagValid = isMagDataValid;
          isMagValid = rtu_magData->isMagDataValid;

          // '<S5>:62:22' latLonAltOut = latLonAltIn;
          bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
          stateEstimator_DW.bodyRatesOut_radps[0] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
          stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
          stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
          bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
          stateEstimator_DW.bodyRatesOut_radps[1] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
          stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
          stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
          bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
          stateEstimator_DW.bodyRatesOut_radps[2] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
          stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
          stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

          // '<S5>:62:23' isGpsValid = isGpsDataValid;
          isGpsValid = rtu_gpsData->isGpsDataValid;

          // '<S5>:62:24' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
          baroAltOut_m = stateEstimator_DW.Divide1 -
            stateEstimator_DW.baroInitAltMean;

          // '<S5>:62:25' isBaroValid = isBaroDataValid;
          isBaroValid = rtu_baroData->isBaroDataValid;

          //
          // '<S5>:62:27' if(isGpsValid)
          if (rtu_gpsData->isGpsDataValid) {
            // '<S5>:62:28' gpsValidCount = gpsValidCount + 1;
            stateEstimator_DW.gpsValidCount = static_cast<uint16_T>
              (stateEstimator_DW.gpsValidCount + 1U);
          }
        }
      }
      break;

     default:
      mode = enumStateEstimateMode::RUN_INIT_GPS;

      // During 'RUN_INIT_GPS': '<S5>:64'
      // '<S5>:66:1' sf_internal_predicateOutput = isPosInitialized;
      if (stateEstimator_DW.isPosInitialized) {
        // Transition: '<S5>:66'
        stateEstimator_DW.durationCounter_1 = 0;
        stateEstimator_DW.is_c3_stateEstimator = stateEstimator_IN_RUN;

        // Entry 'RUN': '<S5>:43'
        // FULL EKF WITH GPS IS RUNNING
        // '<S5>:43:4' resetStates = false;
        stateEstimator_DW.resetStates = false;

        // '<S5>:43:5' mode = enumStateEstimateMode.RUN;
        mode = enumStateEstimateMode::RUN;

        // '<S5>:43:6' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
        // '<S5>:43:7' bodyRatesOut_radps = filtBodyRatesIn_radps;
        // '<S5>:43:8' normMagVecOut_nd = normMagVecIn_nd;
        // '<S5>:43:9' isMagValid = isMagDataValid;
        isMagValid = rtu_magData->isMagDataValid;

        // '<S5>:43:10' latLonAltOut = latLonAltIn;
        bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
        stateEstimator_DW.bodyRatesOut_radps[0] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
        stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
        stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
        bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
        stateEstimator_DW.bodyRatesOut_radps[1] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
        stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
        stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
        bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
        stateEstimator_DW.bodyRatesOut_radps[2] =
          stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
        stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
        stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

        // '<S5>:43:11' isGpsValid = isGpsDataValid;
        isGpsValid = rtu_gpsData->isGpsDataValid;

        // '<S5>:43:12' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
        baroAltOut_m = stateEstimator_DW.Divide1 -
          stateEstimator_DW.baroInitAltMean;

        // '<S5>:43:13' isBaroValid = isBaroDataValid;
        isBaroValid = rtu_baroData->isBaroDataValid;

        //
      } else {
        // '<S5>:70:1' sf_internal_predicateOutput = duration(~isGpsDataValid) >=  ... 
        // '<S5>:70:2' stateEstSmParams.gpsLossCheckDuration_s;
        if (rtu_gpsData->isGpsDataValid) {
          stateEstimator_DW.durationCounter_1_c = 0;
        }

        if (static_cast<real_T>(stateEstimator_DW.durationCounter_1_c) >=
            rtu_stateEstSmParams->gpsLossCheckDuration_s * 250.0F) {
          // Transition: '<S5>:70'
          stateEstimator_DW.is_c3_stateEstimator =
            stateEstimator_IN_RUN_GPS_LOST;

          // Entry 'RUN_GPS_LOST': '<S5>:67'
          // GPS WAS LOST
          // '<S5>:67:4' gpsValidCount = 0;
          stateEstimator_DW.gpsValidCount = 0U;

          // '<S5>:67:5' mode = enumStateEstimateMode.RUN_GPS_LOST;
          mode = enumStateEstimateMode::RUN_GPS_LOST;

          // '<S5>:67:6' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
          // '<S5>:67:7' bodyRatesOut_radps = filtBodyRatesIn_radps;
          // '<S5>:67:8' normMagVecOut_nd = normMagVecIn_nd;
          // '<S5>:67:9' isMagValid = isMagDataValid;
          isMagValid = rtu_magData->isMagDataValid;

          // '<S5>:67:10' latLonAltOut = latLonAltIn;
          bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
          stateEstimator_DW.bodyRatesOut_radps[0] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
          stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
          stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
          bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
          stateEstimator_DW.bodyRatesOut_radps[1] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
          stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
          stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
          bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
          stateEstimator_DW.bodyRatesOut_radps[2] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
          stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
          stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

          // '<S5>:67:11' isGpsValid = isGpsDataValid;
          isGpsValid = rtu_gpsData->isGpsDataValid;

          // '<S5>:67:12' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
          baroAltOut_m = stateEstimator_DW.Divide1 -
            stateEstimator_DW.baroInitAltMean;

          // '<S5>:67:13' isBaroValid = isBaroDataValid;
          isBaroValid = rtu_baroData->isBaroDataValid;

          //
        } else {
          // EKF RUNNING WITHOUT GPS BUT WE HAVE HEALTHY GPS SIGNAL
          // START INITIALIZING GPS
          // '<S5>:64:22' bodyAccelsOut_mps2 = filtBodyAccelsIn_mps2;
          // '<S5>:64:23' bodyRatesOut_radps = filtBodyRatesIn_radps;
          // '<S5>:64:24' normMagVecOut_nd = normMagVecIn_nd;
          // '<S5>:64:25' isMagValid = isMagDataValid;
          isMagValid = rtu_magData->isMagDataValid;

          // '<S5>:64:26' latLonAltOut = latLonAltIn;
          bodyAccelsOut_mps2[0] = stateEstimator_DW.Product[0];
          stateEstimator_DW.bodyRatesOut_radps[0] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[0];
          stateEstimator_DW.normMagVecOut_nd[0] = stateEstimator_DW.Divide[0];
          stateEstimator_DW.latLonAltOut[0] = rtu_gpsData->latLonAlt[0];
          bodyAccelsOut_mps2[1] = stateEstimator_DW.Product[1];
          stateEstimator_DW.bodyRatesOut_radps[1] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[1];
          stateEstimator_DW.normMagVecOut_nd[1] = stateEstimator_DW.Divide[1];
          stateEstimator_DW.latLonAltOut[1] = rtu_gpsData->latLonAlt[1];
          bodyAccelsOut_mps2[2] = stateEstimator_DW.Product[2];
          stateEstimator_DW.bodyRatesOut_radps[2] =
            stateEstimator_DW.TmpSignalConversionAtSFunctionI[2];
          stateEstimator_DW.normMagVecOut_nd[2] = stateEstimator_DW.Divide[2];
          stateEstimator_DW.latLonAltOut[2] = rtu_gpsData->latLonAlt[2];

          // '<S5>:64:27' isGpsValid = isGpsDataValid;
          isGpsValid = rtu_gpsData->isGpsDataValid;

          // '<S5>:64:28' baroAltOut_m = baroPressAlt_m - baroInitAltMean;
          baroAltOut_m = stateEstimator_DW.Divide1 -
            stateEstimator_DW.baroInitAltMean;

          // '<S5>:64:29' isBaroValid = isBaroDataValid;
          isBaroValid = rtu_baroData->isBaroDataValid;

          //
          // Compute running mean of GPS data for NED origin Lat, Lon and Alt
          // '<S5>:64:32' if (isGpsDataValid)
          if (rtu_gpsData->isGpsDataValid) {
            // '<S5>:64:33' if( gpsIdx < max(stateEstSmParams.gpsInitCount, 1) ) 
            if (stateEstimator_DW.gpsIdx < std::fmax
                (rtu_stateEstSmParams->gpsInitCount, 1.0F)) {
              // '<S5>:64:34' gpsIdx = gpsIdx + 1;
              stateEstimator_DW.gpsIdx++;

              // '<S5>:64:35' llhDelta = (latLonAltIn - refLatLonAlt);
              // '<S5>:64:36' refLatLonAlt = refLatLonAlt + llhDelta/gpsIdx;
              stateEstimator_DW.refLatLonAlt[0] += (rtu_gpsData->latLonAlt[0] -
                stateEstimator_DW.refLatLonAlt[0]) / stateEstimator_DW.gpsIdx;
              stateEstimator_DW.refLatLonAlt[1] += (rtu_gpsData->latLonAlt[1] -
                stateEstimator_DW.refLatLonAlt[1]) / stateEstimator_DW.gpsIdx;
              stateEstimator_DW.refLatLonAlt[2] += (rtu_gpsData->latLonAlt[2] -
                stateEstimator_DW.refLatLonAlt[2]) / stateEstimator_DW.gpsIdx;
            }

            //
            // '<S5>:64:39' if(gpsIdx >= stateEstSmParams.gpsInitCount)
            if (stateEstimator_DW.gpsIdx >= rtu_stateEstSmParams->gpsInitCount)
            {
              // '<S5>:64:40' isPosInitialized = true;
              stateEstimator_DW.isPosInitialized = true;
            }
          }
        }
      }
      break;
    }
  }

  if (static_cast<boolean_T>(rtu_gpsData->isGpsDataValid ^ 1)) {
    stateEstimator_DW.durationCounter_1++;
    stateEstimator_DW.durationCounter_1_c++;
    stateEstimator_DW.durationCounter_1_p++;
  } else {
    stateEstimator_DW.durationCounter_1 = 0;
    stateEstimator_DW.durationCounter_1_c = 0;
    stateEstimator_DW.durationCounter_1_p = 0;
  }

  // Delay: '<S1>/Delay'
  stateEstimator_DW.icLoad = stateEstimator_DW.resetStates |
    stateEstimator_DW.icLoad;
  if (stateEstimator_DW.icLoad) {
    std::memcpy(&stateEstimator_DW.Delay_DSTATE[0],
                &stateEstimator_DW.initialStates[0], 23U * sizeof(real32_T));
  }

  // Delay: '<S1>/Delay2'
  stateEstimator_DW.icLoad_j = stateEstimator_DW.resetStates |
    stateEstimator_DW.icLoad_j;
  if (stateEstimator_DW.icLoad_j) {
    stateEstimator_DW.Delay2_DSTATE[0] = stateEstimator_DW.initialDcmBodyToNed[0];
    stateEstimator_DW.Delay2_DSTATE[1] = stateEstimator_DW.initialDcmBodyToNed[1];
    stateEstimator_DW.Delay2_DSTATE[2] = stateEstimator_DW.initialDcmBodyToNed[2];
    stateEstimator_DW.Delay2_DSTATE[3] = stateEstimator_DW.initialDcmBodyToNed[3];
    stateEstimator_DW.Delay2_DSTATE[4] = stateEstimator_DW.initialDcmBodyToNed[4];
    stateEstimator_DW.Delay2_DSTATE[5] = stateEstimator_DW.initialDcmBodyToNed[5];
    stateEstimator_DW.Delay2_DSTATE[6] = stateEstimator_DW.initialDcmBodyToNed[6];
    stateEstimator_DW.Delay2_DSTATE[7] = stateEstimator_DW.initialDcmBodyToNed[7];
    stateEstimator_DW.Delay2_DSTATE[8] = stateEstimator_DW.initialDcmBodyToNed[8];
  }

  // RelationalOperator: '<S25>/Compare' incorporates:
  //   Constant: '<S25>/Constant'

  rtb_Compare = (mode == enumStateEstimateMode::RUN);

  // MATLAB Function: '<S7>/convertLlhToNedPos'
  // MATLAB Function 'latLonAltToNedPos/convertLlhToNedPos': '<S27>:1'
  // '<S27>:1:3' nedPos_m = convertLlhToNedPos_function(latLonAlt, refLatLonAlt, isGpsValid); 
  //  CONVERTLLHTONEDPOS_FUNCTION converts Latitude, Longitude, and Height (LLH) to 
  //  North-East-Down (nedPos_m) coordinates with the origin at refLlh.
  //
  //  Inputs:
  //  llh: [latRad, lonRad, height] in radians and meters
  //  refLlh: reference [latRad, lonRad, height] in radians and meters
  //  isGpsValid: Flag that indicates if gps is valid or not
  //
  //  Outputs:
  //  nedPos_m: [north, east, down] position in meters
  // 'convertLlhToNedPos_function:13' if(~isGpsValid)
  if (static_cast<boolean_T>(rtb_Compare ^ 1)) {
    // 'convertLlhToNedPos_function:14' nedPos_m = double([0; 0; 0]);
    rtb_nedPos_m_idx_0 = 0.0;
    rtb_nedPos_m_idx_1 = 0.0;
    rtb_nedPos_m_idx_2 = 0.0;
  } else {
    real_T n;
    real_T nRef;
    real_T n_idx_0;
    real_T n_idx_0_tmp;
    real_T rtb_nedPos_m_tmp;
    real_T rtb_nedPos_m_tmp_0;

    //  Constants
    // 'convertLlhToNedPos_function:19' semiMajorAxis = double(6378137.0);
    //  WGS84 semi-major axis (meters)
    // 'convertLlhToNedPos_function:20' flattening = double(1 / 298.257223563);
    //  WGS84 flattening
    // 'convertLlhToNedPos_function:21' eccentricitySquared = 2 * flattening - flattening^2; 
    //  Square of the eccentricity
    //  Extract input values
    // 'convertLlhToNedPos_function:24' latRad = llh(1);
    // 'convertLlhToNedPos_function:25' lonRad = llh(2);
    // 'convertLlhToNedPos_function:26' height = llh(3);
    // 'convertLlhToNedPos_function:28' refLatRad = refLlh(1);
    // 'convertLlhToNedPos_function:29' refLonRad = refLlh(2);
    // 'convertLlhToNedPos_function:30' refHeight = refLlh(3);
    //  Calculate the prime vertical radius of curvature at the reference point
    // 'convertLlhToNedPos_function:33' nRef = semiMajorAxis / sqrt(1 - eccentricitySquared * sin(refLatRad)^2); 
    nRef = std::sin(stateEstimator_DW.refLatLonAlt[0]);
    nRef = 6.378137E+6 / std::sqrt(1.0 - nRef * nRef * 0.0066943799901413165);

    //  Calculate ECEF coordinates of the reference point
    // 'convertLlhToNedPos_function:36' refX = (nRef + refHeight) * cos(refLatRad) * cos(refLonRad); 
    // 'convertLlhToNedPos_function:37' refY = (nRef + refHeight) * cos(refLatRad) * sin(refLonRad); 
    // 'convertLlhToNedPos_function:38' refZ = ((nRef * (1 - eccentricitySquared)) + refHeight) * sin(refLatRad); 
    //  Calculate the prime vertical radius of curvature at the current point
    // 'convertLlhToNedPos_function:41' n = semiMajorAxis / sqrt(1 - eccentricitySquared * sin(latRad)^2); 
    rtb_nedPos_m_idx_1 = std::sin(stateEstimator_DW.latLonAltOut[0]);
    n = 6.378137E+6 / std::sqrt(1.0 - rtb_nedPos_m_idx_1 * rtb_nedPos_m_idx_1 *
      0.0066943799901413165);

    //  Calculate ECEF coordinates of the current point
    // 'convertLlhToNedPos_function:44' x = (n + height) * cos(latRad) * cos(lonRad); 
    // 'convertLlhToNedPos_function:45' y = (n + height) * cos(latRad) * sin(lonRad); 
    // 'convertLlhToNedPos_function:46' z = ((n * (1 - eccentricitySquared)) + height) * sin(latRad); 
    //  Compute the ECEF displacement from the reference point
    // 'convertLlhToNedPos_function:49' dx = x - refX;
    // 'convertLlhToNedPos_function:50' dy = y - refY;
    // 'convertLlhToNedPos_function:51' dz = z - refZ;
    //  Compute the rotation matrix from ECEF to nedPos_m
    // 'convertLlhToNedPos_function:54' R = [-sin(refLatRad) * cos(refLonRad), -sin(refLatRad) * sin(refLonRad), cos(refLatRad); 
    // 'convertLlhToNedPos_function:55'      -sin(refLonRad), cos(refLonRad), 0; 
    // 'convertLlhToNedPos_function:56'      -cos(refLatRad) * cos(refLonRad), -cos(refLatRad) * sin(refLonRad), -sin(refLatRad)]; 
    //  Rotate the displacement vector to nedPos_m coordinates
    // 'convertLlhToNedPos_function:59' nedPos_m = R * [dx; dy; dz];
    rtb_nedPos_m_idx_2 = std::sin(stateEstimator_DW.refLatLonAlt[1]);
    rtb_nedPos_m_tmp = std::cos(stateEstimator_DW.refLatLonAlt[1]);
    rtb_nedPos_m_tmp_0 = std::cos(stateEstimator_DW.refLatLonAlt[0]);
    rtb_nedPos_m_idx_0 = std::sin(stateEstimator_DW.refLatLonAlt[0]);
    n_idx_0_tmp = (n + stateEstimator_DW.latLonAltOut[2]) * std::cos
      (stateEstimator_DW.latLonAltOut[0]);
    n_idx_0 = n_idx_0_tmp * std::cos(stateEstimator_DW.latLonAltOut[1]) - (nRef
      + stateEstimator_DW.refLatLonAlt[2]) * rtb_nedPos_m_tmp_0 *
      rtb_nedPos_m_tmp;
    n_idx_0_tmp = n_idx_0_tmp * std::sin(stateEstimator_DW.latLonAltOut[1]) -
      (nRef + stateEstimator_DW.refLatLonAlt[2]) * std::cos
      (stateEstimator_DW.refLatLonAlt[0]) * rtb_nedPos_m_idx_2;
    nRef = (n * 0.99330562000985867 + stateEstimator_DW.latLonAltOut[2]) *
      rtb_nedPos_m_idx_1 - (nRef * 0.99330562000985867 +
      stateEstimator_DW.refLatLonAlt[2]) * rtb_nedPos_m_idx_0;

    // 'convertLlhToNedPos_function:61' for idx = 1:3
    rtb_nedPos_m_idx_0 = (-rtb_nedPos_m_idx_0 * rtb_nedPos_m_tmp * n_idx_0 +
                          -std::sin(stateEstimator_DW.refLatLonAlt[0]) *
                          rtb_nedPos_m_idx_2 * n_idx_0_tmp) + rtb_nedPos_m_tmp_0
      * nRef;

    // 'convertLlhToNedPos_function:62' if( abs(nedPos_m(idx)) < 1e-7 )
    if (std::abs(rtb_nedPos_m_idx_0) < 1.0E-7) {
      // 'convertLlhToNedPos_function:63' nedPos_m(idx) = double(0);
      rtb_nedPos_m_idx_0 = 0.0;
    }

    rtb_nedPos_m_idx_1 = -rtb_nedPos_m_idx_2 * n_idx_0 + rtb_nedPos_m_tmp *
      n_idx_0_tmp;

    // 'convertLlhToNedPos_function:62' if( abs(nedPos_m(idx)) < 1e-7 )
    if (std::abs(rtb_nedPos_m_idx_1) < 1.0E-7) {
      // 'convertLlhToNedPos_function:63' nedPos_m(idx) = double(0);
      rtb_nedPos_m_idx_1 = 0.0;
    }

    rtb_nedPos_m_idx_2 = (-rtb_nedPos_m_tmp_0 * rtb_nedPos_m_tmp * n_idx_0 +
                          -std::cos(stateEstimator_DW.refLatLonAlt[0]) *
                          rtb_nedPos_m_idx_2 * n_idx_0_tmp) + -std::sin
      (stateEstimator_DW.refLatLonAlt[0]) * nRef;

    // 'convertLlhToNedPos_function:62' if( abs(nedPos_m(idx)) < 1e-7 )
    if (std::abs(rtb_nedPos_m_idx_2) < 1.0E-7) {
      // 'convertLlhToNedPos_function:63' nedPos_m(idx) = double(0);
      rtb_nedPos_m_idx_2 = 0.0;
    }
  }

  // End of MATLAB Function: '<S7>/convertLlhToNedPos'

  // Product: '<S8>/Product' incorporates:
  //   Trigonometry: '<S8>/Cos'
  //   Trigonometry: '<S8>/Cos2'
  //   UnitDelay: '<Root>/Unit Delay1'

  rtb_Product2 = std::cos(stateEstimator_DW.UnitDelay1_DSTATE[0]) * std::cos
    (stateEstimator_DW.UnitDelay1_DSTATE[1]);

  // Sum: '<S8>/Sum1' incorporates:
  //   Gain: '<S8>/Gain'
  //   Product: '<S8>/Product1'
  //   Product: '<S8>/Product2'
  //   Product: '<S8>/Product3'
  //   Sum: '<S8>/Sum'
  //   Trigonometry: '<S8>/Cos1'
  //   UnitDelay: '<Root>/Unit Delay1'

  rtb_Product2 = rtb_Product2 * rtu_lidarData->range_m -
    (rtu_lidarParams->yMntOff_m * -std::sin(stateEstimator_DW.UnitDelay1_DSTATE
      [0]) + rtu_lidarParams->zMntOff_m * rtb_Product2);

  // Logic: '<S8>/AND1' incorporates:
  //   Logic: '<S8>/AND'
  //   Logic: '<S8>/NOT'
  //   Logic: '<S8>/OR'
  //   RelationalOperator: '<S8>/Less Than'
  //   RelationalOperator: '<S8>/Less Than1'

  rtb_AND1 = static_cast<boolean_T>((rtu_lidarData->range_m >=
    rtu_lidarParams->validRange_m[0]) & (rtu_lidarData->range_m <=
    rtu_lidarParams->validRange_m[1])) & static_cast<boolean_T>
    (rtu_lidarData->isLidarDataValid & rtu_lidarData->isLidarInitialized);

  // MATLAB Function: '<S1>/EKF' incorporates:
  //   Delay: '<S1>/Delay'
  //   Delay: '<S1>/Delay2'

  rtb_dcmBodyToNed_k_idx_0 = stateEstimator_DW.Delay2_DSTATE[0];
  rtb_dcmBodyToNed_k_idx_1 = stateEstimator_DW.Delay2_DSTATE[1];
  rtb_dcmBodyToNed_k_idx_2 = stateEstimator_DW.Delay2_DSTATE[2];
  rtb_dcmBodyToNed_k_idx_3 = stateEstimator_DW.Delay2_DSTATE[3];
  rtb_dcmBodyToNed_k_idx_4 = stateEstimator_DW.Delay2_DSTATE[4];
  rtb_dcmBodyToNed_k_idx_5 = stateEstimator_DW.Delay2_DSTATE[5];
  rtb_dcmBodyToNed_k_idx_6 = stateEstimator_DW.Delay2_DSTATE[6];
  rtb_dcmBodyToNed_k_idx_7 = stateEstimator_DW.Delay2_DSTATE[7];
  rtb_dcmBodyToNed_k_idx_8 = stateEstimator_DW.Delay2_DSTATE[8];

  // MATLAB Function 'EKF/EKF': '<S11>:1'
  //   function [states, covP, dcmBodyToNed] = ekf(bodyAccelsIn_mps2, bodyRatesIn_radps,  ... 
  //      normMagVec_nd, isMagValid, nedPosAndVel, isGpsValid, baroAlt_m, isBaroValid, ... 
  //      lidarAgl_m, isLidarValid, prevStates, prevCovP, prevDcmBodyToNed, estSmMode, ... 
  //      processNoiseQ, measNoiseR, gEarth_mps2, sampleTime_s)
  //
  //  [states, covP, dcmBodyToNed] = ekf_function(bodyAccelsIn_mps2, bodyRatesIn_radps,  ... 
  //      normMagVec_nd, isMagValid, nedPosAndVel, isGpsValid, baroAlt_m, isBaroValid, ... 
  //      lidarAgl_m, isLidarValid, prevStates, prevCovP, prevDcmBodyToNed, estSmMode, ... 
  //      processNoiseQ, measNoiseR, gEarth_mps2, sampleTime_s);
  //   end
  // '<S11>:1:18' if isempty(covP) || resetStates
  if (static_cast<boolean_T>(static_cast<boolean_T>
       (stateEstimator_DW.covP_not_empty_p ^ 1) | stateEstimator_DW.resetStates))
  {
    // '<S11>:1:19' covP = initCovP;
    std::memcpy(&stateEstimator_DW.covP_a[0], &rtu_initCovP[0], 529U * sizeof
                (real32_T));
    stateEstimator_DW.covP_not_empty_p = true;
  }

  // '<S11>:1:22' [states, covP, dcmBodyToNed] = ekf_function(bodyAccelsIn_mps2, bodyRatesIn_radps,  ... 
  // '<S11>:1:23'     normMagVec_nd, isMagValid, nedPosAndVel, isGpsValid, baroAlt_m, isBaroValid, ... 
  // '<S11>:1:24'     lidarAgl_m, isLidarValid, prevStates, covP, prevDcmBodyToNed, estSmMode, ... 
  // '<S11>:1:25'     processNoiseQ, measNoiseR, gEarth_mps2, sampleTime_s);
  // EKF runs an EKF to estimate required states
  //
  // Inputs:
  // bodyAccels_mps2:           Body accels measured usig accelerometer
  // bodyRates_radps:           Body angular rates measured using Gyro
  // normMagVec_nd:             Mag data from onboard magnetometer as unit
  // vector
  // isMagValid:                Boolean flag to indicate if mag data is valid
  // nedPosAndVel:              Measured Vehicle NED position and velocity
  // isGpsValid:                Boolean flag to indicate if GPS position is
  // valid
  // baroAlt_m:                 Baro altitude
  // isBaroValid:               Boolean flag to indicate if baro data is valid
  // lidarAgl_m:                Lidar AGL
  // isLidarValid:              Boolean flag to indicate if lidar data is valid
  // prevStates:                Previous state estimate
  // prevCovP:                  Previous covariance
  // prevDcmBodyToNed:          Body to NED DCM computed using prevStates
  // quaternion
  // estSmMode:                 Estimator state machine mode
  // processNoiseQ:             Process Noise Matrix
  // measNoiseR:                Meas Noise Matrix
  // gEarth_mps2:               Acceleration due to gravity
  // sampleTime_s:              Sample Time
  // ekfParams:                 EKF params struct
  //
  // Ouputs:
  // states:                    Current states
  // cov:                       Current covariance
  // dcmBodyToNed:              Body to NED DCM computed using current states
  // quaternion
  // Propagate state
  // 'ekf_function:42' if (estSmMode == enumStateEstimateMode.INITIALIZE || ...
  // 'ekf_function:43'     estSmMode ~= enumStateEstimateMode.RUN)
  if (static_cast<boolean_T>((mode == enumStateEstimateMode::INITIALIZE) | (mode
        != enumStateEstimateMode::RUN))) {
    // 'ekf_function:44' states = prevStates;
    std::memcpy(&rty_states[0], &stateEstimator_DW.Delay_DSTATE[0], 23U * sizeof
                (real32_T));
  } else {
    // 'ekf_function:48' states = updateEkfStates(prevStates, bodyAccels_mps2, bodyRates_radps, ... 
    // 'ekf_function:49'     dcmBodyToNed, sampleTime_s, gEarth_mps2);
    // UPDATEEKFSTATES propogates the state
    //
    // Inputs:
    // states:                Previous state estimate
    // bodyAccels_mps2:           Body accels measured usig accelerometer
    // bodyRates_radps:           Body angular rates measured using Gyro
    // dcmBodyToNed:              Body to NED DCM computed using states
    // quaternion
    // baroData:                  Baro altitude
    // sampleTime_s:              Sample time for integration
    // gEarth_mps2:               Accel due to gravity
    //
    // Ouputs:
    // states:                    Current states
    // Quaternions
    // 'updateEkfStates:19' q0 = states(1);
    // 'updateEkfStates:20' q1 = states(2);
    // 'updateEkfStates:21' q2 = states(3);
    // 'updateEkfStates:22' q3 = states(4);
    // 'updateEkfStates:24' vN = states(8);
    // 'updateEkfStates:25' vE = states(9);
    // 'updateEkfStates:26' vD = states(10);
    // Gyro biases
    // 'updateEkfStates:29' bwx = states(11);
    // 'updateEkfStates:30' bwy = states(12);
    // 'updateEkfStates:31' bwz = states(13);
    // Accel biases
    // 'updateEkfStates:34' bax = states(14);
    // 'updateEkfStates:35' bay = states(15);
    // 'updateEkfStates:36' baz = states(16);
    // Angular rate measurement from gyros
    // 'updateEkfStates:39' wx = bodyRates_radps(1);
    //  rad/s
    // 'updateEkfStates:40' wy = bodyRates_radps(2);
    // 'updateEkfStates:41' wz = bodyRates_radps(3);
    // Specific force measurement from accelerometers
    // 'updateEkfStates:44' ax = bodyAccels_mps2(1);
    //  m/s^2
    // 'updateEkfStates:45' ay = bodyAccels_mps2(2);
    // 'updateEkfStates:46' az = bodyAccels_mps2(3);
    // d(quaternion)/dt = C_bodyrate2qdot*[wx; wy; wz]
    // 'updateEkfStates:49' C_bodyrates2qdot  = 0.5*[-q1  -q2  -q3; ...
    // 'updateEkfStates:50'     q0  -q3   q2; ...
    // 'updateEkfStates:51'     q3   q0  -q1; ...
    // 'updateEkfStates:52'     -q2   q1   q0];
    // 'updateEkfStates:54' stateDot = [ ...
    // 'updateEkfStates:55'     C_bodyrates2qdot*([wx;wy;wz]-[bwx;bwy;bwz]); ...                            % Derivative of [q0; q1; q2; q3] 
    // 'updateEkfStates:56'     [vN; vE; vD]; ...                                                           % Derivative of [pN; pE; pD] 
    // 'updateEkfStates:57'     dcmBodyToNed * ([ax;ay;az]-[bax;bay;baz]) + [0; 0; gEarth_mps2]; ...        % Derivative of [vN; vE; vD] 
    // 'updateEkfStates:58'     [0;0;0]; ...                                                                % Derivative of [bwx; bwy; bwz] 
    // 'updateEkfStates:59'     [0;0;0]; ...                                                                % Derivative of [bax; bay; baz] 
    // 'updateEkfStates:60'     [0;0;0]; ...                                                                % Derivative of normalized NED mag  
    // 'updateEkfStates:61'     [0;0;0]; ...                                                                % Derivative of [bmx; bmy; bmz] 
    // 'updateEkfStates:62'      0; ...                                                                     % Derivative of Baro bias 
    // 'updateEkfStates:63'     ];
    //                             % Derivative of [q0; q1; q2; q3]
    //                                                            % Derivative of [pN; pE; pD] 
    //         % Derivative of [vN; vE; vD]
    //                                                                 % Derivative of [bwx; bwy; bwz] 
    //                                                                 % Derivative of [bax; bay; baz] 
    //                                                                 % Derivative of normalized NED mag 
    //                                                                 % Derivative of [bmx; bmy; bmz] 
    //                                                                      % Derivative of Baro bias 
    // 'updateEkfStates:65' states = states + sampleTime_s * stateDot;
    rtb_XAxis = 0.5F * -stateEstimator_DW.Delay_DSTATE[1];
    tmp[0] = rtb_XAxis;
    rtb_XAxis1 = 0.5F * -stateEstimator_DW.Delay_DSTATE[2];
    tmp[4] = rtb_XAxis1;
    rtb_XAxis2 = 0.5F * -stateEstimator_DW.Delay_DSTATE[3];
    tmp[8] = rtb_XAxis2;
    tmp[1] = 0.5F * stateEstimator_DW.Delay_DSTATE[0];
    tmp[5] = rtb_XAxis2;
    tmp[9] = 0.5F * stateEstimator_DW.Delay_DSTATE[2];
    tmp[2] = 0.5F * stateEstimator_DW.Delay_DSTATE[3];
    tmp[6] = 0.5F * stateEstimator_DW.Delay_DSTATE[0];
    tmp[10] = rtb_XAxis;
    tmp[3] = rtb_XAxis1;
    tmp[7] = 0.5F * stateEstimator_DW.Delay_DSTATE[1];
    tmp[11] = 0.5F * stateEstimator_DW.Delay_DSTATE[0];
    rtb_XAxis1 = stateEstimator_DW.bodyRatesOut_radps[0] -
      stateEstimator_DW.Delay_DSTATE[10];
    rtb_XAxis2 = stateEstimator_DW.bodyRatesOut_radps[1] -
      stateEstimator_DW.Delay_DSTATE[11];
    rtb_UnitDelay_e = stateEstimator_DW.bodyRatesOut_radps[2] -
      stateEstimator_DW.Delay_DSTATE[12];
    for (i_0 = 0; i_0 < 4; i_0++) {
      tmp_0[i_0] = 0.0F;
      tmp_0[i_0] += tmp[i_0] * rtb_XAxis1;
      tmp_0[i_0] += tmp[i_0 + 4] * rtb_XAxis2;
      tmp_0[i_0] += tmp[i_0 + 8] * rtb_UnitDelay_e;
    }

    rtb_XAxis1 = bodyAccelsOut_mps2[0] - stateEstimator_DW.Delay_DSTATE[13];
    rtb_XAxis2 = bodyAccelsOut_mps2[1] - stateEstimator_DW.Delay_DSTATE[14];
    rtb_UnitDelay_e = bodyAccelsOut_mps2[2] - stateEstimator_DW.Delay_DSTATE[15];
    rty_states[0] = 0.004F * tmp_0[0] + stateEstimator_DW.Delay_DSTATE[0];
    rty_states[1] = 0.004F * tmp_0[1] + stateEstimator_DW.Delay_DSTATE[1];
    rty_states[2] = 0.004F * tmp_0[2] + stateEstimator_DW.Delay_DSTATE[2];
    rty_states[3] = 0.004F * tmp_0[3] + stateEstimator_DW.Delay_DSTATE[3];
    rty_states[4] = 0.004F * stateEstimator_DW.Delay_DSTATE[7] +
      stateEstimator_DW.Delay_DSTATE[4];
    rty_states[5] = 0.004F * stateEstimator_DW.Delay_DSTATE[8] +
      stateEstimator_DW.Delay_DSTATE[5];
    rty_states[6] = 0.004F * stateEstimator_DW.Delay_DSTATE[9] +
      stateEstimator_DW.Delay_DSTATE[6];
    rty_states[7] = ((stateEstimator_DW.Delay2_DSTATE[0] * rtb_XAxis1 +
                      stateEstimator_DW.Delay2_DSTATE[3] * rtb_XAxis2) +
                     stateEstimator_DW.Delay2_DSTATE[6] * rtb_UnitDelay_e) *
      0.004F + stateEstimator_DW.Delay_DSTATE[7];
    rty_states[10] = stateEstimator_DW.Delay_DSTATE[10];
    rty_states[13] = stateEstimator_DW.Delay_DSTATE[13];
    rty_states[16] = stateEstimator_DW.Delay_DSTATE[16];
    rty_states[19] = stateEstimator_DW.Delay_DSTATE[19];
    rty_states[8] = ((stateEstimator_DW.Delay2_DSTATE[1] * rtb_XAxis1 +
                      stateEstimator_DW.Delay2_DSTATE[4] * rtb_XAxis2) +
                     stateEstimator_DW.Delay2_DSTATE[7] * rtb_UnitDelay_e) *
      0.004F + stateEstimator_DW.Delay_DSTATE[8];
    rty_states[11] = stateEstimator_DW.Delay_DSTATE[11];
    rty_states[14] = stateEstimator_DW.Delay_DSTATE[14];
    rty_states[17] = stateEstimator_DW.Delay_DSTATE[17];
    rty_states[20] = stateEstimator_DW.Delay_DSTATE[20];
    rty_states[9] = (((stateEstimator_DW.Delay2_DSTATE[2] * rtb_XAxis1 +
                       stateEstimator_DW.Delay2_DSTATE[5] * rtb_XAxis2) +
                      stateEstimator_DW.Delay2_DSTATE[8] * rtb_UnitDelay_e) +
                     *rtu_gEarth_mps2) * 0.004F +
      stateEstimator_DW.Delay_DSTATE[9];
    rty_states[12] = stateEstimator_DW.Delay_DSTATE[12];
    rty_states[15] = stateEstimator_DW.Delay_DSTATE[15];
    rty_states[18] = stateEstimator_DW.Delay_DSTATE[18];
    rty_states[21] = stateEstimator_DW.Delay_DSTATE[21];
    rty_states[22] = stateEstimator_DW.Delay_DSTATE[22];

    // Normalize the quaternion
    // 'ekf_function:52' states(1:4) = states(1:4)/norm(states(1:4));
    rtb_XAxis = norm_7MzYkgry(&rty_states[0]);
    rty_states[0] /= rtb_XAxis;
    rty_states[1] /= rtb_XAxis;
    rty_states[2] /= rtb_XAxis;
    rty_states[3] /= rtb_XAxis;

    // Propogate covariances
    // 'ekf_function:55' stateJac = computeStateJac(prevStates, bodyAccels_mps2, bodyRates_radps, sampleTime_s); 
    //  covP = stateJac * covP * stateJac' + processNoiseQ;
    // 'ekf_function:58' covP = updateCovP(covP, stateJac, processNoiseQ);
    computeStateJac_YJuzBPAK(stateEstimator_DW.Delay_DSTATE, bodyAccelsOut_mps2,
      stateEstimator_DW.bodyRatesOut_radps, 0.004F, tmp_1);
    updateCovP_UFEbdQNU(stateEstimator_DW.covP_a, tmp_1, rtu_processNoiseQ);

    //  covP = propagateCovInplace(prevCovP, stateJac, processNoiseQ);
    // Fuse Mag data if it is valid
    // 'ekf_function:62' if(isMagValid)
    if (isMagValid) {
      real32_T tmp19;

      // 'ekf_function:63' measJac = computeMagMeasJac(states);
      // COMPUTEMAGMEASJAC computes the measurement jacobian for magnetometer
      // states in the estimator
      //
      // Inputs:
      // states:                EKF states
      //
      // Outputs:
      // magMeasJac:            3x23 Mag Meas Jacobian
      // 'computeMagMeasJac:11' magMeasJac = zeros(3, 23, 'single');
      std::memset(&measJac[0], 0, 69U * sizeof(real32_T));

      // Extract quat states
      // 'computeMagMeasJac:14' q0 = states(1);
      // 'computeMagMeasJac:15' q1 = states(2);
      // 'computeMagMeasJac:16' q2 = states(3);
      // 'computeMagMeasJac:17' q3 = states(4);
      // Extract Est NED Mag
      // 'computeMagMeasJac:19' magN = states(17);
      // 'computeMagMeasJac:20' magE = states(18);
      // 'computeMagMeasJac:21' magD = states(19);
      // 'computeMagMeasJac:22' tmp1 = 2*magD;
      tmp1 = 2.0F * rty_states[18];

      // 'computeMagMeasJac:23' tmp2 = q2*tmp1;
      tmp2 = rty_states[2] * tmp1;

      // 'computeMagMeasJac:24' tmp3 = q3*tmp1;
      rtb_Product1 = rty_states[3] * tmp1;

      // 'computeMagMeasJac:25' tmp4 = 2*magE;
      rtb_XAxis = 2.0F * rty_states[17];

      // 'computeMagMeasJac:26' tmp5 = q2*tmp4;
      rtb_XAxis1 = rty_states[2] * rtb_XAxis;

      // 'computeMagMeasJac:27' tmp6 = q0*tmp1;
      tmp6 = rty_states[0] * tmp1;

      // 'computeMagMeasJac:28' tmp7 = q1*tmp4;
      rtb_XAxis2 = rty_states[1] * rtb_XAxis;

      // 'computeMagMeasJac:29' tmp8 = 4*magN;
      tmp8 = 4.0F * rty_states[16];

      // 'computeMagMeasJac:30' tmp9 = q1*tmp1;
      tmp9 = rty_states[1] * tmp1;

      // 'computeMagMeasJac:31' tmp10 = q0*tmp4;
      rtb_UnitDelay_e = rty_states[0] * rtb_XAxis;

      // 'computeMagMeasJac:32' tmp11 = 2*q2^2;
      tmp11 = rty_states[2] * rty_states[2] * 2.0F;

      // 'computeMagMeasJac:33' tmp12 = 2*q3^2 - 1;
      rtb_dcmBodyToNed_k_idx_5 = rty_states[3] * rty_states[3] * 2.0F - 1.0F;

      // 'computeMagMeasJac:34' tmp13 = 2*q0;
      rtb_dcmBodyToNed_k_idx_3 = 2.0F * rty_states[0];

      // 'computeMagMeasJac:35' tmp14 = q3*tmp13;
      rtb_dcmBodyToNed_k_idx_8 = rty_states[3] * rtb_dcmBodyToNed_k_idx_3;

      // 'computeMagMeasJac:36' tmp15 = 2*q1;
      tmp1 = 2.0F * rty_states[1];

      // 'computeMagMeasJac:37' tmp16 = q2*tmp13;
      rtb_dcmBodyToNed_k_idx_6 = rty_states[2] * rtb_dcmBodyToNed_k_idx_3;

      // 'computeMagMeasJac:38' tmp17 = 2*magN;
      tmp17 = 2.0F * rty_states[16];

      // 'computeMagMeasJac:39' tmp18 = -q3*tmp17;
      rtb_dcmBodyToNed_k_idx_1 = -rty_states[3] * tmp17;

      // 'computeMagMeasJac:40' tmp19 = 4*magE;
      tmp19 = 4.0F * rty_states[17];

      // 'computeMagMeasJac:41' tmp20 = magN*tmp15;
      rtb_dcmBodyToNed_k_idx_4 = rty_states[16] * tmp1;

      // 'computeMagMeasJac:42' tmp21 = magN*tmp13;
      rtb_dcmBodyToNed_k_idx_7 = rty_states[16] * rtb_dcmBodyToNed_k_idx_3;

      // 'computeMagMeasJac:43' tmp22 = 2*q1^2;
      rtb_dcmBodyToNed_k_idx_0 = rty_states[1] * rty_states[1];
      rtb_dcmBodyToNed_k_idx_2 = rtb_dcmBodyToNed_k_idx_0 * 2.0F;

      // 'computeMagMeasJac:44' tmp23 = q1*tmp13;
      rtb_dcmBodyToNed_k_idx_3 *= rty_states[1];

      // 'computeMagMeasJac:45' magMeasJac(1, 1) = 2*magE*q3 - tmp2;
      measJac[0] = 2.0F * rty_states[17] * rty_states[3] - tmp2;

      // 'computeMagMeasJac:46' magMeasJac(1, 2) = tmp3 + tmp5;
      measJac[3] = rtb_Product1 + rtb_XAxis1;

      // 'computeMagMeasJac:47' magMeasJac(1, 3) = -q2*tmp8 - tmp6 + tmp7;
      measJac[6] = (-rty_states[2] * tmp8 - tmp6) + rtb_XAxis2;

      // 'computeMagMeasJac:48' magMeasJac(1, 4) = -q3*tmp8 + tmp10 + tmp9;
      measJac[9] = (-rty_states[3] * tmp8 + rtb_UnitDelay_e) + tmp9;

      // 'computeMagMeasJac:49' magMeasJac(1, 17) = -tmp11 - tmp12;
      measJac[48] = -tmp11 - rtb_dcmBodyToNed_k_idx_5;

      // 'computeMagMeasJac:50' magMeasJac(1, 18) = q2*tmp15 + tmp14;
      measJac[51] = rty_states[2] * tmp1 + rtb_dcmBodyToNed_k_idx_8;

      // 'computeMagMeasJac:51' magMeasJac(1, 19) = 2*q1*q3 - tmp16;
      measJac[54] = 2.0F * rty_states[1] * rty_states[3] -
        rtb_dcmBodyToNed_k_idx_6;

      // 'computeMagMeasJac:52' magMeasJac(1, 20) = 1;
      measJac[57] = 1.0F;

      // 'computeMagMeasJac:54' magMeasJac(2, 1) = tmp18 + tmp9;
      measJac[1] = rtb_dcmBodyToNed_k_idx_1 + tmp9;

      // 'computeMagMeasJac:55' magMeasJac(2, 2) = -q1*tmp19 + q2*tmp17 + tmp6;
      measJac[4] = (-rty_states[1] * tmp19 + rty_states[2] * tmp17) + tmp6;

      // 'computeMagMeasJac:56' magMeasJac(2, 3) = tmp20 + tmp3;
      measJac[7] = rtb_dcmBodyToNed_k_idx_4 + rtb_Product1;

      // 'computeMagMeasJac:57' magMeasJac(2, 4) = -q3*tmp19 + tmp2 - tmp21;
      measJac[10] = (-rty_states[3] * tmp19 + tmp2) - rtb_dcmBodyToNed_k_idx_7;

      // 'computeMagMeasJac:58' magMeasJac(2, 17) = 2*q1*q2 - tmp14;
      measJac[49] = 2.0F * rty_states[1] * rty_states[2] -
        rtb_dcmBodyToNed_k_idx_8;

      // 'computeMagMeasJac:59' magMeasJac(2, 18) = -tmp12 - tmp22;
      measJac[52] = -rtb_dcmBodyToNed_k_idx_5 - rtb_dcmBodyToNed_k_idx_2;

      // 'computeMagMeasJac:60' magMeasJac(2, 19) = 2*q2*q3 + tmp23;
      rtb_dcmBodyToNed_k_idx_5 = 2.0F * rty_states[2] * rty_states[3];
      measJac[55] = rtb_dcmBodyToNed_k_idx_5 + rtb_dcmBodyToNed_k_idx_3;

      // 'computeMagMeasJac:61' magMeasJac(2, 21) = 1;
      measJac[61] = 1.0F;

      // 'computeMagMeasJac:63' tmp24 = 4*magD;
      tmp2 = 4.0F * rty_states[18];

      // 'computeMagMeasJac:64' magMeasJac(3, 1) = 2*magN*q2 - tmp7;
      measJac[2] = 2.0F * rty_states[16] * rty_states[2] - rtb_XAxis2;

      // 'computeMagMeasJac:65' magMeasJac(3, 2) = -q1*tmp24 - tmp10 - tmp18;
      measJac[5] = (-rty_states[1] * tmp2 - rtb_UnitDelay_e) -
        rtb_dcmBodyToNed_k_idx_1;

      // 'computeMagMeasJac:66' magMeasJac(3, 3) = -q2*tmp24 + q3*tmp4 + tmp21;
      measJac[8] = (-rty_states[2] * tmp2 + rty_states[3] * rtb_XAxis) +
        rtb_dcmBodyToNed_k_idx_7;

      // 'computeMagMeasJac:67' magMeasJac(3, 4) = tmp20 + tmp5;
      measJac[11] = rtb_dcmBodyToNed_k_idx_4 + rtb_XAxis1;

      // 'computeMagMeasJac:68' magMeasJac(3, 17) = q3*tmp15 + tmp16;
      measJac[50] = rty_states[3] * tmp1 + rtb_dcmBodyToNed_k_idx_6;

      // 'computeMagMeasJac:69' magMeasJac(3, 18) = 2*q2*q3 - tmp23;
      measJac[53] = rtb_dcmBodyToNed_k_idx_5 - rtb_dcmBodyToNed_k_idx_3;

      // 'computeMagMeasJac:70' magMeasJac(3, 19) = -tmp11 - tmp22 + 1;
      measJac[56] = (-tmp11 - rtb_dcmBodyToNed_k_idx_2) + 1.0F;

      // 'computeMagMeasJac:71' magMeasJac(3, 22) = 1;
      measJac[65] = 1.0F;

      // 'ekf_function:64' [states, covP] = fuseMagData(states, covP, normMagVec_nd, measJac, measNoiseR); 
      // FUSEMAGDATA fuses Magnetometer data in the EKF
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // bodyMagUnitVec:    Normalized magnetometer measurement in body frame
      // measJac:           Measurement Jacobian
      // measNoiseR:        Magnetometer measurement noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // Quaternions
      // 'fuseMagData:16' q0 = states(1);
      // 'fuseMagData:17' q1 = states(2);
      // 'fuseMagData:18' q2 = states(3);
      // 'fuseMagData:19' q3 = states(4);
      // NED Mag Norm and biases
      // 'fuseMagData:22' nedMagUnitVec = states(17:19);
      //  Direction Cosine Matrix (DCM) from NED coordinates to body cooridinates 
      //  expressed using quaternions, using the current state estimate.
      // 'fuseMagData:26' C_ned2b  = [ 1-2*(q2^2+q3^2)    2*(q1*q2+q3*q0)     2*(q1*q3-q2*q0); ... 
      // 'fuseMagData:27'              2*(q1*q2-q3*q0)    1-2*(q1^2+q3^2)     2*(q2*q3+q1*q0); ... 
      // 'fuseMagData:28'              2*(q1*q3+q2*q0)    2*(q2*q3-q1*q0)     1-2*(q1^2+q2^2)]; 
      // Rotate propogated mag states and add bias to estimate measurements
      // 'fuseMagData:31' estBodyMagUnitVec = C_ned2b*nedMagUnitVec + states(20:22); 
      // 'fuseMagData:33' tmp1 = covP*measJac';
      for (i_0 = 0; i_0 < 23; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 3; i++) {
          b_tmp1_tmp = f_y_tmp_0 + i_0;
          b_tmp1[b_tmp1_tmp] = 0.0F;
          i_1 = 0;
          covP_tmp = 0;
          for (d_tmp1_tmp = 0; d_tmp1_tmp < 23; d_tmp1_tmp++) {
            b_tmp1[b_tmp1_tmp] += stateEstimator_DW.covP_a[i_1 + i_0] *
              measJac[covP_tmp + i];
            i_1 += 23;
            covP_tmp += 3;
          }

          f_y_tmp_0 += 23;
        }
      }

      // 'fuseMagData:35' K = tmp1/(measJac * tmp1 + measNoiseR(1:3, 1:3));
      for (i_0 = 0; i_0 < 3; i_0++) {
        f_y_tmp_0 = 0;
        i = 0;
        i_1 = 0;
        for (covP_tmp = 0; covP_tmp < 3; covP_tmp++) {
          rtb_XAxis = 0.0F;
          d_tmp1_tmp = 0;
          for (b_tmp1_tmp = 0; b_tmp1_tmp < 23; b_tmp1_tmp++) {
            rtb_XAxis += measJac[d_tmp1_tmp + i_0] * b_tmp1[b_tmp1_tmp + i_1];
            d_tmp1_tmp += 3;
          }

          measJac_1[f_y_tmp_0 + i_0] = rtu_measNoiseR[i + i_0] + rtb_XAxis;
          f_y_tmp_0 += 3;
          i += 11;
          i_1 += 23;
        }
      }

      mrdiv_8PldARhB(b_tmp1, measJac_1, K);

      // 'fuseMagData:36' states = states + K*(bodyMagUnitVec - estBodyMagUnitVec); 
      rtb_XAxis2 = rty_states[3] * rty_states[3];
      tmp6 = rty_states[1] * rty_states[2];
      tmp8 = rty_states[0] * rty_states[3];
      rtb_XAxis = rty_states[1] * rty_states[3];
      rtb_UnitDelay_e = rty_states[0] * rty_states[2];
      tmp1 = rty_states[2] * rty_states[2];
      rtb_XAxis1 = stateEstimator_DW.normMagVecOut_nd[0] - ((((1.0F - (tmp1 +
        rtb_XAxis2) * 2.0F) * rty_states[16] + (tmp6 + tmp8) * 2.0F *
        rty_states[17]) + (rtb_XAxis - rtb_UnitDelay_e) * 2.0F * rty_states[18])
        + rty_states[19]);
      tmp2 = rty_states[2] * rty_states[3];
      rtb_Product1 = rty_states[0] * rty_states[1];
      rtb_XAxis2 = stateEstimator_DW.normMagVecOut_nd[1] - ((((1.0F -
        (rtb_dcmBodyToNed_k_idx_0 + rtb_XAxis2) * 2.0F) * rty_states[17] + (tmp6
        - tmp8) * 2.0F * rty_states[16]) + (tmp2 + rtb_Product1) * 2.0F *
        rty_states[18]) + rty_states[20]);
      rtb_UnitDelay_e = stateEstimator_DW.normMagVecOut_nd[2] - ((((rtb_XAxis +
        rtb_UnitDelay_e) * 2.0F * rty_states[16] + (tmp2 - rtb_Product1) * 2.0F *
        rty_states[17]) + (1.0F - (rtb_dcmBodyToNed_k_idx_0 + tmp1) * 2.0F) *
        rty_states[18]) + rty_states[21]);

      // 'fuseMagData:37' covP = covP - K * measJac*covP;
      for (i_0 = 0; i_0 < 23; i_0++) {
        rtb_dcmBodyToNed_k_idx_0 = K[i_0 + 23];
        rtb_dcmBodyToNed_k_idx_3 = K[i_0 + 46];
        c_tmp1[i_0] = ((rtb_dcmBodyToNed_k_idx_0 * rtb_XAxis2 + K[i_0] *
                        rtb_XAxis1) + rtb_dcmBodyToNed_k_idx_3 * rtb_UnitDelay_e)
          + rty_states[i_0];
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          i = 23 * f_y_tmp_0 + i_0;
          K_1[i] = 0.0F;
          K_1[i] += measJac[3 * f_y_tmp_0] * K[i_0];
          K_1[i] += measJac[3 * f_y_tmp_0 + 1] * rtb_dcmBodyToNed_k_idx_0;
          K_1[i] += measJac[3 * f_y_tmp_0 + 2] * rtb_dcmBodyToNed_k_idx_3;
        }

        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          rtb_XAxis = 0.0F;
          for (i = 0; i < 23; i++) {
            rtb_XAxis += K_1[23 * i + i_0] * stateEstimator_DW.covP_a[23 *
              f_y_tmp_0 + i];
          }

          covP_tmp = 23 * f_y_tmp_0 + i_0;
          covP[covP_tmp] = stateEstimator_DW.covP_a[covP_tmp] - rtb_XAxis;
        }
      }

      // 'fuseMagData:38' covP = (covP + covP')/2;
      i = 0;
      for (i_1 = 0; i_1 < 23; i_1++) {
        i_0 = 0;
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          covP_tmp = f_y_tmp_0 + i;
          stateEstimator_DW.covP_a[covP_tmp] = (covP[i_0 + i_1] + covP[covP_tmp])
            / 2.0F;
          i_0 += 23;
        }

        rty_states[i_1] = c_tmp1[i_1];
        i += 23;
      }

      // Normalize the quaternion
      // 'ekf_function:66' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_XAxis = norm_7MzYkgry(&c_tmp1[0]);
      rty_states[0] = c_tmp1[0] / rtb_XAxis;
      rty_states[1] = c_tmp1[1] / rtb_XAxis;
      rty_states[2] = c_tmp1[2] / rtb_XAxis;
      rty_states[3] = c_tmp1[3] / rtb_XAxis;
    }

    // Fuse GPS data if it is valid
    // 'ekf_function:70' if(isGpsValid)
    if (isGpsValid) {
      // 'ekf_function:71' measJac = zeros(6, 23, 'single');
      // 'ekf_function:73' measJac(1, 5) = single(1);
      // 'ekf_function:74' measJac(2, 6) = single(1);
      // 'ekf_function:75' measJac(3, 7) = single(1);
      // 'ekf_function:76' measJac(4, 8) = single(1);
      // 'ekf_function:77' measJac(5, 9) = single(1);
      // 'ekf_function:78' measJac(6, 10) = single(1);
      // 'ekf_function:80' [states, covP] = fuseGpsData(states, covP, nedPosAndVel, measJac, measNoiseR); 
      // FUSEGPSDATA Fuses GPS data in EKF
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // nedPosVel:         NED Position and velocity from GPS
      // measJac:           Measurement Jacobian
      // measNoiseR:        GPS measurement noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // 'fuseGpsData:15' gpsMeasIdx = [4, 5, 6, 7, 8, 9];
      // 'fuseGpsData:16' gpsStateIdx = [5, 6, 7, 8, 9, 10];
      // 'fuseGpsData:17' gpsSnsIdx = [1, 2, 3, 4, 5, 6];
      // 'fuseGpsData:19' tmp1 = covP * measJac';
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 6; f_y_tmp_0++) {
        for (i = 0; i < 23; i++) {
          d_tmp1_tmp = i + i_0;
          d_tmp1[d_tmp1_tmp] = 0.0F;
          i_1 = 0;
          for (covP_tmp = 0; covP_tmp < 23; covP_tmp++) {
            d_tmp1[d_tmp1_tmp] += stateEstimator_DW.covP_a[i_1 + i] *
              static_cast<real32_T>(e[covP_tmp + i_0]);
            i_1 += 23;
          }
        }

        i_0 += 23;
      }

      // 'fuseGpsData:20' K = tmp1/(measJac * tmp1 + measNoiseR(gpsMeasIdx, gpsMeasIdx)); 
      b_K_tmp_0 = &d_measJac[0];
      for (i_0 = 0; i_0 < 6; i_0++) {
        f_y_tmp_0 = 0;
        i = 0;
        i_1 = 0;
        for (covP_tmp = 0; covP_tmp < 6; covP_tmp++) {
          rtb_XAxis = 0.0F;
          d_tmp1_tmp = 0;
          for (b_tmp1_tmp = 0; b_tmp1_tmp < 23; b_tmp1_tmp++) {
            rtb_XAxis += static_cast<real32_T>(b_K_tmp_0[d_tmp1_tmp + i_0]) *
              d_tmp1[b_tmp1_tmp + i_1];
            d_tmp1_tmp += 6;
          }

          b_K_tmp[f_y_tmp_0 + i_0] = rtu_measNoiseR[(i + i_0) + 36] + rtb_XAxis;
          f_y_tmp_0 += 6;
          i += 11;
          i_1 += 23;
        }
      }

      mrdiv_0UppHIxu(d_tmp1, b_K_tmp, b_K);

      // SignalConversion generated from: '<S11>/ SFunction ' incorporates:
      //   DataTypeConversion: '<S7>/Cast To Single'
      //   Sum: '<S7>/Sum'
      //   UnitDelay: '<S7>/Unit Delay'

      // 'fuseGpsData:21' states = states + K*(nedPosVel(gpsSnsIdx) - states(gpsStateIdx)); 
      tmp1 = static_cast<real32_T>(rtb_nedPos_m_idx_0) - rty_states[4];
      rtb_XAxis1 = static_cast<real32_T>(rtb_nedPos_m_idx_1) - rty_states[5];
      tmp2 = (static_cast<real32_T>(rtb_nedPos_m_idx_2) +
              stateEstimator_DW.UnitDelay_DSTATE_c) - rty_states[6];
      rtb_dcmBodyToNed_k_idx_6 = rtu_gpsData->nedVel_mps[0] - rty_states[7];
      rtb_dcmBodyToNed_k_idx_1 = rtu_gpsData->nedVel_mps[1] - rty_states[8];
      rtb_dcmBodyToNed_k_idx_4 = rtu_gpsData->nedVel_mps[2] - rty_states[9];

      // 'fuseGpsData:22' covP = covP - K*measJac*covP;
      for (i_0 = 0; i_0 < 23; i_0++) {
        rtb_dcmBodyToNed_k_idx_0 = b_K[i_0 + 23];
        rtb_dcmBodyToNed_k_idx_3 = b_K[i_0 + 46];
        rtb_dcmBodyToNed_k_idx_7 = b_K[i_0 + 69];
        rtb_dcmBodyToNed_k_idx_2 = b_K[i_0 + 92];
        rtb_dcmBodyToNed_k_idx_5 = b_K[i_0 + 115];
        c_tmp1[i_0] = (((((rtb_dcmBodyToNed_k_idx_0 * rtb_XAxis1 + b_K[i_0] *
                           tmp1) + rtb_dcmBodyToNed_k_idx_3 * tmp2) +
                         rtb_dcmBodyToNed_k_idx_7 * rtb_dcmBodyToNed_k_idx_6) +
                        rtb_dcmBodyToNed_k_idx_2 * rtb_dcmBodyToNed_k_idx_1) +
                       rtb_dcmBodyToNed_k_idx_5 * rtb_dcmBodyToNed_k_idx_4) +
          rty_states[i_0];
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          i = 23 * f_y_tmp_0 + i_0;
          K_1[i] = 0.0F;
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0]) * b_K[i_0];
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0 + 1]) *
            rtb_dcmBodyToNed_k_idx_0;
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0 + 2]) *
            rtb_dcmBodyToNed_k_idx_3;
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0 + 3]) *
            rtb_dcmBodyToNed_k_idx_7;
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0 + 4]) *
            rtb_dcmBodyToNed_k_idx_2;
          K_1[i] += static_cast<real32_T>(b_K_tmp_0[6 * f_y_tmp_0 + 5]) *
            rtb_dcmBodyToNed_k_idx_5;
        }

        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          rtb_XAxis = 0.0F;
          for (i = 0; i < 23; i++) {
            rtb_XAxis += K_1[23 * i + i_0] * stateEstimator_DW.covP_a[23 *
              f_y_tmp_0 + i];
          }

          covP_tmp = 23 * f_y_tmp_0 + i_0;
          covP[covP_tmp] = stateEstimator_DW.covP_a[covP_tmp] - rtb_XAxis;
        }
      }

      // 'fuseGpsData:23' covP = (covP + covP')/2;
      i = 0;
      for (i_1 = 0; i_1 < 23; i_1++) {
        i_0 = 0;
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
          covP_tmp = f_y_tmp_0 + i;
          stateEstimator_DW.covP_a[covP_tmp] = (covP[i_0 + i_1] + covP[covP_tmp])
            / 2.0F;
          i_0 += 23;
        }

        rty_states[i_1] = c_tmp1[i_1];
        i += 23;
      }

      // Normalize the quaternion
      // 'ekf_function:82' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&c_tmp1[0]);
      rty_states[0] = c_tmp1[0] / rtb_Product1;
      rty_states[1] = c_tmp1[1] / rtb_Product1;
      rty_states[2] = c_tmp1[2] / rtb_Product1;
      rty_states[3] = c_tmp1[3] / rtb_Product1;
    }

    // Fuse Baro data if it is valid
    // 'ekf_function:86' if(isBaroValid)
    if (isBaroValid) {
      // 'ekf_function:87' measJac = zeros(1, 23, 'single');
      // 'ekf_function:88' measJac(1, 7) = -1;
      // 'ekf_function:89' measJac(1, 23) = 1;
      // 'ekf_function:90' [states, covP] = fuseBaroData(states, covP, baroAlt_m, measJac, measNoiseR); 
      // FUSEBARODATA Fuses BARO data in EKF
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // baroAlt_m:         Baro alt data
      // measJac:           Measurement Jacobian
      // measNoiseR:        Baro alt meas noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // 'fuseBaroData:15' tmp1 = covP * measJac';
      // 'fuseBaroData:16' K = tmp1/(measJac * tmp1 + measNoiseR(10, 10));
      rtb_dcmBodyToNed_k_idx_0 = 0.0F;
      for (i_0 = 0; i_0 < 23; i_0++) {
        c_tmp1[i_0] = 0.0F;
        f_y_tmp_0 = 0;
        for (i = 0; i < 23; i++) {
          c_tmp1[i_0] += stateEstimator_DW.covP_a[f_y_tmp_0 + i_0] *
            static_cast<real32_T>(b_measJac[i]);
          f_y_tmp_0 += 23;
        }

        f_y_tmp_0 = b_measJac[i_0];
        rtb_dcmBodyToNed_k_idx_0 += static_cast<real32_T>(f_y_tmp_0) *
          c_tmp1[i_0];
        f_y_tmp[i_0] = static_cast<int8_T>(f_y_tmp_0);
      }

      rtb_Product1 = rtb_dcmBodyToNed_k_idx_0 + rtu_measNoiseR[108];
      for (i_0 = 0; i_0 < 23; i_0++) {
        c_tmp1[i_0] /= rtb_Product1;
      }

      // 'fuseBaroData:17' states = states + K*(baroAlt_m  + states(7) - states(23)); 
      rtb_Product1 = (baroAltOut_m + rty_states[6]) - rty_states[22];

      // 'fuseBaroData:18' covP = covP - K*measJac*covP;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
        rty_states[f_y_tmp_0] += c_tmp1[f_y_tmp_0] * rtb_Product1;
        for (i = 0; i < 23; i++) {
          K_1[i + i_0] = c_tmp1[i] * static_cast<real32_T>(f_y_tmp[f_y_tmp_0]);
        }

        i_0 += 23;
      }

      for (i_0 = 0; i_0 < 23; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 23; i++) {
          rtb_XAxis = 0.0F;
          i_1 = 0;
          for (covP_tmp = 0; covP_tmp < 23; covP_tmp++) {
            rtb_XAxis += K_1[i_1 + i_0] * stateEstimator_DW.covP_a[covP_tmp +
              f_y_tmp_0];
            i_1 += 23;
          }

          covP_tmp = f_y_tmp_0 + i_0;
          covP[covP_tmp] = stateEstimator_DW.covP_a[covP_tmp] - rtb_XAxis;
          f_y_tmp_0 += 23;
        }
      }

      // 'fuseBaroData:19' covP = (covP + covP')/2;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
        i = 0;
        for (i_1 = 0; i_1 < 23; i_1++) {
          covP_tmp = i_1 + i_0;
          stateEstimator_DW.covP_a[covP_tmp] = (covP[i + f_y_tmp_0] +
            covP[covP_tmp]) / 2.0F;
          i += 23;
        }

        i_0 += 23;
      }

      // Normalize the quaternion
      // 'ekf_function:92' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&rty_states[0]);
      rty_states[0] /= rtb_Product1;
      rty_states[1] /= rtb_Product1;
      rty_states[2] /= rtb_Product1;
      rty_states[3] /= rtb_Product1;
    }

    // Fuse Lidar data if it is valid
    // 'ekf_function:96' if(isLidarValid)
    if (rtb_AND1) {
      // 'ekf_function:97' measJac = zeros(1, 23, 'single');
      // 'ekf_function:98' measJac(1, 7) = -1;
      // 'ekf_function:99' [states, covP] = fuseLidarData(states, covP, lidarAgl_m, measJac, measNoiseR); 
      // FUSELIDARDATA Fuses LIDAR data in EKF
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // lidarAgl_m:        Lidar agl data
      // measJac:           Measurement Jacobian
      // measNoiseR:        Baro alt meas noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // 'fuseLidarData:15' tmp1 = covP*measJac';
      // 'fuseLidarData:16' K = tmp1/(measJac * tmp1 + measNoiseR(11, 11));
      rtb_dcmBodyToNed_k_idx_0 = 0.0F;
      for (i_0 = 0; i_0 < 23; i_0++) {
        c_tmp1[i_0] = 0.0F;
        f_y_tmp_0 = 0;
        for (i = 0; i < 23; i++) {
          c_tmp1[i_0] += stateEstimator_DW.covP_a[f_y_tmp_0 + i_0] *
            static_cast<real32_T>(c_measJac[i]);
          f_y_tmp_0 += 23;
        }

        f_y_tmp_0 = c_measJac[i_0];
        rtb_dcmBodyToNed_k_idx_0 += static_cast<real32_T>(f_y_tmp_0) *
          c_tmp1[i_0];
        f_y_tmp[i_0] = static_cast<int8_T>(f_y_tmp_0);
      }

      rtb_Product1 = rtb_dcmBodyToNed_k_idx_0 + rtu_measNoiseR[120];
      for (i_0 = 0; i_0 < 23; i_0++) {
        c_tmp1[i_0] /= rtb_Product1;
      }

      // 'fuseLidarData:17' states = states + K*(lidarAgl_m  + states(7));
      rtb_Product1 = rtb_Product2 + rty_states[6];

      // 'fuseLidarData:18' covP = covP - K*measJac*covP;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
        rty_states[f_y_tmp_0] += c_tmp1[f_y_tmp_0] * rtb_Product1;
        for (i = 0; i < 23; i++) {
          K_1[i + i_0] = c_tmp1[i] * static_cast<real32_T>(f_y_tmp[f_y_tmp_0]);
        }

        i_0 += 23;
      }

      for (i_0 = 0; i_0 < 23; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 23; i++) {
          rtb_XAxis = 0.0F;
          i_1 = 0;
          for (covP_tmp = 0; covP_tmp < 23; covP_tmp++) {
            rtb_XAxis += K_1[i_1 + i_0] * stateEstimator_DW.covP_a[covP_tmp +
              f_y_tmp_0];
            i_1 += 23;
          }

          covP_tmp = f_y_tmp_0 + i_0;
          covP[covP_tmp] = stateEstimator_DW.covP_a[covP_tmp] - rtb_XAxis;
          f_y_tmp_0 += 23;
        }
      }

      // 'fuseLidarData:19' covP = (covP + covP')/2;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 23; f_y_tmp_0++) {
        i = 0;
        for (i_1 = 0; i_1 < 23; i_1++) {
          covP_tmp = i_1 + i_0;
          stateEstimator_DW.covP_a[covP_tmp] = (covP[i + f_y_tmp_0] +
            covP[covP_tmp]) / 2.0F;
          i += 23;
        }

        i_0 += 23;
      }

      // Normalize the quaternion
      // 'ekf_function:101' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&rty_states[0]);
      rty_states[0] /= rtb_Product1;
      rty_states[1] /= rtb_Product1;
      rty_states[2] /= rtb_Product1;
      rty_states[3] /= rtb_Product1;
    }

    // Compute Body To NED DCM
    // 'ekf_function:105' dcmBodyToNed = quatToDcm_function(states(1:4));
    // Quaternions
    // 'quatToDcm_function:3' q0 = quat(1);
    // 'quatToDcm_function:4' q1 = quat(2);
    // 'quatToDcm_function:5' q2 = quat(3);
    // 'quatToDcm_function:6' q3 = quat(4);
    //  Direction Cosine Matrix (DCM) from body cooridinates to NED coordinates
    //  expressed using quaternions.
    // 'quatToDcm_function:10' dcmBodyToNed = [1-2*(q2^2+q3^2), 2*(q1*q2-q3*q0), 2*(q1*q3+q2*q0); 
    // 'quatToDcm_function:11'     2*(q1*q2+q3*q0), 1-2*(q1^2+q3^2), 2*(q2*q3-q1*q0); 
    // 'quatToDcm_function:12'     2*(q1*q3-q2*q0), 2*(q2*q3+q1*q0), 1-2*(q1^2+q2^2)]; 
    rtb_dcmBodyToNed_k_idx_4 = rty_states[3] * rty_states[3];
    rtb_dcmBodyToNed_k_idx_8 = rty_states[2] * rty_states[2];
    rtb_dcmBodyToNed_k_idx_0 = 1.0F - (rtb_dcmBodyToNed_k_idx_8 +
      rtb_dcmBodyToNed_k_idx_4) * 2.0F;
    rtb_dcmBodyToNed_k_idx_1 = rty_states[1] * rty_states[2];
    rtb_dcmBodyToNed_k_idx_7 = rty_states[0] * rty_states[3];
    rtb_dcmBodyToNed_k_idx_3 = (rtb_dcmBodyToNed_k_idx_1 -
      rtb_dcmBodyToNed_k_idx_7) * 2.0F;
    rtb_dcmBodyToNed_k_idx_2 = rty_states[1] * rty_states[3];
    rtb_dcmBodyToNed_k_idx_5 = rty_states[0] * rty_states[2];
    rtb_dcmBodyToNed_k_idx_6 = (rtb_dcmBodyToNed_k_idx_2 +
      rtb_dcmBodyToNed_k_idx_5) * 2.0F;
    rtb_dcmBodyToNed_k_idx_1 = (rtb_dcmBodyToNed_k_idx_1 +
      rtb_dcmBodyToNed_k_idx_7) * 2.0F;
    rtb_XAxis1 = rty_states[1] * rty_states[1];
    rtb_dcmBodyToNed_k_idx_4 = 1.0F - (rtb_XAxis1 + rtb_dcmBodyToNed_k_idx_4) *
      2.0F;
    rtb_XAxis = rty_states[2] * rty_states[3];
    rtb_XAxis2 = rty_states[0] * rty_states[1];
    rtb_dcmBodyToNed_k_idx_7 = (rtb_XAxis - rtb_XAxis2) * 2.0F;
    rtb_dcmBodyToNed_k_idx_2 = (rtb_dcmBodyToNed_k_idx_2 -
      rtb_dcmBodyToNed_k_idx_5) * 2.0F;
    rtb_dcmBodyToNed_k_idx_5 = (rtb_XAxis + rtb_XAxis2) * 2.0F;
    rtb_dcmBodyToNed_k_idx_8 = 1.0F - (rtb_XAxis1 + rtb_dcmBodyToNed_k_idx_8) *
      2.0F;
  }

  // End of MATLAB Function: '<S1>/EKF'

  // MATLAB Function: '<S1>/EKF NO GPS' incorporates:
  //   Delay: '<S1>/Delay'
  //   Delay: '<S1>/Delay2'

  rtb_XAxis1 = stateEstimator_DW.Delay2_DSTATE[0];
  rtb_UnitDelay_e = stateEstimator_DW.Delay2_DSTATE[1];
  rtb_Product1 = stateEstimator_DW.Delay2_DSTATE[2];
  rtb_XAxis = stateEstimator_DW.Delay2_DSTATE[3];
  tmp1 = stateEstimator_DW.Delay2_DSTATE[4];
  tmp8 = stateEstimator_DW.Delay2_DSTATE[5];
  rtb_XAxis2 = stateEstimator_DW.Delay2_DSTATE[6];
  tmp2 = stateEstimator_DW.Delay2_DSTATE[7];
  tmp6 = stateEstimator_DW.Delay2_DSTATE[8];

  // MATLAB Function 'EKF/EKF NO GPS': '<S12>:1'
  // '<S12>:1:7' if isempty(covP) || resetStates
  if (static_cast<boolean_T>(static_cast<boolean_T>
       (stateEstimator_DW.covP_not_empty ^ 1) | stateEstimator_DW.resetStates))
  {
    // '<S12>:1:8' covP = initCovP;
    std::memcpy(&stateEstimator_DW.covP[0], &rtu_initCovNoGpsP[0], 324U * sizeof
                (real32_T));
    stateEstimator_DW.covP_not_empty = true;
  }

  // '<S12>:1:11' noGpsEkfStatesIdx = [1, 2, 3, 4, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; 
  // '<S12>:1:12' [states, covP, dcmBodyToNed] = ekfNoGps_function(bodyAccelsIn_mps2, bodyRatesIn_radps,  ... 
  // '<S12>:1:13'     normMagVec_nd, isMagValid, baroAlt_m, isBaroValid, lidarAgl_m, isLidarValid, ... 
  // '<S12>:1:14'     prevStates(noGpsEkfStatesIdx), covP, prevDcmBodyToNed, estSmMode, processNoiseQ, ... 
  // '<S12>:1:15'     measNoiseR, gEarth_mps2, sampleTime_s);
  for (i_0 = 0; i_0 < 18; i_0++) {
    d_tmp1_0[i_0] = stateEstimator_DW.Delay_DSTATE[c[i_0]];
  }

  // EKF runs an EKF to estimate required states
  //
  // Inputs:
  // bodyAccels_mps2:           Body accels measured usig accelerometer
  // bodyRates_radps:           Body angular rates measured using Gyro
  // normMagVec_nd:             Mag data from onboard magnetometer as unit
  // vector
  // isMagValid:                Boolean flag to indicate if mag data is valid
  // valid
  // baroAlt_m:                 Baro altitude
  // isBaroValid:               Boolean flag to indicate if baro data is valid
  // lidarAgl_m:                Lidar AGL
  // isLidarValid:              Boolean flag to indicate if lidar data is valid
  // prevStates:                Previous state estimate
  // covP:                      Previous covariance
  // dcmBodyToNed:              Body to NED DCM computed using prevStates
  // quaternion
  // estSmMode:                 Estimator state machine mode
  // processNoiseQ:             Process Noise Matrix
  // measNoiseR:                Meas Noise Matrix
  // gEarth_mps2:               Acceleration due to gravity
  // sampleTime_s:              Sample Time
  // ekfParams:                 EKF params struct
  //
  // Ouputs:
  // states:                    Current states
  // cov:                       Current covariance
  // dcmBodyToNed:              Body to NED DCM computed using current states
  // quaternion
  // Propagate state
  // 'ekfNoGps_function:37' if (estSmMode == enumStateEstimateMode.INITIALIZE || ... 
  // 'ekfNoGps_function:38'     estSmMode == enumStateEstimateMode.RUN)
  if (static_cast<boolean_T>((mode == enumStateEstimateMode::INITIALIZE) | (mode
        == enumStateEstimateMode::RUN))) {
    // 'ekfNoGps_function:39' states = prevStates;
    std::memcpy(&rtb_states_a[0], &d_tmp1_0[0], 18U * sizeof(real32_T));
  } else {
    // 'ekfNoGps_function:43' states = updateEkfStatesNoGps(prevStates, bodyRates_radps, ... 
    // 'ekfNoGps_function:44'                               sampleTime_s);
    // UPDATEEKFSTATESNOGPS propogates the state for an EKF that doesn't have
    // access to the GPS
    //
    // Inputs:
    // states:                Previous state estimate
    // bodyRates_radps:           Body angular rates measured using Gyro
    // sampleTime_s:              Sample time for integration
    //
    // Ouputs:
    // states:                    Current states
    // Quaternions
    // 'updateEkfStatesNoGps:15' q0 = states(1);
    // 'updateEkfStatesNoGps:16' q1 = states(2);
    // 'updateEkfStatesNoGps:17' q2 = states(3);
    // 'updateEkfStatesNoGps:18' q3 = states(4);
    // Gyro biases
    // 'updateEkfStatesNoGps:21' bwx = states(6);
    // 'updateEkfStatesNoGps:22' bwy = states(7);
    // 'updateEkfStatesNoGps:23' bwz = states(8);
    // Angular rate measurement from gyros
    // 'updateEkfStatesNoGps:26' wx = bodyRates_radps(1);
    //  rad/s
    // 'updateEkfStatesNoGps:27' wy = bodyRates_radps(2);
    // 'updateEkfStatesNoGps:28' wz = bodyRates_radps(3);
    // d(quaternion)/dt = C_bodyrate2qdot*[wx; wy; wz]
    // 'updateEkfStatesNoGps:31' C_bodyrates2qdot  = 0.5*[-q1  -q2  -q3; ...
    // 'updateEkfStatesNoGps:32'     q0  -q3   q2; ...
    // 'updateEkfStatesNoGps:33'     q3   q0  -q1; ...
    // 'updateEkfStatesNoGps:34'     -q2   q1   q0];
    // 'updateEkfStatesNoGps:36' stateDot = [ ...
    // 'updateEkfStatesNoGps:37'     C_bodyrates2qdot*([wx;wy;wz]-[bwx;bwy;bwz]); ...                            % Derivative of [q0; q1; q2; q3] 
    // 'updateEkfStatesNoGps:38'     0; ...                                                                      % Down Position 
    // 'updateEkfStatesNoGps:39'     [0;0;0]; ...                                                                % Derivative of [bwx; bwy; bwz] 
    // 'updateEkfStatesNoGps:40'     [0;0;0]; ...                                                                % Derivative of [bax; bay; baz] 
    // 'updateEkfStatesNoGps:41'     [0;0;0]; ...                                                                % Derivative of normalized NED mag  
    // 'updateEkfStatesNoGps:42'     [0;0;0]; ...                                                                % Derivative of [bmx; bmy; bmz] 
    // 'updateEkfStatesNoGps:43'      0; ...                                                                     % Derivative of Baro bias 
    // 'updateEkfStatesNoGps:44'     ];
    //                             % Derivative of [q0; q1; q2; q3]
    //                                                                       % Down Position 
    //                                                                 % Derivative of [bwx; bwy; bwz] 
    //                                                                 % Derivative of [bax; bay; baz] 
    //                                                                 % Derivative of normalized NED mag 
    //                                                                 % Derivative of [bmx; bmy; bmz] 
    //                                                                      % Derivative of Baro bias 
    // 'updateEkfStatesNoGps:46' states = states + sampleTime_s * stateDot;
    rtb_XAxis1 = 0.5F * -stateEstimator_DW.Delay_DSTATE[1];
    tmp[0] = rtb_XAxis1;
    rtb_UnitDelay_e = 0.5F * -stateEstimator_DW.Delay_DSTATE[2];
    tmp[4] = rtb_UnitDelay_e;
    tmp2 = 0.5F * -stateEstimator_DW.Delay_DSTATE[3];
    tmp[8] = tmp2;
    rtb_Product1 = 0.5F * stateEstimator_DW.Delay_DSTATE[0];
    tmp[1] = rtb_Product1;
    tmp[5] = tmp2;
    tmp[9] = 0.5F * stateEstimator_DW.Delay_DSTATE[2];
    tmp[2] = 0.5F * stateEstimator_DW.Delay_DSTATE[3];
    tmp[6] = rtb_Product1;
    tmp[10] = rtb_XAxis1;
    tmp[3] = rtb_UnitDelay_e;
    tmp[7] = 0.5F * stateEstimator_DW.Delay_DSTATE[1];
    tmp[11] = rtb_Product1;
    rtb_XAxis1 = stateEstimator_DW.bodyRatesOut_radps[0] -
      stateEstimator_DW.Delay_DSTATE[10];
    rtb_XAxis2 = stateEstimator_DW.bodyRatesOut_radps[1] -
      stateEstimator_DW.Delay_DSTATE[11];
    rtb_UnitDelay_e = stateEstimator_DW.bodyRatesOut_radps[2] -
      stateEstimator_DW.Delay_DSTATE[12];
    for (i_0 = 0; i_0 < 4; i_0++) {
      tmp_0[i_0] = 0.0F;
      tmp_0[i_0] += tmp[i_0] * rtb_XAxis1;
      tmp_0[i_0] += tmp[i_0 + 4] * rtb_XAxis2;
      tmp_0[i_0] += tmp[i_0 + 8] * rtb_UnitDelay_e;
      rtb_states_a[i_0] = 0.004F * tmp_0[i_0] + d_tmp1_0[i_0];
    }

    rtb_states_a[4] = stateEstimator_DW.Delay_DSTATE[6];
    rtb_states_a[5] = d_tmp1_0[5];
    rtb_states_a[8] = d_tmp1_0[8];
    rtb_states_a[11] = d_tmp1_0[11];
    rtb_states_a[14] = d_tmp1_0[14];
    rtb_states_a[6] = d_tmp1_0[6];
    rtb_states_a[9] = d_tmp1_0[9];
    rtb_states_a[12] = d_tmp1_0[12];
    rtb_states_a[15] = d_tmp1_0[15];
    rtb_states_a[7] = d_tmp1_0[7];
    rtb_states_a[10] = d_tmp1_0[10];
    rtb_states_a[13] = d_tmp1_0[13];
    rtb_states_a[16] = d_tmp1_0[16];
    rtb_states_a[17] = stateEstimator_DW.Delay_DSTATE[22];

    // Normalize the quaternion
    // 'ekfNoGps_function:47' states(1:4) = states(1:4)/norm(states(1:4));
    rtb_XAxis = norm_7MzYkgry(&rtb_states_a[0]);
    rtb_states_a[0] /= rtb_XAxis;
    rtb_states_a[1] /= rtb_XAxis;
    rtb_states_a[2] /= rtb_XAxis;
    tmp17 = rtb_states_a[3] / rtb_XAxis;
    rtb_states_a[3] = tmp17;

    // Propogate covariances
    // 'ekfNoGps_function:50' stateJac = computeStateJacNoGps(prevStates, bodyRates_radps, sampleTime_s); 
    // COMPUTESTATEJACNOGPS Computes the NO GPS EKF state Jacobian
    //
    // Inputs:
    // states:             Time propogated states
    // bodyRates_radps:    Measured body rates
    // sampleTime_s:       Sample time
    //
    // Outputs:
    // stateJac:           State Jacobian
    // Initialize the state jacobian to zero
    //  stateJac = zeros(18, 18, 'single');
    // 'computeStateJacNoGps:14' stateJac = zeros(28, 1, 'single');
    // Extract quat states
    // 'computeStateJacNoGps:17' q0 = states(1);
    // 'computeStateJacNoGps:18' q1 = states(2);
    // 'computeStateJacNoGps:19' q2 = states(3);
    // 'computeStateJacNoGps:20' q3 = states(4);
    // Extract gyro biases from states
    // 'computeStateJacNoGps:22' bwx = states(6);
    // 'computeStateJacNoGps:23' bwy = states(7);
    // 'computeStateJacNoGps:24' bwz = states(8);
    // Extract gyro inputs
    // 'computeStateJacNoGps:26' wx = bodyRates_radps(1);
    // 'computeStateJacNoGps:27' wy = bodyRates_radps(2);
    // 'computeStateJacNoGps:28' wz = bodyRates_radps(3);
    // 'computeStateJacNoGps:30' tmp1 = sampleTime_s*(bwx/2 - wx/2);
    tmp1 = (stateEstimator_DW.Delay_DSTATE[10] / 2.0F -
            stateEstimator_DW.bodyRatesOut_radps[0] / 2.0F) * 0.004F;

    // 'computeStateJacNoGps:31' tmp2 = sampleTime_s*(bwy/2 - wy/2);
    tmp2 = (stateEstimator_DW.Delay_DSTATE[11] / 2.0F -
            stateEstimator_DW.bodyRatesOut_radps[1] / 2.0F) * 0.004F;

    // 'computeStateJacNoGps:32' tmp3 = sampleTime_s*(bwz/2 - wz/2);
    rtb_Product1 = (stateEstimator_DW.Delay_DSTATE[12] / 2.0F -
                    stateEstimator_DW.bodyRatesOut_radps[2] / 2.0F) * 0.004F;

    // 'computeStateJacNoGps:33' tmp4 = sampleTime_s/2;
    // 'computeStateJacNoGps:34' tmp5 = q1*tmp4;
    rtb_XAxis1 = stateEstimator_DW.Delay_DSTATE[1] * 0.002F;

    // 'computeStateJacNoGps:35' tmp6 = q2*tmp4;
    tmp6 = stateEstimator_DW.Delay_DSTATE[2] * 0.002F;

    // 'computeStateJacNoGps:36' tmp7 = q3*tmp4;
    rtb_XAxis2 = stateEstimator_DW.Delay_DSTATE[3] * 0.002F;

    // 'computeStateJacNoGps:37' tmp8 = -tmp4*(bwx - wx);
    tmp8 = (stateEstimator_DW.Delay_DSTATE[10] -
            stateEstimator_DW.bodyRatesOut_radps[0]) * -0.002F;

    // 'computeStateJacNoGps:38' tmp9 = -tmp4*(bwz - wz);
    tmp9 = (stateEstimator_DW.Delay_DSTATE[12] -
            stateEstimator_DW.bodyRatesOut_radps[2]) * -0.002F;

    // 'computeStateJacNoGps:39' tmp10 = -q0*tmp4;
    rtb_UnitDelay_e = -stateEstimator_DW.Delay_DSTATE[0] * 0.002F;

    // 'computeStateJacNoGps:40' tmp11 = -tmp4*(bwy - wy);
    tmp11 = (stateEstimator_DW.Delay_DSTATE[11] -
             stateEstimator_DW.bodyRatesOut_radps[1]) * -0.002F;

    // stateJac(1, 1) = 1;
    // 'computeStateJacNoGps:43' stateJac(1) = 1;
    stateJac[0] = 1.0F;

    // stateJac(1, 2) = tmp1;
    // 'computeStateJacNoGps:46' stateJac(2) = tmp1;
    stateJac[1] = tmp1;

    //  stateJac(1, 3) = tmp2;
    // 'computeStateJacNoGps:49' stateJac(3) = tmp2;
    stateJac[2] = tmp2;

    //  stateJac(1, 4) = tmp3;
    // 'computeStateJacNoGps:52' stateJac(4) = tmp3;
    stateJac[3] = rtb_Product1;

    //  stateJac(1, 6) = tmp5;
    // 'computeStateJacNoGps:55' stateJac(5) = tmp5;
    stateJac[4] = rtb_XAxis1;

    //  stateJac(1, 7) = tmp6;
    // 'computeStateJacNoGps:58' stateJac(6) = tmp6;
    stateJac[5] = tmp6;

    //  stateJac(1, 8) = tmp7;
    // 'computeStateJacNoGps:61' stateJac(7) = tmp7;
    stateJac[6] = rtb_XAxis2;

    //  stateJac(2, 1) = tmp8;
    // 'computeStateJacNoGps:64' stateJac(8) = tmp8;
    stateJac[7] = tmp8;

    //  stateJac(2, 2) = 1;
    // 'computeStateJacNoGps:67' stateJac(9) = 1;
    stateJac[8] = 1.0F;

    //  stateJac(2, 3) = tmp9;
    // 'computeStateJacNoGps:70' stateJac(10) = tmp9;
    stateJac[9] = tmp9;

    //  stateJac(2, 4) = tmp2;
    // 'computeStateJacNoGps:73' stateJac(11) = tmp2;
    stateJac[10] = tmp2;

    //  stateJac(2, 6) = tmp10;
    // 'computeStateJacNoGps:76' stateJac(12) = tmp10;
    stateJac[11] = rtb_UnitDelay_e;

    //  stateJac(2, 7) = tmp7;
    // 'computeStateJacNoGps:79' stateJac(13) = tmp7;
    stateJac[12] = rtb_XAxis2;

    //  stateJac(2, 8) = -tmp6;
    // 'computeStateJacNoGps:82' stateJac(14) = -tmp6;
    stateJac[13] = -tmp6;

    //  stateJac(3, 1) = tmp11;
    // 'computeStateJacNoGps:85' stateJac(15) = tmp11;
    stateJac[14] = tmp11;

    //  stateJac(3, 2) = tmp3;
    // 'computeStateJacNoGps:88' stateJac(16) = tmp3;
    stateJac[15] = rtb_Product1;

    //  stateJac(3, 3) = 1;
    // 'computeStateJacNoGps:91' stateJac(17) = 1;
    stateJac[16] = 1.0F;

    //  stateJac(3, 4) = tmp8;
    // 'computeStateJacNoGps:94' stateJac(18) = tmp8;
    stateJac[17] = tmp8;

    //  stateJac(3, 6) = -tmp7;
    // 'computeStateJacNoGps:97' stateJac(19) = -tmp7;
    stateJac[18] = -rtb_XAxis2;

    //  stateJac(3, 7) = tmp10;
    // 'computeStateJacNoGps:100' stateJac(20) = tmp10;
    stateJac[19] = rtb_UnitDelay_e;

    //  stateJac(3, 8) = tmp5;
    // 'computeStateJacNoGps:103' stateJac(21) = tmp5;
    stateJac[20] = rtb_XAxis1;

    //  stateJac(4, 1) = tmp9;
    // 'computeStateJacNoGps:106' stateJac(22) = tmp9;
    stateJac[21] = tmp9;

    //  stateJac(4, 2) = tmp11;
    // 'computeStateJacNoGps:109' stateJac(23) = tmp11;
    stateJac[22] = tmp11;

    //  stateJac(4, 3) = tmp1;
    // 'computeStateJacNoGps:112' stateJac(24) = tmp1;
    stateJac[23] = tmp1;

    //  stateJac(4, 4) = 1;
    // 'computeStateJacNoGps:115' stateJac(25) = 1;
    stateJac[24] = 1.0F;

    //  stateJac(4, 6) = tmp6;
    // 'computeStateJacNoGps:118' stateJac(26) = tmp6;
    stateJac[25] = tmp6;

    //  stateJac(4, 7) = -tmp5;
    // 'computeStateJacNoGps:121' stateJac(27) = -tmp5;
    stateJac[26] = -rtb_XAxis1;

    //  stateJac(4, 8) = tmp10;
    // 'computeStateJacNoGps:124' stateJac(28) = tmp10;
    stateJac[27] = rtb_UnitDelay_e;

    //  idx = 1;
    //
    //  % stateJac(1, 1) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(1, 2) = sampleTime_s*(bwx/2 - wx/2);
    //  stateJac(idx) = sampleTime_s*(bwx/2 - wx/2);
    //  idx = idx + 1;
    //
    //  % stateJac(1, 3) = sampleTime_s*(bwy/2 - wy/2);
    //  stateJac(idx) = sampleTime_s*(bwy/2 - wy/2);
    //  idx = idx + 1;
    //
    //  % stateJac(1, 4) = sampleTime_s*(bwz/2 - wz/2);
    //  stateJac(idx) = sampleTime_s*(bwz/2 - wz/2);
    //  idx = idx + 1;
    //
    //  % stateJac(1, 6) = (q1*sampleTime_s)/2;
    //  stateJac(idx) = (q1*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(1, 7) = (q2*sampleTime_s)/2;
    //  stateJac(idx) = (q2*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(1, 8) = (q3*sampleTime_s)/2;
    //  stateJac(idx) = (q3*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 1) = -(sampleTime_s*(bwx - wx))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwx - wx))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 2) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 3) = -(sampleTime_s*(bwz - wz))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwz - wz))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 4) = sampleTime_s*(bwy/2 - wy/2);
    //  stateJac(idx) = sampleTime_s*(bwy/2 - wy/2);
    //  idx = idx + 1;
    //
    //  % stateJac(2, 6) = -(q0*sampleTime_s)/2;
    //  stateJac(idx) = -(q0*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 7) = (q3*sampleTime_s)/2;
    //  stateJac(idx) = (q3*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(2, 8) = -(q2*sampleTime_s)/2;
    //  stateJac(idx) = -(q2*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 1) = -(sampleTime_s*(bwy - wy))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwy - wy))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 2) = sampleTime_s*(bwz/2 - wz/2);
    //  stateJac(idx) = sampleTime_s*(bwz/2 - wz/2);
    //  idx = idx + 1;
    //
    //  % stateJac(3, 3) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 4) = -(sampleTime_s*(bwx - wx))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwx - wx))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 6) = -(q3*sampleTime_s)/2;
    //  stateJac(idx) = -(q3*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 7) = -(q0*sampleTime_s)/2;
    //  stateJac(idx) = -(q0*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(3, 8) = (q1*sampleTime_s)/2;
    //  stateJac(idx) = (q1*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 1) = -(sampleTime_s*(bwz - wz))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwz - wz))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 2) = -(sampleTime_s*(bwy - wy))/2;
    //  stateJac(idx) = -(sampleTime_s*(bwy - wy))/2;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 3) = sampleTime_s*(bwx/2 - wx/2);
    //  stateJac(idx) = sampleTime_s*(bwx/2 - wx/2);
    //  idx = idx + 1;
    //
    //  % stateJac(4, 4) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 6) = (q2*sampleTime_s)/2;
    //  stateJac(idx) = (q2*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 7) = -(q1*sampleTime_s)/2;
    //  stateJac(idx) = -(q1*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(4, 8) = -(q0*sampleTime_s)/2;
    //  stateJac(idx) = -(q0*sampleTime_s)/2;
    //  idx = idx + 1;
    //
    //  % stateJac(5, 5) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(6, 6) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(7, 7) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(8, 8) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(9, 9) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(10, 10) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(11, 11) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(12, 12) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(13, 13) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(14, 14) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(15, 15) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(16, 16) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(17, 17) = 1;
    //  stateJac(idx) = 1;
    //  idx = idx + 1;
    //
    //  % stateJac(18, 18) = 1;
    //  stateJac(idx) = 1;
    //  covP = stateJac * covP * stateJac' + processNoiseQ;
    // 'ekfNoGps_function:53' covP = updateCovPNoGps(covP, stateJac, processNoiseQ); 
    updateCovPNoGps_vX01OT1j(stateEstimator_DW.covP, stateJac,
      rtu_processNoiseNoGpsQ);

    //  covP = propagateCovInplace(covP, stateJac,processNoiseQ);
    // 'ekfNoGps_function:56' accelMeasJac = computeAccelMeasJacNoGps(states, gEarth_mps2); 
    // COMPUTEACCELaccelMeasJacNOGPS Computes Meas Jacobian for accelerometer
    // measurements
    //
    // Inputs:
    // states:                EKF states
    //
    // Outputs:
    // accelaccelMeasJac:            3x18 Accel Meas Jacobian
    // Initialize the Meas jacobian to zero
    // 'computeAccelMeasJacNoGps:12' accelMeasJac = zeros(3, 18, 'single');
    std::memset(&measJac_0[0], 0, 54U * sizeof(real32_T));

    // Extract quat states
    // 'computeAccelMeasJacNoGps:15' q0 = states(1);
    // 'computeAccelMeasJacNoGps:16' q1 = states(2);
    // 'computeAccelMeasJacNoGps:17' q2 = states(3);
    // 'computeAccelMeasJacNoGps:18' q3 = states(4);
    // 'computeAccelMeasJacNoGps:20' tmp1 = gEarth_mps2*2*q2;
    rtb_Product1 = *rtu_gEarth_mps2 * 2.0F * rtb_states_a[2];

    // 'computeAccelMeasJacNoGps:21' tmp2 = -2*gEarth_mps2*q3;
    rtb_XAxis = -2.0F * *rtu_gEarth_mps2 * tmp17;

    // 'computeAccelMeasJacNoGps:22' tmp3 = gEarth_mps2*2*q0;
    rtb_XAxis1 = *rtu_gEarth_mps2 * 2.0F * rtb_states_a[0];

    // 'computeAccelMeasJacNoGps:23' tmp4 = -gEarth_mps2*2*q1;
    rtb_XAxis2 = -*rtu_gEarth_mps2 * 2.0F * rtb_states_a[1];

    // 'computeAccelMeasJacNoGps:24' tmp5 = 4*gEarth_mps2;
    rtb_UnitDelay_e = 4.0F * *rtu_gEarth_mps2;

    // 'computeAccelMeasJacNoGps:26' accelMeasJac(1, 1) = tmp1;
    measJac_0[0] = rtb_Product1;

    // 'computeAccelMeasJacNoGps:27' accelMeasJac(1, 2) = tmp2;
    measJac_0[3] = rtb_XAxis;

    // 'computeAccelMeasJacNoGps:28' accelMeasJac(1, 3) = tmp3;
    measJac_0[6] = rtb_XAxis1;

    // 'computeAccelMeasJacNoGps:29' accelMeasJac(1, 4) = tmp4;
    measJac_0[9] = rtb_XAxis2;

    // 'computeAccelMeasJacNoGps:30' accelMeasJac(1, 9) = 1;
    measJac_0[24] = 1.0F;

    // 'computeAccelMeasJacNoGps:32' accelMeasJac(2, 1) = tmp4;
    measJac_0[1] = rtb_XAxis2;

    // 'computeAccelMeasJacNoGps:33' accelMeasJac(2, 2) = -tmp3;
    measJac_0[4] = -rtb_XAxis1;

    // 'computeAccelMeasJacNoGps:34' accelMeasJac(2, 3) = tmp2;
    measJac_0[7] = rtb_XAxis;

    // 'computeAccelMeasJacNoGps:35' accelMeasJac(2, 4) = -tmp1;
    measJac_0[10] = -rtb_Product1;

    // 'computeAccelMeasJacNoGps:36' accelMeasJac(2, 10) = 1;
    measJac_0[28] = 1.0F;

    // 'computeAccelMeasJacNoGps:38' accelMeasJac(3, 2) = q1*tmp5;
    measJac_0[5] = rtb_states_a[1] * rtb_UnitDelay_e;

    // 'computeAccelMeasJacNoGps:39' accelMeasJac(3, 3) = q2*tmp5;
    measJac_0[8] = rtb_states_a[2] * rtb_UnitDelay_e;

    // 'computeAccelMeasJacNoGps:40' accelMeasJac(3, 11) = 1;
    measJac_0[32] = 1.0F;

    // Fuse Accel data if it is valid
    // 'ekfNoGps_function:59' [states, covP] = fuseAccelData(states, covP, bodyAccels_mps2, gEarth_mps2, accelMeasJac, measNoiseR); 
    // FUSEACCELDATA fuses Accelerometer data in the EKF when NO GPS
    //
    // Inputs:
    // states:            EKF states
    // covP:              State Covariance
    // bodyAccels_mps2:   Accelerometer measurements in body frame
    // gEarth_mps2:       Acceleration due to gravity
    // measJac:           Measurement Jacobian
    // measNoiseR:        Magnetometer measurement noise
    //
    // Outputs:
    // states:            Corrected states after fusing
    // covP:              Corrected covariance after fusing
    // Quaternions
    // 'fuseAccelData:18' q0 = states(1);
    // 'fuseAccelData:19' q1 = states(2);
    // 'fuseAccelData:20' q2 = states(3);
    // 'fuseAccelData:21' q3 = states(4);
    //  Direction Cosine Matrix (DCM) from NED coordinates to body cooridinates
    //  expressed using quaternions, using the current state estimate.
    // 'fuseAccelData:25' C_ned2b  = [ 1-2*(q2^2+q3^2)    2*(q1*q2+q3*q0)     2*(q1*q3-q2*q0); ... 
    // 'fuseAccelData:26'              2*(q1*q2-q3*q0)    1-2*(q1^2+q3^2)     2*(q2*q3+q1*q0); ... 
    // 'fuseAccelData:27'              2*(q1*q3+q2*q0)    2*(q2*q3-q1*q0)     1-2*(q1^2+q2^2)]; 
    // NED gravity in body frame
    // 'fuseAccelData:30' estGravityInBodyFrame = C_ned2b*[0; 0; -gEarth_mps2] + states(9:11); 
    // 'fuseAccelData:32' tmp1 = covP * measJac';
    for (i_0 = 0; i_0 < 18; i_0++) {
      f_y_tmp_0 = 0;
      for (i = 0; i < 3; i++) {
        b_tmp1_tmp = f_y_tmp_0 + i_0;
        c_tmp1_0[b_tmp1_tmp] = 0.0F;
        i_1 = 0;
        covP_tmp = 0;
        for (d_tmp1_tmp = 0; d_tmp1_tmp < 18; d_tmp1_tmp++) {
          c_tmp1_0[b_tmp1_tmp] += stateEstimator_DW.covP[i_1 + i_0] *
            measJac_0[covP_tmp + i];
          i_1 += 18;
          covP_tmp += 3;
        }

        f_y_tmp_0 += 18;
      }
    }

    // 'fuseAccelData:33' K = tmp1/(measJac * tmp1 + measNoiseR(4:6, 4:6));
    for (i_0 = 0; i_0 < 3; i_0++) {
      f_y_tmp_0 = 0;
      i = 0;
      i_1 = 0;
      for (covP_tmp = 0; covP_tmp < 3; covP_tmp++) {
        rtb_XAxis = 0.0F;
        d_tmp1_tmp = 0;
        for (b_tmp1_tmp = 0; b_tmp1_tmp < 18; b_tmp1_tmp++) {
          rtb_XAxis += measJac_0[d_tmp1_tmp + i_0] * c_tmp1_0[b_tmp1_tmp + i_1];
          d_tmp1_tmp += 3;
        }

        measJac_1[f_y_tmp_0 + i_0] = rtu_measNoiseNoGpsR[(i + i_0) + 27] +
          rtb_XAxis;
        f_y_tmp_0 += 3;
        i += 8;
        i_1 += 18;
      }
    }

    mrdiv_2CyhSJN4(c_tmp1_0, measJac_1, K_0);

    // 'fuseAccelData:34' states = states + K*(bodyAccels_mps2 - estGravityInBodyFrame); 
    rtb_XAxis1 = bodyAccelsOut_mps2[0] - ((rtb_states_a[1] * tmp17 -
      rtb_states_a[0] * rtb_states_a[2]) * 2.0F * -*rtu_gEarth_mps2 + d_tmp1_0[8]);
    rtb_XAxis2 = bodyAccelsOut_mps2[1] - ((rtb_states_a[2] * tmp17 +
      rtb_states_a[0] * rtb_states_a[1]) * 2.0F * -*rtu_gEarth_mps2 + d_tmp1_0[9]);
    rtb_UnitDelay_e = bodyAccelsOut_mps2[2] - ((1.0F - (rtb_states_a[1] *
      rtb_states_a[1] + rtb_states_a[2] * rtb_states_a[2]) * 2.0F) *
      -*rtu_gEarth_mps2 + d_tmp1_0[10]);

    // 'fuseAccelData:35' covP = covP - K * measJac*covP;
    for (i_0 = 0; i_0 < 18; i_0++) {
      rtb_XAxis = K_0[i_0 + 18];
      tmp1 = K_0[i_0 + 36];
      d_tmp1_0[i_0] = ((rtb_XAxis * rtb_XAxis2 + K_0[i_0] * rtb_XAxis1) + tmp1 *
                       rtb_UnitDelay_e) + rtb_states_a[i_0];
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        i = 18 * f_y_tmp_0 + i_0;
        K_2[i] = 0.0F;
        K_2[i] += measJac_0[3 * f_y_tmp_0] * K_0[i_0];
        K_2[i] += measJac_0[3 * f_y_tmp_0 + 1] * rtb_XAxis;
        K_2[i] += measJac_0[3 * f_y_tmp_0 + 2] * tmp1;
      }

      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        rtb_XAxis = 0.0F;
        for (i = 0; i < 18; i++) {
          rtb_XAxis += K_2[18 * i + i_0] * stateEstimator_DW.covP[18 * f_y_tmp_0
            + i];
        }

        i = 18 * f_y_tmp_0 + i_0;
        covP_0[i] = stateEstimator_DW.covP[i] - rtb_XAxis;
      }
    }

    std::memcpy(&stateEstimator_DW.covP[0], &covP_0[0], 324U * sizeof(real32_T));

    // 'fuseAccelData:36' covP = (covP + covP')/2;
    i_0 = 0;
    for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
      i = 0;
      for (i_1 = 0; i_1 < 18; i_1++) {
        covP_tmp = i_1 + i_0;
        covP_0[covP_tmp] = (stateEstimator_DW.covP[i + f_y_tmp_0] +
                            stateEstimator_DW.covP[covP_tmp]) / 2.0F;
        i += 18;
      }

      i_0 += 18;
    }

    std::memcpy(&stateEstimator_DW.covP[0], &covP_0[0], 324U * sizeof(real32_T));
    std::memcpy(&rtb_states_a[0], &d_tmp1_0[0], 18U * sizeof(real32_T));

    // Normalize the quaternion
    // 'ekfNoGps_function:61' states(1:4) = states(1:4)/norm(states(1:4));
    rtb_XAxis = norm_7MzYkgry(&d_tmp1_0[0]);
    rtb_states_a[0] = d_tmp1_0[0] / rtb_XAxis;
    rtb_states_a[1] = d_tmp1_0[1] / rtb_XAxis;
    rtb_states_a[2] = d_tmp1_0[2] / rtb_XAxis;
    rtb_states_a[3] = d_tmp1_0[3] / rtb_XAxis;

    // Fuse Mag data if it is valid
    // 'ekfNoGps_function:64' if(isMagValid)
    if (isMagValid) {
      // 'ekfNoGps_function:65' measJac = computeMagMeasJacNoGps(states);
      computeMagMeasJacNoGps_mcxOPfaA(rtb_states_a, measJac_0);

      // 'ekfNoGps_function:66' [states, covP] = fuseMagDataNoGps(states, covP, normMagVec_nd, measJac, measNoiseR); 
      // FUSEMAGDATANOGPS fuses Magnetometer data in the EKF NO GPS
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // bodyMagUnitVec:    Normalized magnetometer measurement in body frame
      // measJac:           Measurement Jacobian
      // measNoiseR:        Magnetometer measurement noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // Quaternions
      // 'fuseMagDataNoGps:16' q0 = states(1);
      // 'fuseMagDataNoGps:17' q1 = states(2);
      // 'fuseMagDataNoGps:18' q2 = states(3);
      // 'fuseMagDataNoGps:19' q3 = states(4);
      // NED Mag Norm and biases
      // 'fuseMagDataNoGps:22' nedMagUnitVec = states(12:14);
      //  Direction Cosine Matrix (DCM) from NED coordinates to body cooridinates 
      //  expressed using quaternions, using the current state estimate.
      // 'fuseMagDataNoGps:26' C_ned2b  = [ 1-2*(q2^2+q3^2)    2*(q1*q2+q3*q0)     2*(q1*q3-q2*q0); ... 
      // 'fuseMagDataNoGps:27'              2*(q1*q2-q3*q0)    1-2*(q1^2+q3^2)     2*(q2*q3+q1*q0); ... 
      // 'fuseMagDataNoGps:28'              2*(q1*q3+q2*q0)    2*(q2*q3-q1*q0)     1-2*(q1^2+q2^2)]; 
      // Rotate propogated mag states and add bias to estimate measurements
      // 'fuseMagDataNoGps:31' estBodyMagUnitVec = C_ned2b*nedMagUnitVec + states(15:17); 
      // 'fuseMagDataNoGps:33' tmp1 = covP * measJac';
      for (i_0 = 0; i_0 < 18; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 3; i++) {
          b_tmp1_tmp = f_y_tmp_0 + i_0;
          c_tmp1_0[b_tmp1_tmp] = 0.0F;
          i_1 = 0;
          covP_tmp = 0;
          for (d_tmp1_tmp = 0; d_tmp1_tmp < 18; d_tmp1_tmp++) {
            c_tmp1_0[b_tmp1_tmp] += stateEstimator_DW.covP[i_1 + i_0] *
              measJac_0[covP_tmp + i];
            i_1 += 18;
            covP_tmp += 3;
          }

          f_y_tmp_0 += 18;
        }
      }

      // 'fuseMagDataNoGps:34' K = tmp1/(measJac * tmp1 + measNoiseR(1:3, 1:3)); 
      for (i_0 = 0; i_0 < 3; i_0++) {
        f_y_tmp_0 = 0;
        i = 0;
        i_1 = 0;
        for (covP_tmp = 0; covP_tmp < 3; covP_tmp++) {
          rtb_XAxis = 0.0F;
          d_tmp1_tmp = 0;
          for (b_tmp1_tmp = 0; b_tmp1_tmp < 18; b_tmp1_tmp++) {
            rtb_XAxis += measJac_0[d_tmp1_tmp + i_0] * c_tmp1_0[b_tmp1_tmp + i_1];
            d_tmp1_tmp += 3;
          }

          measJac_1[f_y_tmp_0 + i_0] = rtu_measNoiseNoGpsR[i + i_0] + rtb_XAxis;
          f_y_tmp_0 += 3;
          i += 8;
          i_1 += 18;
        }
      }

      mrdiv_2CyhSJN4(c_tmp1_0, measJac_1, K_0);

      // 'fuseMagDataNoGps:35' states = states + K*(bodyMagUnitVec - estBodyMagUnitVec); 
      rtb_XAxis2 = rtb_states_a[3] * rtb_states_a[3];
      tmp6 = rtb_states_a[1] * rtb_states_a[2];
      tmp8 = rtb_states_a[0] * rtb_states_a[3];
      rtb_XAxis = rtb_states_a[1] * rtb_states_a[3];
      rtb_UnitDelay_e = rtb_states_a[0] * rtb_states_a[2];
      tmp1 = rtb_states_a[2] * rtb_states_a[2];
      rtb_XAxis1 = stateEstimator_DW.normMagVecOut_nd[0] - ((((1.0F - (tmp1 +
        rtb_XAxis2) * 2.0F) * rtb_states_a[11] + (tmp6 + tmp8) * 2.0F *
        rtb_states_a[12]) + (rtb_XAxis - rtb_UnitDelay_e) * 2.0F * rtb_states_a
        [13]) + rtb_states_a[14]);
      tmp2 = rtb_states_a[2] * rtb_states_a[3];
      rtb_Product1 = rtb_states_a[0] * rtb_states_a[1];
      tmp9 = rtb_states_a[1] * rtb_states_a[1];
      rtb_XAxis2 = stateEstimator_DW.normMagVecOut_nd[1] - ((((1.0F - (tmp9 +
        rtb_XAxis2) * 2.0F) * rtb_states_a[12] + (tmp6 - tmp8) * 2.0F *
        rtb_states_a[11]) + (tmp2 + rtb_Product1) * 2.0F * rtb_states_a[13]) +
        rtb_states_a[15]);
      rtb_UnitDelay_e = stateEstimator_DW.normMagVecOut_nd[2] - ((((rtb_XAxis +
        rtb_UnitDelay_e) * 2.0F * rtb_states_a[11] + (tmp2 - rtb_Product1) *
        2.0F * rtb_states_a[12]) + (1.0F - (tmp9 + tmp1) * 2.0F) * rtb_states_a
        [13]) + rtb_states_a[16]);

      // 'fuseMagDataNoGps:36' covP = covP - K * measJac*covP;
      for (i_0 = 0; i_0 < 18; i_0++) {
        rtb_XAxis = K_0[i_0 + 18];
        tmp1 = K_0[i_0 + 36];
        d_tmp1_0[i_0] = ((rtb_XAxis * rtb_XAxis2 + K_0[i_0] * rtb_XAxis1) + tmp1
                         * rtb_UnitDelay_e) + rtb_states_a[i_0];
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
          i = 18 * f_y_tmp_0 + i_0;
          K_2[i] = 0.0F;
          K_2[i] += measJac_0[3 * f_y_tmp_0] * K_0[i_0];
          K_2[i] += measJac_0[3 * f_y_tmp_0 + 1] * rtb_XAxis;
          K_2[i] += measJac_0[3 * f_y_tmp_0 + 2] * tmp1;
        }

        for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
          rtb_XAxis = 0.0F;
          for (i = 0; i < 18; i++) {
            rtb_XAxis += K_2[18 * i + i_0] * stateEstimator_DW.covP[18 *
              f_y_tmp_0 + i];
          }

          covP_tmp = 18 * f_y_tmp_0 + i_0;
          covP_0[covP_tmp] = stateEstimator_DW.covP[covP_tmp] - rtb_XAxis;
        }
      }

      // 'fuseMagDataNoGps:37' covP = (covP + covP')/2;
      i = 0;
      for (i_1 = 0; i_1 < 18; i_1++) {
        i_0 = 0;
        for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
          covP_tmp = f_y_tmp_0 + i;
          stateEstimator_DW.covP[covP_tmp] = (covP_0[i_0 + i_1] +
            covP_0[covP_tmp]) / 2.0F;
          i_0 += 18;
        }

        rtb_states_a[i_1] = d_tmp1_0[i_1];
        i += 18;
      }

      // Normalize the quaternion
      // 'ekfNoGps_function:68' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&d_tmp1_0[0]);
      rtb_states_a[0] = d_tmp1_0[0] / rtb_Product1;
      rtb_states_a[1] = d_tmp1_0[1] / rtb_Product1;
      rtb_states_a[2] = d_tmp1_0[2] / rtb_Product1;
      rtb_states_a[3] = d_tmp1_0[3] / rtb_Product1;
    }

    // Fuse Baro data if it is valid
    // 'ekfNoGps_function:72' if(isBaroValid)
    if (isBaroValid) {
      // 'ekfNoGps_function:73' measJac = zeros(1, 18, 'single');
      // 'ekfNoGps_function:74' measJac(1, 5) = -1;
      // 'ekfNoGps_function:75' measJac(1, 18) = 1;
      // 'ekfNoGps_function:76' [states, covP] = fuseBaroDataNoGps(states, covP, baroAlt_m, measJac, measNoiseR); 
      // FUSEBARODATANOGPS Fuses BARO data in EKF No GPS
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // baroAltm:          Baro alt data
      // measJac:           Measurement Jacobian
      // measNoiseR:        Baro alt meas noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // 'fuseBaroDataNoGps:15' tmp1 = covP * measJac';
      // 'fuseBaroDataNoGps:16' K = tmp1/(measJac * tmp1 + measNoiseR(7, 7));
      rtb_XAxis1 = 0.0F;
      for (i_0 = 0; i_0 < 18; i_0++) {
        d_tmp1_0[i_0] = 0.0F;
        f_y_tmp_0 = 0;
        for (i = 0; i < 18; i++) {
          d_tmp1_0[i_0] += stateEstimator_DW.covP[f_y_tmp_0 + i_0] *
            static_cast<real32_T>(b_measJac_0[i]);
          f_y_tmp_0 += 18;
        }

        f_y_tmp_0 = b_measJac_0[i_0];
        rtb_XAxis1 += static_cast<real32_T>(f_y_tmp_0) * d_tmp1_0[i_0];
        g_y_tmp[i_0] = static_cast<int8_T>(f_y_tmp_0);
      }

      rtb_Product1 = rtb_XAxis1 + rtu_measNoiseNoGpsR[54];
      for (i_0 = 0; i_0 < 18; i_0++) {
        d_tmp1_0[i_0] /= rtb_Product1;
      }

      // 'fuseBaroDataNoGps:17' states = states + K*(baroAlt_m  + states(5) - states(18)); 
      rtb_Product1 = (baroAltOut_m + rtb_states_a[4]) - rtb_states_a[17];

      // 'fuseBaroDataNoGps:18' covP = covP - K*measJac*covP;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        for (i = 0; i < 18; i++) {
          K_2[i + i_0] = d_tmp1_0[i] * static_cast<real32_T>(g_y_tmp[f_y_tmp_0]);
        }

        rtb_states_a[f_y_tmp_0] += d_tmp1_0[f_y_tmp_0] * rtb_Product1;
        i_0 += 18;
      }

      for (i_0 = 0; i_0 < 18; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 18; i++) {
          rtb_XAxis = 0.0F;
          i_1 = 0;
          for (covP_tmp = 0; covP_tmp < 18; covP_tmp++) {
            rtb_XAxis += K_2[i_1 + i_0] * stateEstimator_DW.covP[covP_tmp +
              f_y_tmp_0];
            i_1 += 18;
          }

          covP_tmp = f_y_tmp_0 + i_0;
          covP_0[covP_tmp] = stateEstimator_DW.covP[covP_tmp] - rtb_XAxis;
          f_y_tmp_0 += 18;
        }
      }

      // 'fuseBaroDataNoGps:19' covP = (covP + covP')/2;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        i = 0;
        for (i_1 = 0; i_1 < 18; i_1++) {
          covP_tmp = i_1 + i_0;
          stateEstimator_DW.covP[covP_tmp] = (covP_0[i + f_y_tmp_0] +
            covP_0[covP_tmp]) / 2.0F;
          i += 18;
        }

        i_0 += 18;
      }

      // Normalize the quaternion
      // 'ekfNoGps_function:78' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&rtb_states_a[0]);
      rtb_states_a[0] /= rtb_Product1;
      rtb_states_a[1] /= rtb_Product1;
      rtb_states_a[2] /= rtb_Product1;
      rtb_states_a[3] /= rtb_Product1;
    }

    // Fuse Lidar data if it is valid
    // 'ekfNoGps_function:82' if(isLidarValid)
    if (rtb_AND1) {
      // 'ekfNoGps_function:83' measJac = zeros(1, 18, 'single');
      // 'ekfNoGps_function:84' measJac(1, 5) = -1;
      // 'ekfNoGps_function:85' [states, covP] = fuseLidarDataNoGps(states, covP, lidarAgl_m, measJac, measNoiseR); 
      // FUSELIDARDATANOGPS Fuses LIDAR data in EKF No GPS
      //
      // Inputs:
      // states:            EKF states
      // covP:              State Covariance
      // lidarAgl_m:        Lidar agl data
      // measJac:           Measurement Jacobian
      // measNoiseR:        Baro alt meas noise
      //
      // Outputs:
      // states:            Corrected states after fusing
      // covP:              Corrected covariance after fusing
      // 'fuseLidarDataNoGps:15' tmp1 = covP*measJac';
      // 'fuseLidarDataNoGps:16' K = tmp1/(measJac * tmp1 + measNoiseR(8, 8));
      rtb_XAxis1 = 0.0F;
      for (i_0 = 0; i_0 < 18; i_0++) {
        d_tmp1_0[i_0] = 0.0F;
        f_y_tmp_0 = 0;
        for (i = 0; i < 18; i++) {
          d_tmp1_0[i_0] += stateEstimator_DW.covP[f_y_tmp_0 + i_0] *
            static_cast<real32_T>(c_measJac_0[i]);
          f_y_tmp_0 += 18;
        }

        f_y_tmp_0 = c_measJac_0[i_0];
        rtb_XAxis1 += static_cast<real32_T>(f_y_tmp_0) * d_tmp1_0[i_0];
        g_y_tmp[i_0] = static_cast<int8_T>(f_y_tmp_0);
      }

      rtb_Product1 = rtb_XAxis1 + rtu_measNoiseNoGpsR[63];
      for (i_0 = 0; i_0 < 18; i_0++) {
        d_tmp1_0[i_0] /= rtb_Product1;
      }

      // 'fuseLidarDataNoGps:17' states = states + K*(lidarAgl_m  + states(5));
      rtb_Product1 = rtb_Product2 + rtb_states_a[4];

      // 'fuseLidarDataNoGps:18' covP = covP - K*measJac*covP;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        for (i = 0; i < 18; i++) {
          K_2[i + i_0] = d_tmp1_0[i] * static_cast<real32_T>(g_y_tmp[f_y_tmp_0]);
        }

        rtb_states_a[f_y_tmp_0] += d_tmp1_0[f_y_tmp_0] * rtb_Product1;
        i_0 += 18;
      }

      for (i_0 = 0; i_0 < 18; i_0++) {
        f_y_tmp_0 = 0;
        for (i = 0; i < 18; i++) {
          rtb_XAxis = 0.0F;
          i_1 = 0;
          for (covP_tmp = 0; covP_tmp < 18; covP_tmp++) {
            rtb_XAxis += K_2[i_1 + i_0] * stateEstimator_DW.covP[covP_tmp +
              f_y_tmp_0];
            i_1 += 18;
          }

          covP_tmp = f_y_tmp_0 + i_0;
          covP_0[covP_tmp] = stateEstimator_DW.covP[covP_tmp] - rtb_XAxis;
          f_y_tmp_0 += 18;
        }
      }

      // 'fuseLidarDataNoGps:19' covP = (covP + covP')/2;
      i_0 = 0;
      for (f_y_tmp_0 = 0; f_y_tmp_0 < 18; f_y_tmp_0++) {
        i = 0;
        for (i_1 = 0; i_1 < 18; i_1++) {
          covP_tmp = i_1 + i_0;
          stateEstimator_DW.covP[covP_tmp] = (covP_0[i + f_y_tmp_0] +
            covP_0[covP_tmp]) / 2.0F;
          i += 18;
        }

        i_0 += 18;
      }

      // Normalize the quaternion
      // 'ekfNoGps_function:87' states(1:4) = states(1:4)/norm(states(1:4));
      rtb_Product1 = norm_7MzYkgry(&rtb_states_a[0]);
      rtb_states_a[0] /= rtb_Product1;
      rtb_states_a[1] /= rtb_Product1;
      rtb_states_a[2] /= rtb_Product1;
      rtb_states_a[3] /= rtb_Product1;
    }

    // Compute Body To NED DCM
    // 'ekfNoGps_function:91' dcmBodyToNed = quatToDcm_function(states(1:4));
    // Quaternions
    // 'quatToDcm_function:3' q0 = quat(1);
    // 'quatToDcm_function:4' q1 = quat(2);
    // 'quatToDcm_function:5' q2 = quat(3);
    // 'quatToDcm_function:6' q3 = quat(4);
    //  Direction Cosine Matrix (DCM) from body cooridinates to NED coordinates
    //  expressed using quaternions.
    // 'quatToDcm_function:10' dcmBodyToNed = [1-2*(q2^2+q3^2), 2*(q1*q2-q3*q0), 2*(q1*q3+q2*q0); 
    // 'quatToDcm_function:11'     2*(q1*q2+q3*q0), 1-2*(q1^2+q3^2), 2*(q2*q3-q1*q0); 
    // 'quatToDcm_function:12'     2*(q1*q3-q2*q0), 2*(q2*q3+q1*q0), 1-2*(q1^2+q2^2)]; 
    tmp1 = rtb_states_a[3] * rtb_states_a[3];
    rtb_Product2 = rtb_states_a[2] * rtb_states_a[2];
    rtb_XAxis1 = 1.0F - (rtb_Product2 + tmp1) * 2.0F;
    rtb_UnitDelay_e = rtb_states_a[1] * rtb_states_a[2];
    tmp2 = rtb_states_a[0] * rtb_states_a[3];
    rtb_XAxis = (rtb_UnitDelay_e - tmp2) * 2.0F;
    baroAltOut_m = rtb_states_a[1] * rtb_states_a[3];
    rtb_Product1 = rtb_states_a[0] * rtb_states_a[2];
    rtb_XAxis2 = (baroAltOut_m + rtb_Product1) * 2.0F;
    rtb_UnitDelay_e = (rtb_UnitDelay_e + tmp2) * 2.0F;
    tmp6 = rtb_states_a[1] * rtb_states_a[1];
    tmp1 = 1.0F - (tmp6 + tmp1) * 2.0F;
    tmp8 = rtb_states_a[2] * rtb_states_a[3];
    tmp9 = rtb_states_a[0] * rtb_states_a[1];
    tmp2 = (tmp8 - tmp9) * 2.0F;
    rtb_Product1 = (baroAltOut_m - rtb_Product1) * 2.0F;
    tmp8 = (tmp8 + tmp9) * 2.0F;
    tmp6 = 1.0F - (tmp6 + rtb_Product2) * 2.0F;
  }

  // End of MATLAB Function: '<S1>/EKF NO GPS'

  // MATLAB Function: '<S1>/State Selector'
  // MATLAB Function 'EKF/State Selector': '<S13>:1'
  // '<S13>:1:4' if(smMode == enumStateEstimateMode.INITIALIZE || ...
  // '<S13>:1:5'    smMode == enumStateEstimateMode.RUN)
  if ((static_cast<boolean_T>((mode != enumStateEstimateMode::INITIALIZE) &
        (mode != enumStateEstimateMode::RUN))) && (mode != enumStateEstimateMode::
       RUN)) {
    // '<S13>:1:11' if(smMode ~= enumStateEstimateMode.RUN)
    // '<S13>:1:12' states = ekfStates;
    // '<S13>:1:13' states(1:4) = noGpsEkfStates(1:4);
    rty_states[0] = rtb_states_a[0];
    rty_states[1] = rtb_states_a[1];
    rty_states[2] = rtb_states_a[2];
    rty_states[3] = rtb_states_a[3];

    // '<S13>:1:14' states(7) = noGpsEkfStates(5);
    rty_states[6] = rtb_states_a[4];

    // '<S13>:1:15' states(11:23) = noGpsEkfStates(6:18);
    for (i = 0; i < 13; i++) {
      rty_states[i + 10] = rtb_states_a[i + 5];
    }

    // '<S13>:1:16' dcmBodyToNed = noGpsDcmBodyToNed;
    rtb_dcmBodyToNed_k_idx_0 = rtb_XAxis1;
    rtb_dcmBodyToNed_k_idx_1 = rtb_UnitDelay_e;
    rtb_dcmBodyToNed_k_idx_2 = rtb_Product1;
    rtb_dcmBodyToNed_k_idx_3 = rtb_XAxis;
    rtb_dcmBodyToNed_k_idx_4 = tmp1;
    rtb_dcmBodyToNed_k_idx_5 = tmp8;
    rtb_dcmBodyToNed_k_idx_6 = rtb_XAxis2;
    rtb_dcmBodyToNed_k_idx_7 = tmp2;
    rtb_dcmBodyToNed_k_idx_8 = tmp6;
  } else {
    // '<S13>:1:6' states = ekfStates;
    // '<S13>:1:7' dcmBodyToNed = ekfDcmBodyToNed;
    // '<S13>:1:20' states = ekfStates;
    // '<S13>:1:21' dcmBodyToNed = ekfDcmBodyToNed;
  }

  // End of MATLAB Function: '<S1>/State Selector'

  // Sqrt: '<S20>/sqrt' incorporates:
  //   Product: '<S21>/Product'
  //   Product: '<S21>/Product1'
  //   Product: '<S21>/Product2'
  //   Product: '<S21>/Product3'
  //   Sum: '<S21>/Sum'

  rtb_XAxis = std::sqrt(((rty_states[0] * rty_states[0] + rty_states[1] *
    rty_states[1]) + rty_states[2] * rty_states[2]) + rty_states[3] *
                        rty_states[3]);

  // Product: '<S15>/Product'
  rtb_XAxis1 = rty_states[0] / rtb_XAxis;

  // Product: '<S15>/Product1'
  rtb_Product1 = rty_states[1] / rtb_XAxis;

  // Product: '<S15>/Product2'
  rtb_Product2 = rty_states[2] / rtb_XAxis;

  // Product: '<S15>/Product3'
  rtb_XAxis = rty_states[3] / rtb_XAxis;

  // Fcn: '<S2>/fcn2' incorporates:
  //   Fcn: '<S2>/fcn5'

  baroAltOut_m = rtb_XAxis1 * rtb_XAxis1;
  rtb_XAxis2 = rtb_Product1 * rtb_Product1;
  rtb_UnitDelay_e = rtb_Product2 * rtb_Product2;
  tmp2 = rtb_XAxis * rtb_XAxis;

  // Trigonometry: '<S14>/Trigonometric Function1' incorporates:
  //   Fcn: '<S2>/fcn1'
  //   Fcn: '<S2>/fcn2'

  tmp1 = std::atan2((rtb_Product1 * rtb_Product2 + rtb_XAxis1 * rtb_XAxis) *
                    2.0F, ((baroAltOut_m + rtb_XAxis2) - rtb_UnitDelay_e) - tmp2);

  // Trigonometry: '<S14>/Trigonometric Function3' incorporates:
  //   Fcn: '<S2>/fcn4'
  //   Fcn: '<S2>/fcn5'

  tmp2 = std::atan2((rtb_Product2 * rtb_XAxis + rtb_XAxis1 * rtb_Product1) *
                    2.0F, ((baroAltOut_m - rtb_XAxis2) - rtb_UnitDelay_e) + tmp2);

  // Fcn: '<S2>/fcn3'
  rtb_Product2 = (rtb_Product1 * rtb_XAxis - rtb_XAxis1 * rtb_Product2) * -2.0F;

  // If: '<S16>/If' incorporates:
  //   Constant: '<S17>/Constant'
  //   Constant: '<S18>/Constant'

  if (rtb_Product2 > 1.0F) {
    // Outputs for IfAction SubSystem: '<S16>/If Action Subsystem' incorporates:
    //   ActionPort: '<S17>/Action Port'

    rtb_Product2 = 1.0F;

    // End of Outputs for SubSystem: '<S16>/If Action Subsystem'
  } else if (rtb_Product2 < -1.0F) {
    // Outputs for IfAction SubSystem: '<S16>/If Action Subsystem1' incorporates:
    //   ActionPort: '<S18>/Action Port'

    rtb_Product2 = 1.0F;

    // End of Outputs for SubSystem: '<S16>/If Action Subsystem1'
  }

  // End of If: '<S16>/If'

  // Trigonometry: '<S14>/trigFcn'
  rtb_XAxis1 = std::asin(rtb_Product2);

  // MATLAB Function: '<Root>/eulToDcm' incorporates:
  //   Gain: '<Root>/ZeroOutRollAndPitch'

  //  precalculate trignometric values
  // MATLAB Function 'eulToDcm': '<S6>:1'
  // '<S6>:1:4' s_phi = sin(eul_rad(1));
  // '<S6>:1:5' c_phi = cos(eul_rad(1));
  // '<S6>:1:7' s_theta = sin(eul_rad(2));
  // '<S6>:1:8' c_theta = cos(eul_rad(2));
  // '<S6>:1:10' s_psi = sin(eul_rad(3));
  rtb_XAxis2 = std::sin(tmp1);

  // '<S6>:1:11' c_psi = cos(eul_rad(3));
  rtb_UnitDelay_e = std::cos(tmp1);

  // '<S6>:1:13' dcmFromNed = [c_psi*c_theta, c_theta*s_psi, -s_theta;
  // '<S6>:1:14'     c_psi*s_phi*s_theta - c_phi*s_psi, c_phi*c_psi + s_phi*s_psi*s_theta, c_theta*s_phi; 
  // '<S6>:1:15'     s_phi*s_psi + c_phi*c_psi*s_theta, c_phi*s_psi*s_theta - c_psi*s_phi, c_phi*c_theta]; 
  rty_dcmNedToFep[0] = rtb_UnitDelay_e;
  rty_dcmNedToFep[3] = rtb_XAxis2;
  rty_dcmNedToFep[6] = -0.0F;
  rty_dcmNedToFep[1] = 0.0F - rtb_XAxis2;
  rty_dcmNedToFep[4] = rtb_UnitDelay_e;
  rty_dcmNedToFep[7] = 0.0F;
  rty_dcmNedToFep[2] = 0.0F;
  rty_dcmNedToFep[5] = 0.0F;
  rty_dcmNedToFep[8] = 1.0F;

  // SignalConversion generated from: '<Root>/eulAng_rad'
  rty_eulAng_rad[0] = tmp2;
  rty_eulAng_rad[1] = rtb_XAxis1;
  rty_eulAng_rad[2] = tmp1;

  // Math: '<Root>/Transpose'
  rty_dcmNedToBody[0] = rtb_dcmBodyToNed_k_idx_0;
  rty_dcmNedToBody[1] = rtb_dcmBodyToNed_k_idx_3;
  rty_dcmNedToBody[2] = rtb_dcmBodyToNed_k_idx_6;

  // Sum: '<Root>/Sum' incorporates:
  //   Product: '<S4>/Product'

  rty_bodyAccels_mps2[0] = stateEstimator_DW.Product[0] - rty_states[13];

  // Math: '<Root>/Transpose'
  rty_dcmNedToBody[3] = rtb_dcmBodyToNed_k_idx_1;
  rty_dcmNedToBody[4] = rtb_dcmBodyToNed_k_idx_4;
  rty_dcmNedToBody[5] = rtb_dcmBodyToNed_k_idx_7;

  // Sum: '<Root>/Sum' incorporates:
  //   Product: '<S4>/Product'

  rty_bodyAccels_mps2[1] = stateEstimator_DW.Product[1] - rty_states[14];

  // Math: '<Root>/Transpose'
  rty_dcmNedToBody[6] = rtb_dcmBodyToNed_k_idx_2;
  rty_dcmNedToBody[7] = rtb_dcmBodyToNed_k_idx_5;
  rty_dcmNedToBody[8] = rtb_dcmBodyToNed_k_idx_8;

  // Sum: '<Root>/Sum' incorporates:
  //   Product: '<S4>/Product'

  rty_bodyAccels_mps2[2] = stateEstimator_DW.Product[2] - rty_states[15];

  // BusCreator: '<Root>/Bus Creator'
  rty_stateEstimatorDebug->stateEstInitPct = stateEstimator_DW.stateEstInitPct;
  rty_stateEstimatorDebug->smMode = mode;

  // Update for DiscreteTransferFcn: '<S23>/X Axis'
  stateEstimator_DW.XAxis_states[1] = stateEstimator_DW.XAxis_states[0];
  stateEstimator_DW.XAxis_states[0] = stateEstimator_DW.XAxis_tmp;

  // Update for DiscreteTransferFcn: '<S23>/X Axis1'
  stateEstimator_DW.XAxis1_states[1] = stateEstimator_DW.XAxis1_states[0];
  stateEstimator_DW.XAxis1_states[0] = stateEstimator_DW.XAxis1_tmp;

  // Update for DiscreteTransferFcn: '<S23>/X Axis2'
  stateEstimator_DW.XAxis2_states[1] = stateEstimator_DW.XAxis2_states[0];
  stateEstimator_DW.XAxis2_states[0] = stateEstimator_DW.XAxis2_tmp;

  // Update for DiscreteTransferFcn: '<S24>/X Axis'
  stateEstimator_DW.XAxis_states_b[1] = stateEstimator_DW.XAxis_states_b[0];
  stateEstimator_DW.XAxis_states_b[0] = stateEstimator_DW.XAxis_tmp_k;

  // Update for DiscreteTransferFcn: '<S24>/X Axis1'
  stateEstimator_DW.XAxis1_states_m[1] = stateEstimator_DW.XAxis1_states_m[0];
  stateEstimator_DW.XAxis1_states_m[0] = stateEstimator_DW.XAxis1_tmp_m;

  // Update for DiscreteTransferFcn: '<S24>/X Axis2'
  stateEstimator_DW.XAxis2_states_h[1] = stateEstimator_DW.XAxis2_states_h[0];
  stateEstimator_DW.XAxis2_states_h[0] = stateEstimator_DW.XAxis2_tmp_g;

  // Update for Delay: '<S1>/Delay'
  stateEstimator_DW.icLoad = false;

  // Update for UnitDelay: '<Root>/Unit Delay'
  std::memcpy(&stateEstimator_DW.UnitDelay_DSTATE[0], &rty_states[0], 23U *
              sizeof(real32_T));

  // Update for Delay: '<S1>/Delay' incorporates:
  //   UnitDelay: '<Root>/Unit Delay'

  std::memcpy(&stateEstimator_DW.Delay_DSTATE[0], &rty_states[0], 23U * sizeof
              (real32_T));

  // Update for Delay: '<S1>/Delay2'
  stateEstimator_DW.icLoad_j = false;
  stateEstimator_DW.Delay2_DSTATE[0] = rtb_dcmBodyToNed_k_idx_0;
  stateEstimator_DW.Delay2_DSTATE[1] = rtb_dcmBodyToNed_k_idx_1;
  stateEstimator_DW.Delay2_DSTATE[2] = rtb_dcmBodyToNed_k_idx_2;
  stateEstimator_DW.Delay2_DSTATE[3] = rtb_dcmBodyToNed_k_idx_3;
  stateEstimator_DW.Delay2_DSTATE[4] = rtb_dcmBodyToNed_k_idx_4;
  stateEstimator_DW.Delay2_DSTATE[5] = rtb_dcmBodyToNed_k_idx_5;
  stateEstimator_DW.Delay2_DSTATE[6] = rtb_dcmBodyToNed_k_idx_6;
  stateEstimator_DW.Delay2_DSTATE[7] = rtb_dcmBodyToNed_k_idx_7;
  stateEstimator_DW.Delay2_DSTATE[8] = rtb_dcmBodyToNed_k_idx_8;

  // Switch: '<S7>/Switch' incorporates:
  //   RelationalOperator: '<S26>/FixPt Relational Operator'
  //   RelationalOperator: '<S28>/Compare'
  //   UnitDelay: '<S26>/Delay Input1'
  //
  //  Block description for '<S26>/Delay Input1':
  //
  //   Store in Global RAM

  if (static_cast<int32_T>(rtb_Compare) > static_cast<int32_T>
      (stateEstimator_DW.DelayInput1_DSTATE)) {
    // Update for UnitDelay: '<S7>/Unit Delay'
    stateEstimator_DW.UnitDelay_DSTATE_c = rtb_UnitDelay[6];
  }

  // End of Switch: '<S7>/Switch'

  // Update for UnitDelay: '<Root>/Unit Delay1'
  stateEstimator_DW.UnitDelay1_DSTATE[0] = tmp2;
  stateEstimator_DW.UnitDelay1_DSTATE[1] = rtb_XAxis1;
  stateEstimator_DW.UnitDelay1_DSTATE[2] = tmp1;

  // Update for UnitDelay: '<S26>/Delay Input1' incorporates:
  //   RelationalOperator: '<S28>/Compare'
  //
  //  Block description for '<S26>/Delay Input1':
  //
  //   Store in Global RAM

  stateEstimator_DW.DelayInput1_DSTATE = rtb_Compare;
}

// Constructor
stateEstimator::stateEstimator():
  stateEstimator_DW()
{
  // Currently there is no constructor body generated.
}

// Destructor
stateEstimator::~stateEstimator()
{
  // Currently there is no destructor body generated.
}

//
// File trailer for generated code.
//
// [EOF]
//
