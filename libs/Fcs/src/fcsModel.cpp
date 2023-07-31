//
// File: fcsModel.cpp
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
#include "rtwtypes.h"
#include "fcsModel_types.h"
#include <cmath>
#include <cstring>
#include <array>

// Named constants for Chart: '<S4>/State Machine'
const uint8_T fcsModel_IN_inActive{ 1U };

const uint8_T fcsModel_IN_inFlight{ 2U };

const uint8_T fcsModel_IN_mtrArmed{ 3U };

//
// Output and update for atomic system:
//    '<S11>/pidWithDebug'
//    '<S52>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, real_T rtu_integratorReset, const busPidParams
  *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd, busPidDebug
  *rty_pidDebug, real_T rtp_sampleTime_s, DW_pidWithDebug_fcsModel_T *localDW)
{
  real_T normalizer;
  real_T rtb_DiscreteTransferFcn_k;
  real_T rtb_Product2;
  real_T rtb_Sum;
  real_T rtb_Switch2;
  real_T rtb_UkYk1;
  real_T rtb_UnitDelay_i;

  // Product: '<S46>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S46>/sample time'
  //
  //  About '<S46>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Product2 = rtu_pidParamBus->outputRateLimits[1] * 0.004;

  // Sum: '<S14>/Sum'
  rtb_Sum = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S14>/Discrete First Order Deriv Filter'
  // MATLAB Function: '<S45>/Compute Deriv Filter Numerator And Denominator'
  //  Call the main function
  // MATLAB Function 'Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator': '<S48>:1' 
  // '<S48>:1:4' [num, den] = computeFirstOrderDerivFilterNumAndDen_function(filterBandwidth_radps, sampleTime_s); 
  //  This function computes the numerator and denominator of the discrete
  //  first order derivative filter
  //
  // Inputs:
  // filterBandwidth_radps: Bandwidth of the filter
  // sampleTime_s: sampling time
  //
  // Outputs:
  // num: Numerator array for the discrete transfer function
  // den: Denominator array for the discrete transfer function
  // 'computeFirstOrderDerivFilterNumAndDen_function:13' B1 = filterBandwidth_radps; 
  // 'computeFirstOrderDerivFilterNumAndDen_function:14' B0 = 0;
  // 'computeFirstOrderDerivFilterNumAndDen_function:16' A0 = B1;
  // 'computeFirstOrderDerivFilterNumAndDen_function:17' A1 = 1;
  // 'computeFirstOrderDerivFilterNumAndDen_function:18' K = 2/sampleTime_s;
  rtb_UnitDelay_i = 2.0 / rtp_sampleTime_s;

  // 'computeFirstOrderDerivFilterNumAndDen_function:20' [num, den] = computeDiscreteTFNumAndDen_function([B0, B1], [A0, A1], K); 
  // COMPUTEDISCRETETFNUMANDDEN_FUNCTION computes the numerator and denominator
  //  for a first and second order discrete transfer function from it's
  //  continuous counterpart
  //
  //  Inputs:
  //  B: Array of coefficients of continuous transfer function numerator arranged 
  //  in ascending power of s
  //  A: Array of coefficients of continuous transfer function denominator arranged 
  //  in ascending power of s
  //  K: 2/sampling time in sec
  //
  //  Outputs:
  // num: numerator of the equivalent discrete transfer function in descending power of z 
  // den: denominator of the equivalent discrete transfer function in descending power of z 
  //  get the length of coefficient array to determine the order of transfer
  //  function
  // 'computeDiscreteTFNumAndDen_function:19' nArray = length(B);
  // 'computeDiscreteTFNumAndDen_function:21' if (nArray == 2)
  //  For 1st order system
  // 'computeDiscreteTFNumAndDen_function:23' normalizer = A(1) + A(2)*K;
  normalizer = rtu_pidParamBus->filterBandwidth_radps + rtb_UnitDelay_i;

  // 'computeDiscreteTFNumAndDen_function:24' b0 = (B(1) + B(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:25' b1 = (B(1) - B(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:27' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:28' a1 = (A(1) - A(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:29' num = [b0, b1];
  rtb_DiscreteTransferFcn_k = rtu_pidParamBus->filterBandwidth_radps *
    rtb_UnitDelay_i;
  localDW->num[0] = rtb_DiscreteTransferFcn_k / normalizer;
  localDW->num[1] = (0.0 - rtb_DiscreteTransferFcn_k) / normalizer;

  // 'computeDiscreteTFNumAndDen_function:30' den = [a0, a1];
  localDW->den[0] = 1.0;
  localDW->den[1] = (rtu_pidParamBus->filterBandwidth_radps - rtb_UnitDelay_i) /
    normalizer;

  // DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn'
  rtb_UnitDelay_i = rtb_Sum - localDW->den[1] *
    localDW->DiscreteTransferFcn_states;
  rtb_DiscreteTransferFcn_k = localDW->num[0] * rtb_UnitDelay_i + localDW->num[1]
    * localDW->DiscreteTransferFcn_states;

  // Update for DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states = rtb_UnitDelay_i;

  // End of Outputs for SubSystem: '<S14>/Discrete First Order Deriv Filter'

  // Product: '<S14>/Product'
  rtb_UnitDelay_i = rtb_DiscreteTransferFcn_k * rtu_pidParamBus->Kd;

  // Product: '<S14>/Product1'
  normalizer = rtb_Sum * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S14>/Discrete-Time Integrator'
  if ((rtu_integratorReset != 0.0) || (localDW->DiscreteTimeIntegrator_PrevRese
       != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = 0.0;
  }

  // Sum: '<S14>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S14>/Discrete-Time Integrator'

  rtb_DiscreteTransferFcn_k = ((rtu_feedForward + rtb_UnitDelay_i) + normalizer)
    + localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S47>/Switch2' incorporates:
  //   RelationalOperator: '<S47>/LowerRelop1'
  //   RelationalOperator: '<S47>/UpperRelop'
  //   Switch: '<S47>/Switch'

  if (rtb_DiscreteTransferFcn_k > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2 = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_DiscreteTransferFcn_k < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S47>/Switch'
    rtb_Switch2 = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2 = rtb_DiscreteTransferFcn_k;
  }

  // End of Switch: '<S47>/Switch2'

  // Sum: '<S46>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S46>/Delay Input2'
  //
  //  Block description for '<S46>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S46>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2 - localDW->DelayInput2_DSTATE;

  // Switch: '<S49>/Switch2' incorporates:
  //   RelationalOperator: '<S49>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Product2) {
    // Product: '<S46>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S46>/sample time'
    //
    //  About '<S46>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Product2 = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S49>/Switch' incorporates:
    //   RelationalOperator: '<S49>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Product2) {
      rtb_Product2 = rtb_UkYk1;
    }

    // End of Switch: '<S49>/Switch'
  }

  // End of Switch: '<S49>/Switch2'

  // Sum: '<S46>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S46>/Delay Input2'
  //
  //  Block description for '<S46>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S46>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Product2 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S14>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S14>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = normalizer;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_UnitDelay_i;

  // Update for DiscreteIntegrator: '<S14>/Discrete-Time Integrator' incorporates:
  //   Product: '<S14>/Product2'
  //   Product: '<S14>/Product3'
  //   Product: '<S14>/Product5'
  //   Sum: '<S14>/Sum2'
  //   Sum: '<S14>/Sum3'
  //   Sum: '<S14>/Sum4'
  //   Sum: '<S14>/Sum5'
  //   UnitDelay: '<S14>/Unit Delay'
  //   UnitDelay: '<S14>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum * rtu_pidParamBus->Ki) * 0.004;
  if (rtu_integratorReset > 0.0) {
    localDW->DiscreteTimeIntegrator_PrevRese = 1;
  } else if (rtu_integratorReset < 0.0) {
    localDW->DiscreteTimeIntegrator_PrevRese = -1;
  } else if (rtu_integratorReset == 0.0) {
    localDW->DiscreteTimeIntegrator_PrevRese = 0;
  } else {
    localDW->DiscreteTimeIntegrator_PrevRese = 2;
  }

  // End of Update for DiscreteIntegrator: '<S14>/Discrete-Time Integrator'

  // Update for UnitDelay: '<S46>/Delay Input2'
  //
  //  Block description for '<S46>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S14>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2;

  // Update for UnitDelay: '<S14>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_DiscreteTransferFcn_k;
}

//
// Output and update for atomic system:
//    '<S30>/Compute Natural Frequency'
//    '<S31>/Compute Natural Frequency'
//    '<S15>/Compute Natural Frequency'
//    '<S16>/Compute Natural Frequency'
//    '<S71>/Compute Natural Frequency'
//    '<S72>/Compute Natural Frequency'
//    '<S56>/Compute Natural Frequency'
//    '<S57>/Compute Natural Frequency'
//
void fcsModel::fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps,
  real_T rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps)
{
  real_T tmp;

  //  call the main function
  // MATLAB Function 'Discrete Second Order Filter/Compute Natural Frequency': '<S38>:1' 
  // '<S38>:1:4' naturalFrequency_radps = computeSecondOrderSystemNaturalFrequency_function(bandwidth_radps, dampingRatio_nd); 
  // COMPUTESECONDORDERSYSTEMNATURALFREQUENCY_FUNCTION computes the natural
  // frequency of a second order system when user provides damping ratio and
  // required bandwith.
  //
  // Input:
  // bandwidth_radps: Desired Bandwidth in rad/s
  // dampingRatio_nd: Damping ration of the system
  // 'computeSecondOrderSystemNaturalFrequency_function:9' naturalFrequency_radps = bandwidth_radps/(sqrt(1 - 2*dampingRatio_nd^2 + sqrt(2 - 4*dampingRatio_nd^2 + 4*dampingRatio_nd^4))); 
  tmp = rtu_dampingRatio_nd * rtu_dampingRatio_nd;
  *rty_naturalFrequency_radps = rtu_bandwidth_radps / std::sqrt(std::sqrt((2.0 -
    tmp * 4.0) + 4.0 * std::pow(rtu_dampingRatio_nd, 4.0)) + (1.0 - tmp * 2.0));
}

//
// Output and update for atomic system:
//    '<S11>/Signal Conditioning Block1'
//    '<S11>/Signal Conditioning Block'
//    '<S52>/Signal Conditioning Block1'
//    '<S52>/Signal Conditioning Block'
//
void fcsModel::fcsMod_SignalConditioningBlock1(real_T rtu_input, const
  busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
  rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW)
{
  real_T B0;
  real_T num_tmp;
  real_T rtb_Switch2;
  real_T rtb_Switch2_m_tmp;
  real_T rtb_UkYk1;

  // MATLAB Function: '<S30>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S31>/Compute Natural Frequency'
  //  call the main function
  // MATLAB Function 'Discrete Second Order Deriv Filter/Compute Numerator And Denominator': '<S39>:1' 
  // '<S39>:1:4' [rateNum, accelNum, den] = computeSecondOrderDerivFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
  // COMPUTESECONDORDERDERIVFILTERNUMANDDEN_FUNCTION % This function computes the numerator and denominator of the dicrete 
  //  second order derivative and double derivative filters
  //
  // Inputs:
  // naturalFrequency_radps: Natural frequency of the filter
  // dampingRation_nd: Damping Ration of the filter
  // sampleTime_s: sampling time
  //
  // Outputs:
  // rateNum: Numerator array for the discrete derivative transfer function
  // accelNum: Numerator array for the discrete double derivative transfer function 
  // den: Denominator array for the discrete derivative transfer functions
  // 'computeSecondOrderDerivFilterNumAndDen_function:15' K = 2/sampleTime_s;
  // 'computeSecondOrderDerivFilterNumAndDen_function:16' A0 = naturalFrequency_radps^2; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:17' A1 = 2*dampingRatio_nd*naturalFrequency_radps; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:18' A2 = 1;
  // 'computeSecondOrderDerivFilterNumAndDen_function:20' B0 = 0;
  // 'computeSecondOrderDerivFilterNumAndDen_function:22' B1 = naturalFrequency_radps^2; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:23' B2 = 0;
  //  compute the rate transfer function numerator and the denominator
  // 'computeSecondOrderDerivFilterNumAndDen_function:26' [rateNum, den] = computeDiscreteTFNumAndDen_function([B0, B1, B2], [A0, A1, A2], K); 
  // COMPUTEDISCRETETFNUMANDDEN_FUNCTION computes the numerator and denominator
  //  for a first and second order discrete transfer function from it's
  //  continuous counterpart
  //
  //  Inputs:
  //  B: Array of coefficients of continuous transfer function numerator arranged 
  //  in ascending power of s
  //  A: Array of coefficients of continuous transfer function denominator arranged 
  //  in ascending power of s
  //  K: 2/sampling time in sec
  //
  //  Outputs:
  // num: numerator of the equivalent discrete transfer function in descending power of z 
  // den: denominator of the equivalent discrete transfer function in descending power of z 
  //  get the length of coefficient array to determine the order of transfer
  //  function
  // 'computeDiscreteTFNumAndDen_function:19' nArray = length(B);
  // 'computeDiscreteTFNumAndDen_function:21' if (nArray == 2)
  // 'computeDiscreteTFNumAndDen_function:32' elseif (nArray == 3)
  //  For 2nd order system
  // 'computeDiscreteTFNumAndDen_function:34' normalizer = A(1) + A(2)*K + A(3)*K^2; 
  // 'computeDiscreteTFNumAndDen_function:35' b0 = (B(1) + B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:36' b1 = (2*B(1) - 2*B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:37' b2 =  (B(1) - B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:39' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:40' a1 = (2*A(1) - 2*A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:41' a2 = (A(1) - A(2)*K + A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:42' num = [b0, b1, b2];
  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
  // 'computeSecondOrderDerivFilterNumAndDen_function:28' B1 = 0;
  // 'computeSecondOrderDerivFilterNumAndDen_function:29' B2 = naturalFrequency_radps^2; 
  //  compute the accel transfer function numerator
  // 'computeSecondOrderDerivFilterNumAndDen_function:31' [accelNum, ~] = computeDiscreteTFNumAndDen_function([B0, B1, B2], [A0, A1, A2], K); 
  // COMPUTEDISCRETETFNUMANDDEN_FUNCTION computes the numerator and denominator
  //  for a first and second order discrete transfer function from it's
  //  continuous counterpart
  //
  //  Inputs:
  //  B: Array of coefficients of continuous transfer function numerator arranged 
  //  in ascending power of s
  //  A: Array of coefficients of continuous transfer function denominator arranged 
  //  in ascending power of s
  //  K: 2/sampling time in sec
  //
  //  Outputs:
  // num: numerator of the equivalent discrete transfer function in descending power of z 
  // den: denominator of the equivalent discrete transfer function in descending power of z 
  //  get the length of coefficient array to determine the order of transfer
  //  function
  // 'computeDiscreteTFNumAndDen_function:19' nArray = length(B);
  // 'computeDiscreteTFNumAndDen_function:21' if (nArray == 2)
  // 'computeDiscreteTFNumAndDen_function:32' elseif (nArray == 3)
  //  For 2nd order system
  // 'computeDiscreteTFNumAndDen_function:34' normalizer = A(1) + A(2)*K + A(3)*K^2; 
  // 'computeDiscreteTFNumAndDen_function:35' b0 = (B(1) + B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:36' b1 = (2*B(1) - 2*B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:37' b2 =  (B(1) - B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:39' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:40' a1 = (2*A(1) - 2*A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:41' a2 = (A(1) - A(2)*K + A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:42' num = [b0, b1, b2];
  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S31>/Compute Filter Numerator And Denominator'
  //  Call the main function
  // MATLAB Function 'Discrete Second Order Filter/Compute Filter Numerator And Denominator': '<S40>:1' 
  // '<S40>:1:4' [num, den] = computeSecondOrderFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
  //  This function computes the numerator and denominator of the dicrete
  //  second order filter
  //
  // Inputs:
  // naturalFrequency_radps: Natural frequency of the filter
  // dampingRation_nd: Damping Ration of the filter
  // sampleTime_s: sampling time
  //
  // Outputs:
  // num: Numerator array for the discrete transfer function
  // den: Denominator array for the discrete transfer function
  // 'computeSecondOrderFilterNumAndDen_function:14' B0 = naturalFrequency_radps^2; 
  B0 = rtb_Switch2 * rtb_Switch2;

  // 'computeSecondOrderFilterNumAndDen_function:15' B1 = 0;
  // 'computeSecondOrderFilterNumAndDen_function:16' B2 = 0;
  // 'computeSecondOrderFilterNumAndDen_function:18' A0 = B0;
  // 'computeSecondOrderFilterNumAndDen_function:19' A1 = 2*dampingRatio_nd*naturalFrequency_radps; 
  // 'computeSecondOrderFilterNumAndDen_function:20' A2 = 1;
  // 'computeSecondOrderFilterNumAndDen_function:22' K = 2/sampleTime_s;
  rtb_UkYk1 = 2.0 / rtp_sampleTime_s;

  // 'computeSecondOrderFilterNumAndDen_function:24' [num, den] = computeDiscreteTFNumAndDen_function([B0, B1, B2], [A0, A1, A2], K); 
  // COMPUTEDISCRETETFNUMANDDEN_FUNCTION computes the numerator and denominator
  //  for a first and second order discrete transfer function from it's
  //  continuous counterpart
  //
  //  Inputs:
  //  B: Array of coefficients of continuous transfer function numerator arranged 
  //  in ascending power of s
  //  A: Array of coefficients of continuous transfer function denominator arranged 
  //  in ascending power of s
  //  K: 2/sampling time in sec
  //
  //  Outputs:
  // num: numerator of the equivalent discrete transfer function in descending power of z 
  // den: denominator of the equivalent discrete transfer function in descending power of z 
  //  get the length of coefficient array to determine the order of transfer
  //  function
  // 'computeDiscreteTFNumAndDen_function:19' nArray = length(B);
  // 'computeDiscreteTFNumAndDen_function:21' if (nArray == 2)
  // 'computeDiscreteTFNumAndDen_function:32' elseif (nArray == 3)
  //  For 2nd order system
  // 'computeDiscreteTFNumAndDen_function:34' normalizer = A(1) + A(2)*K + A(3)*K^2; 
  rtb_Switch2_m_tmp = rtb_UkYk1 * rtb_UkYk1;
  rtb_UkYk1 *= 2.0 * rtu_params->filterParams.dampingRatio_nd * rtb_Switch2;
  rtb_Switch2 = (rtb_UkYk1 + B0) + rtb_Switch2_m_tmp;

  // 'computeDiscreteTFNumAndDen_function:35' b0 = (B(1) + B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:36' b1 = (2*B(1) - 2*B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:37' b2 =  (B(1) - B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:39' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:40' a1 = (2*A(1) - 2*A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:41' a2 = (A(1) - A(2)*K + A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:42' num = [b0, b1, b2];
  num_tmp = B0 / rtb_Switch2;
  localDW->num[0] = num_tmp;
  localDW->num[1] = 2.0 * B0 / rtb_Switch2;
  localDW->num[2] = num_tmp;

  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
  localDW->den[0] = 1.0;
  localDW->den[1] = (2.0 * B0 - rtb_Switch2_m_tmp * 2.0) / rtb_Switch2;
  localDW->den[2] = ((B0 - rtb_UkYk1) + rtb_Switch2_m_tmp) / rtb_Switch2;

  // DiscreteTransferFcn: '<S31>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_Switch2 = (localDW->num[0] * localDW->DiscreteTransferFcn_tmp +
                 localDW->DiscreteTransferFcn_states[0] * localDW->num[1]) +
    localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S35>/Switch2' incorporates:
  //   RelationalOperator: '<S35>/LowerRelop1'
  //   RelationalOperator: '<S35>/UpperRelop'
  //   Switch: '<S35>/Switch'

  if (rtb_Switch2 > rtu_params->filteredInputLimits[1]) {
    rtb_Switch2 = rtu_params->filteredInputLimits[1];
  } else if (rtb_Switch2 < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S35>/Switch'
    rtb_Switch2 = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S35>/Switch2'

  // Sum: '<S32>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S32>/Delay Input2'
  //
  //  Block description for '<S32>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S32>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2 - localDW->DelayInput2_DSTATE;

  // Switch: '<S42>/Switch2' incorporates:
  //   Product: '<S32>/delta rise limit'
  //   SampleTimeMath: '<S32>/sample time'
  //
  //  About '<S32>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S42>/Switch2' incorporates:
  //   RelationalOperator: '<S42>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Switch2) {
    // Product: '<S32>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S32>/sample time'
    //
    //  About '<S32>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S42>/Switch' incorporates:
    //   RelationalOperator: '<S42>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Switch2) {
      // Switch: '<S42>/Switch2'
      rtb_Switch2 = rtb_UkYk1;
    }

    // End of Switch: '<S42>/Switch'
  }

  // End of Switch: '<S42>/Switch2'

  // Sum: '<S32>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S32>/Delay Input2'
  //
  //  Block description for '<S32>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S32>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S31>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S32>/Delay Input2'
  //
  //  Block description for '<S32>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_filteredInput;
}

//
// Function for Chart: '<S4>/State Machine'
// function isStickLow = getStickStates(rcInCmds, rcParamsStruct)
//
boolean_T fcsModel::fcsModel_getStickStates(uint16_T rcInCmds_throttleCmd_nd,
  uint16_T rcInCmds_joystickYCmd_nd, uint16_T rcInCmds_joystickXCmd_nd, uint16_T
  rcInCmds_joystickZCmd_nd, const struct_ek64hRkWPGOrPAEJD7wTRH
  *b_rcParamsStruct)
{
  boolean_T isStickLow;

  // MATLAB Function 'getStickStates': '<S93>:10'
  // '<S93>:10:2' if(rcInCmds.throttleCmd_nd <= rcParamsStruct.lowPwmThreshold(2) && rcInCmds.throttleCmd_nd >= rcParamsStruct.lowPwmThreshold(1) && ... 
  // '<S93>:10:3'        rcInCmds.joystickXCmd_nd <= rcParamsStruct.lowPwmThreshold(2) && rcInCmds.throttleCmd_nd >= rcParamsStruct.lowPwmThreshold(1) && ... 
  // '<S93>:10:4'        rcInCmds.joystickYCmd_nd <= rcParamsStruct.lowPwmThreshold(2) && rcInCmds.throttleCmd_nd >= rcParamsStruct.lowPwmThreshold(1) && ... 
  // '<S93>:10:5'        rcInCmds.joystickZCmd_nd <= rcParamsStruct.lowPwmThreshold(2) && rcInCmds.throttleCmd_nd >= rcParamsStruct.lowPwmThreshold(1)) 
  if (rcInCmds_throttleCmd_nd <= b_rcParamsStruct->lowPwmThreshold[1]) {
    if (rcInCmds_throttleCmd_nd >= b_rcParamsStruct->lowPwmThreshold[0]) {
      if (rcInCmds_joystickXCmd_nd <= b_rcParamsStruct->lowPwmThreshold[1]) {
        if (rcInCmds_throttleCmd_nd >= b_rcParamsStruct->lowPwmThreshold[0]) {
          if (rcInCmds_joystickYCmd_nd <= b_rcParamsStruct->lowPwmThreshold[1])
          {
            if (rcInCmds_throttleCmd_nd >= b_rcParamsStruct->lowPwmThreshold[0])
            {
              if (rcInCmds_joystickZCmd_nd <= b_rcParamsStruct->lowPwmThreshold
                  [1]) {
                if (rcInCmds_throttleCmd_nd >= b_rcParamsStruct->
                    lowPwmThreshold[0]) {
                  // '<S93>:10:6' isStickLow = true;
                  isStickLow = true;
                } else {
                  // '<S93>:10:7' else
                  // '<S93>:10:8' isStickLow = false;
                  isStickLow = false;
                }
              } else {
                // '<S93>:10:7' else
                // '<S93>:10:8' isStickLow = false;
                isStickLow = false;
              }
            } else {
              // '<S93>:10:7' else
              // '<S93>:10:8' isStickLow = false;
              isStickLow = false;
            }
          } else {
            // '<S93>:10:7' else
            // '<S93>:10:8' isStickLow = false;
            isStickLow = false;
          }
        } else {
          // '<S93>:10:7' else
          // '<S93>:10:8' isStickLow = false;
          isStickLow = false;
        }
      } else {
        // '<S93>:10:7' else
        // '<S93>:10:8' isStickLow = false;
        isStickLow = false;
      }
    } else {
      // '<S93>:10:7' else
      // '<S93>:10:8' isStickLow = false;
      isStickLow = false;
    }
  } else {
    // '<S93>:10:7' else
    // '<S93>:10:8' isStickLow = false;
    isStickLow = false;
  }

  return isStickLow;
}

// Model step function
void fcsModel::step()
{
  std::array<real_T, 4> DiscreteTransferFcn_tmp;
  std::array<real_T, 4> rtb_DiscreteTransferFcn;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDeb_b;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDebug;
  std::array<busCtrlInputs, 3> rtb_ctrlInputsArray;
  std::array<real_T, 3> rtb_momCmd_Nm;
  busPidDebug rtb_BusCreator_o;
  real_T pCmd;
  real_T rCmd;
  real_T tCmd;
  real_T tmp;
  real_T tmp_0;
  int32_T yCmd;
  enumStateMachine rtb_state;

  // Chart: '<S4>/State Machine' incorporates:
  //   BusCreator generated from: '<S4>/State Machine'
  //   Inport: '<Root>/rcCmdsIn'

  if (fcsModel_DW.temporalCounter_i1 < 8191U) {
    fcsModel_DW.temporalCounter_i1 = static_cast<uint16_T>
      (fcsModel_DW.temporalCounter_i1 + 1U);
  }

  // Gateway: rcInterpreter/State Machine
  // During: rcInterpreter/State Machine
  if (fcsModel_DW.is_active_c4_rcInterpreter == 0U) {
    // Entry: rcInterpreter/State Machine
    fcsModel_DW.is_active_c4_rcInterpreter = 1U;

    // Entry Internal: rcInterpreter/State Machine
    // Transition: '<S93>:4'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c4_rcInterpreter = fcsModel_IN_inActive;

    // Entry 'inActive': '<S93>:3'
    // '<S93>:3:3' state = enumStateMachine.INACTIVE;
    rtb_state = enumStateMachine::INACTIVE;

    // '<S93>:3:4' integratorReset = 1;
    fcsModel_DW.integratorReset = 1.0;
  } else {
    switch (fcsModel_DW.is_c4_rcInterpreter) {
     case fcsModel_IN_inActive:
      rtb_state = enumStateMachine::INACTIVE;

      // During 'inActive': '<S93>:3'
      // '<S93>:6:1' sf_internal_predicateOutput = duration(getStickStates(rcInCmds, rcParamsStruct), sec) >= 5; 
      if (!fcsModel_getStickStates(fcsModel_U.rcCmdsIn.throttleCmd_nd,
           fcsModel_U.rcCmdsIn.joystickYCmd_nd,
           fcsModel_U.rcCmdsIn.joystickXCmd_nd,
           fcsModel_U.rcCmdsIn.joystickZCmd_nd,
           &fcsModel_ConstP.StateMachine_rcParamsStruct)) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if (fcsModel_DW.durationCounter_1 >= 1250) {
        // Transition: '<S93>:6'
        fcsModel_DW.durationCounter_1_e = 0;
        fcsModel_DW.is_c4_rcInterpreter = fcsModel_IN_mtrArmed;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'mtrArmed': '<S93>:5'
        // '<S93>:5:3' state = enumStateMachine.MTR_ARMED;
        rtb_state = enumStateMachine::MTR_ARMED;

        // '<S93>:5:4' integratorReset = 1;
        fcsModel_DW.integratorReset = 1.0;
      }
      break;

     case fcsModel_IN_inFlight:
      rtb_state = enumStateMachine::INFLIGHT;

      // During 'inFlight': '<S93>:15'
      // '<S93>:19:1' sf_internal_predicateOutput = rcInCmds.throttleCmd_nd < rcParamsStruct.lowPwmThreshold(2); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd < 990) {
        // Transition: '<S93>:19'
        fcsModel_DW.durationCounter_1_e = 0;
        fcsModel_DW.is_c4_rcInterpreter = fcsModel_IN_mtrArmed;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'mtrArmed': '<S93>:5'
        // '<S93>:5:3' state = enumStateMachine.MTR_ARMED;
        rtb_state = enumStateMachine::MTR_ARMED;

        // '<S93>:5:4' integratorReset = 1;
        fcsModel_DW.integratorReset = 1.0;
      }
      break;

     default:
      {
        rtb_state = enumStateMachine::MTR_ARMED;

        // During 'mtrArmed': '<S93>:5'
        // '<S93>:18:1' sf_internal_predicateOutput = rcInCmds.throttleCmd_nd > rcParamsStruct.lowPwmThreshold(2); 
        if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 990) {
          // Transition: '<S93>:18'
          // Exit 'mtrArmed': '<S93>:5'
          // '<S93>:5:6' integratorReset = 0;
          fcsModel_DW.integratorReset = 0.0;
          fcsModel_DW.is_c4_rcInterpreter = fcsModel_IN_inFlight;

          // Entry 'inFlight': '<S93>:15'
          // '<S93>:15:3' state = enumStateMachine.INFLIGHT;
          rtb_state = enumStateMachine::INFLIGHT;
        } else {
          boolean_T rtb_Compare_f;

          // '<S93>:14:1' sf_internal_predicateOutput = after(20, sec) || duration(getStickStates(rcInCmds, rcParamsStruct), sec) >= 3; 
          if (fcsModel_DW.temporalCounter_i1 >= 5000U) {
            rtb_Compare_f = true;
          } else {
            if (!fcsModel_getStickStates(fcsModel_U.rcCmdsIn.throttleCmd_nd,
                 fcsModel_U.rcCmdsIn.joystickYCmd_nd,
                 fcsModel_U.rcCmdsIn.joystickXCmd_nd,
                 fcsModel_U.rcCmdsIn.joystickZCmd_nd,
                 &fcsModel_ConstP.StateMachine_rcParamsStruct)) {
              fcsModel_DW.durationCounter_1_e = 0;
            }

            rtb_Compare_f = (fcsModel_DW.durationCounter_1_e >= 750);
          }

          if (rtb_Compare_f) {
            // Transition: '<S93>:14'
            // Exit 'mtrArmed': '<S93>:5'
            // '<S93>:5:6' integratorReset = 0;
            fcsModel_DW.durationCounter_1 = 0;
            fcsModel_DW.is_c4_rcInterpreter = fcsModel_IN_inActive;

            // Entry 'inActive': '<S93>:3'
            // '<S93>:3:3' state = enumStateMachine.INACTIVE;
            rtb_state = enumStateMachine::INACTIVE;

            // '<S93>:3:4' integratorReset = 1;
            fcsModel_DW.integratorReset = 1.0;
          }
        }
      }
      break;
    }
  }

  if (fcsModel_getStickStates(fcsModel_U.rcCmdsIn.throttleCmd_nd,
       fcsModel_U.rcCmdsIn.joystickYCmd_nd, fcsModel_U.rcCmdsIn.joystickXCmd_nd,
       fcsModel_U.rcCmdsIn.joystickZCmd_nd,
       &fcsModel_ConstP.StateMachine_rcParamsStruct)) {
    fcsModel_DW.durationCounter_1++;
    fcsModel_DW.durationCounter_1_e++;
  } else {
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.durationCounter_1_e = 0;
  }

  // End of Chart: '<S4>/State Machine'

  // MATLAB Function: '<S4>/Interpret RC In Cmds' incorporates:
  //   BusCreator generated from: '<S4>/Interpret RC In Cmds'
  //   Inport: '<Root>/rcCmdsIn'

  // Computes command and flight mode from the rc inputs
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S92>:1'
  // '<S92>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcInCmds, rcParamsStruct); 
  // INTERPRETRCINPUTS_FUNCTION %Computes command and flight mode from the rc inputs 
  // 'interpretRcInputs_function:3' tCoeffs = rcParamsStruct.coeffs.throttle_nd; 
  // 'interpretRcInputs_function:4' tCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:5'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.throttleCmd_nd) ) ); 
  tCmd = std::fmin(2006.0, std::fmax(990.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.throttleCmd_nd)));

  // 'interpretRcInputs_function:6' rcOutCmds.throttleStick = -(tCoeffs(1)*tCmd^3 + tCoeffs(2)*tCmd^2 + ... 
  // 'interpretRcInputs_function:7'                           tCoeffs(3)*tCmd + tCoeffs(4))*rcParamsStruct.cmdLimits.zForce_N(2); 
  // 'interpretRcInputs_function:9' rCoeffs = rcParamsStruct.coeffs.roll_nd;
  // 'interpretRcInputs_function:10' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:11'        max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  rCmd = std::fmin(2006.0, std::fmax(990.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickXCmd_nd)));

  // 'interpretRcInputs_function:12' rcOutCmds.rollStick = (rCoeffs(1)*rCmd^3 + rCoeffs(2)*rCmd^2 + ... 
  // 'interpretRcInputs_function:13'                           rCoeffs(3)*rCmd + rCoeffs(4))*rcParamsStruct.cmdLimits.rollRate_radps(2); 
  // 'interpretRcInputs_function:15' pCoeffs = rcParamsStruct.coeffs.pitch_nd;
  // 'interpretRcInputs_function:16' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:17'        max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickYCmd_nd) ) ); 
  pCmd = std::fmin(2006.0, std::fmax(990.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickYCmd_nd)));

  // 'interpretRcInputs_function:18' rcOutCmds.pitchStick = (pCoeffs(1)*pCmd^3 + pCoeffs(2)*pCmd^2 + ... 
  // 'interpretRcInputs_function:19'                           pCoeffs(3)*pCmd + pCoeffs(4))*rcParamsStruct.cmdLimits.pitchRate_radps(2); 
  // 'interpretRcInputs_function:21' yCoeffs = rcParamsStruct.coeffs.yaw_nd;
  // 'interpretRcInputs_function:22' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:23'        max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  yCmd = static_cast<int32_T>(std::fmin(2006.0, std::fmax(990.0, static_cast<
    real_T>(fcsModel_U.rcCmdsIn.joystickZCmd_nd))));

  // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
  //   Constant: '<S3>/Constant'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'

  // 'interpretRcInputs_function:24' rcOutCmds.yawStick = (yCoeffs(1)*yCmd^3 + yCoeffs(2)*yCmd^2 + ... 
  // 'interpretRcInputs_function:25'                           yCoeffs(3)*yCmd + yCoeffs(4))*rcParamsStruct.cmdLimits.yawRate_radps(2); 
  // if(rcInCmds.rcSwitch1_nd < 1500)
  // 'interpretRcInputs_function:28' flightMode = enumFlightMode.STABILIZE;
  // end
  std::memcpy(&rtb_ctrlInputsArray[0],
              &fcsModel_ConstP.Constant_Value_h.attCtrlInputs.ctrlInputsArray[0],
              3U * sizeof(busCtrlInputs));

  // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S91>:1' 
  // '<S91>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
  //  This is a stop gap setup where we are only assuming that rate control
  //  is active and therefore not setting up attCtrlInputs for Euler angle
  //  control
  // '<S91>:1:6' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
  rtb_ctrlInputsArray[0].cmd = (((rCmd * rCmd * -1.7538875049548508E-5 +
    3.9045645325900015E-9 * std::pow(rCmd, 3.0)) + 0.027221770655502103 * rCmd)
    + -14.548300722291762) * 0.52359877559829882;

  // '<S91>:1:7' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
  rtb_ctrlInputsArray[1].cmd = (((pCmd * pCmd * -1.7538875049548508E-5 +
    3.9045645325900015E-9 * std::pow(pCmd, 3.0)) + 0.027221770655502103 * pCmd)
    + -14.548300722291762) * 0.52359877559829882;

  // '<S91>:1:8' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
  rtb_ctrlInputsArray[2].cmd = (((static_cast<real_T>(yCmd * yCmd) *
    -1.7538875049548508E-5 + 3.9045645325900015E-9 * std::pow(static_cast<real_T>
    (yCmd), 3.0)) + 0.027221770655502103 * static_cast<real_T>(yCmd)) +
    -14.548300722291762) * 0.52359877559829882;

  // Outputs for Iterator SubSystem: '<S10>/For Each Subsystem' incorporates:
  //   ForEach: '<S52>/For Each'

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S52>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[0],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[0].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[0],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[0].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S52>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S52>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S52>/Unit Delay'

  fcsModel_pidWithDebug(0.0, pCmd, rCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[0],
                        fcsModel_DW.CoreSubsys_p[0].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[0].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S52>/pidWithDebug'

  // Update for UnitDelay: '<S52>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[0].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S52>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDebug[0] = rtb_BusCreator_o;

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S52>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[1],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[1].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[1],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[1].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S52>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S52>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S52>/Unit Delay'

  fcsModel_pidWithDebug(0.0, pCmd, rCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[1],
                        fcsModel_DW.CoreSubsys_p[1].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[1].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S52>/pidWithDebug'

  // Update for UnitDelay: '<S52>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[1].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S52>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDebug[1] = rtb_BusCreator_o;

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S52>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[2],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[2].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[2],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[2].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S52>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S52>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S52>/Unit Delay'

  fcsModel_pidWithDebug(0.0, pCmd, rCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[2],
                        fcsModel_DW.CoreSubsys_p[2].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[2].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S52>/pidWithDebug'

  // Update for UnitDelay: '<S52>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[2].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S52>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDebug[2] = rtb_BusCreator_o;

  // End of Outputs for SubSystem: '<S10>/For Each Subsystem'

  // MATLAB Function: '<S9>/EulerRates2BodyRates' incorporates:
  //   Inport: '<Root>/stateEstimate'

  // MATLAB Function 'EulerRates2BodyRates': '<S51>:1'
  // '<S51>:1:3' bodyRates_radps = eulerRates2bodyRates_function(taitBryanRates_radps,shipOrientation_rad); 
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Converts  rate of change of TaitBryan angles in the globle frame to
  // rotational rate of change in the body frame
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Inputs:
  // -taitBryanRates_radps:= {3x1 [omegaX;omegaY;omegaZ] describing the rate
  // of change in the space axes of the TaitBryan angles}
  // -shipOrientation_rad:= {3x1 [phi;theta;psi] vector describing orentation
  // of the body frame in the space frame}
  // Output:
  // -bodyRates_radps:= {3x1 [omegaX;omegaY;omegaZ] vector of rotation
  // rates of the body in the body frame}
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 'eulerRates2bodyRates_function:18' roll = shipOrientation_rad(1);
  // 'eulerRates2bodyRates_function:19' pitch = shipOrientation_rad(2);
  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  rCmd = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(rCmd - 1.5707963267948966) <= 0.004245395477824045) || (std::abs
       (rCmd - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        pCmd = -1.0;
      } else {
        pCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd = (rCmd - 0.004245395477824045) * pCmd;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        pCmd = -1.0;
      } else {
        pCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd = (rCmd - 0.004245395477824045) * pCmd;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        pCmd = -1.0;
      } else {
        pCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd = (rCmd + 0.004245395477824045) * pCmd;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  pCmd = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);

  // End of MATLAB Function: '<S9>/EulerRates2BodyRates'

  // Switch: '<S9>/Switch'
  // 'eulerRates2bodyRates_function:37' conversionMatrix = zeroSmallValues(conversionMatrix,eps); 
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // sets values in the M = zero if abs(values) is below this_eps
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Inputs: M:={Any real valued verable}, eps:={values in M below the abs of this 
  // esp are set to zero}
  // Ouputs: M:={with values below the esp set to zero}
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 'zeroSmallValues:10' for ii=1:size(M,1)
  // Convert rotation rate to change in TaitBryan angles
  // 'eulerRates2bodyRates_function:40' bodyRates_radps =  conversionMatrix * taitBryanRates_radps; 
  rtb_ImpAsg_InsertedFor_angAccel[0] = rtb_ctrlInputsArray[0].cmd;
  rtb_ImpAsg_InsertedFor_angAccel[1] = rtb_ctrlInputsArray[1].cmd;
  rtb_ImpAsg_InsertedFor_angAccel[2] = rtb_ctrlInputsArray[2].cmd;
  for (yCmd = 0; yCmd < 3; yCmd++) {
    // Switch: '<S9>/Switch' incorporates:
    //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    //   Product: '<S2>/Matrix Multiply'
    //
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    rtb_momCmd_Nm[yCmd] = rtb_ImpAsg_InsertedFor_angAccel[yCmd];
  }

  // BusCreator: '<S9>/Bus Creator' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_ctrlInputsArray[0].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[0].cmd = rtb_momCmd_Nm[0];
  rtb_ctrlInputsArray[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_ctrlInputsArray[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S9>/Bus Creator2' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_ctrlInputsArray[1].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[1].cmd = rtb_momCmd_Nm[1];
  rtb_ctrlInputsArray[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_ctrlInputsArray[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S9>/Bus Creator3' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_ctrlInputsArray[2].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[2].cmd = rtb_momCmd_Nm[2];
  rtb_ctrlInputsArray[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_ctrlInputsArray[2].trackingCtrlCmd = 0.0;

  // Outputs for Atomic SubSystem: '<S2>/Angular Rate Controller'
  // Outputs for Iterator SubSystem: '<S8>/For Each Subsystem' incorporates:
  //   ForEach: '<S11>/For Each'

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S11>/ctrlInputs' incorporates:
  //   BusCreator: '<S9>/Bus Creator1'
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[0],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[0].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[0],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[0].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S11>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S11>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S11>/Unit Delay'

  fcsModel_pidWithDebug(0.0, rCmd, pCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        0], fcsModel_DW.CoreSubsys[0].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[0].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S11>/pidWithDebug'

  // Update for UnitDelay: '<S11>/Unit Delay'
  fcsModel_DW.CoreSubsys[0].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S11>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDeb_b[0] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[0] = rCmd;

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S11>/ctrlInputs' incorporates:
  //   BusCreator: '<S9>/Bus Creator1'
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[1],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[1].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[1],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[1].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S11>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S11>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S11>/Unit Delay'

  fcsModel_pidWithDebug(0.0, rCmd, pCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        1], fcsModel_DW.CoreSubsys[1].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[1].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S11>/pidWithDebug'

  // Update for UnitDelay: '<S11>/Unit Delay'
  fcsModel_DW.CoreSubsys[1].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S11>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDeb_b[1] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[1] = rCmd;

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S11>/ctrlInputs' incorporates:
  //   BusCreator: '<S9>/Bus Creator1'
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[2],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[2].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[2],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[2].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S11>/pidWithDebug'
  // ForEachSliceSelector generated from: '<S11>/pidParams' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S11>/Unit Delay'

  fcsModel_pidWithDebug(0.0, rCmd, pCmd, fcsModel_DW.integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        2], fcsModel_DW.CoreSubsys[2].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[2].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S11>/pidWithDebug'

  // Update for UnitDelay: '<S11>/Unit Delay'
  fcsModel_DW.CoreSubsys[2].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S11>/pidDebug'
  rtb_ImpAsg_InsertedFor_pidDeb_b[2] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[2] = rCmd;

  // End of Outputs for SubSystem: '<S8>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'

  for (yCmd = 0; yCmd < 3; yCmd++) {
    rtb_momCmd_Nm[yCmd] = 0.0;
    rtb_momCmd_Nm[yCmd] += fcsModel_ConstP.Constant_Value_n[yCmd] *
      rtb_ImpAsg_InsertedFor_angAccel[0];
    rtb_momCmd_Nm[yCmd] += fcsModel_ConstP.Constant_Value_n[yCmd + 3] *
      rtb_ImpAsg_InsertedFor_angAccel[1];
    rtb_momCmd_Nm[yCmd] += fcsModel_ConstP.Constant_Value_n[yCmd + 6] *
      rtb_ImpAsg_InsertedFor_angAccel[2];
  }

  // End of Product: '<S2>/Matrix Multiply'

  // SignalConversion generated from: '<S1>/Matrix Multiply' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'

  rCmd = -(((tCmd * tCmd * -8.4651009094729818E-6 + 1.8913975371836569E-9 * std::
             pow(tCmd, 3.0)) + 0.013109462896716882 * tCmd) +
           -6.5014878863088388) * 50.0;
  pCmd = rtb_momCmd_Nm[0];
  tmp_0 = rtb_momCmd_Nm[1];
  tmp = rtb_momCmd_Nm[2];

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  for (yCmd = 0; yCmd < 4; yCmd++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   Constant: '<S1>/Constant'

    tCmd = ((fcsModel_ConstP.Constant_Value_c[yCmd + 4] * pCmd +
             fcsModel_ConstP.Constant_Value_c[yCmd] * rCmd) +
            fcsModel_ConstP.Constant_Value_c[yCmd + 8] * tmp_0) +
      fcsModel_ConstP.Constant_Value_c[yCmd + 12] * tmp;

    // Saturate: '<S1>/Saturation'
    if (tCmd > 1.5791367041742974E+6) {
      tCmd = 1.5791367041742974E+6;
    } else if (tCmd < 0.0) {
      tCmd = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    tCmd = 9.5492965855137211 * std::sqrt(tCmd) - -0.92734095767679814 *
      fcsModel_DW.DiscreteTransferFcn_states[yCmd];
    rtb_DiscreteTransferFcn[yCmd] = 0.036329521161600868 * tCmd +
      0.036329521161600868 * fcsModel_DW.DiscreteTransferFcn_states[yCmd];
    DiscreteTransferFcn_tmp[yCmd] = tCmd;
  }

  // Switch: '<S1>/Switch1' incorporates:
  //   Constant: '<S6>/Constant'
  //   Constant: '<S7>/Constant'
  //   RelationalOperator: '<S6>/Compare'
  //   RelationalOperator: '<S7>/Compare'
  //   Switch: '<S1>/Switch'

  if (rtb_state == enumStateMachine::INACTIVE) {
    // Outport: '<Root>/actuatorsCmds' incorporates:
    //   Constant: '<S1>/Constant2'

    fcsModel_Y.actuatorsCmds[0] = -1.0;
    fcsModel_Y.actuatorsCmds[1] = -1.0;
    fcsModel_Y.actuatorsCmds[2] = -1.0;
    fcsModel_Y.actuatorsCmds[3] = -1.0;
  } else if (rtb_state == enumStateMachine::MTR_ARMED) {
    // Switch: '<S1>/Switch' incorporates:
    //   Constant: '<S1>/Constant1'
    //   Outport: '<Root>/actuatorsCmds'

    fcsModel_Y.actuatorsCmds[0] = 0.0;
    fcsModel_Y.actuatorsCmds[1] = 0.0;
    fcsModel_Y.actuatorsCmds[2] = 0.0;
    fcsModel_Y.actuatorsCmds[3] = 0.0;
  } else {
    // Outport: '<Root>/actuatorsCmds' incorporates:
    //   Switch: '<S1>/Switch'

    fcsModel_Y.actuatorsCmds[0] = rtb_DiscreteTransferFcn[0];
    fcsModel_Y.actuatorsCmds[1] = rtb_DiscreteTransferFcn[1];
    fcsModel_Y.actuatorsCmds[2] = rtb_DiscreteTransferFcn[2];
    fcsModel_Y.actuatorsCmds[3] = rtb_DiscreteTransferFcn[3];
  }

  // End of Switch: '<S1>/Switch1'

  // Outport: '<Root>/fcsDebug' incorporates:
  //   BusAssignment: '<Root>/Bus Assignment'
  //   BusCreator: '<S2>/Bus Creator'
  //   ForEachSliceAssignment generated from: '<S11>/pidDebug'
  //   ForEachSliceAssignment generated from: '<S52>/pidDebug'

  std::memcpy(&fcsModel_Y.fcsDebug.innerLoopCtrldebug.angRateCtrlDebug[0],
              &rtb_ImpAsg_InsertedFor_pidDeb_b[0], 3U * sizeof(busPidDebug));
  std::memcpy(&fcsModel_Y.fcsDebug.innerLoopCtrldebug.attCtrlDebug[0],
              &rtb_ImpAsg_InsertedFor_pidDebug[0], 3U * sizeof(busPidDebug));
  fcsModel_Y.fcsDebug.state = rtb_state;

  // Update for DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
  fcsModel_DW.DiscreteTransferFcn_states[0] = DiscreteTransferFcn_tmp[0];
  fcsModel_DW.DiscreteTransferFcn_states[1] = DiscreteTransferFcn_tmp[1];
  fcsModel_DW.DiscreteTransferFcn_states[2] = DiscreteTransferFcn_tmp[2];
  fcsModel_DW.DiscreteTransferFcn_states[3] = DiscreteTransferFcn_tmp[3];
}

// Model initialize function
void fcsModel::initialize()
{
  // (no initialization code required)
}

// Model terminate function
void fcsModel::terminate()
{
  // (no terminate code required)
}

// Constructor
fcsModel::fcsModel():
  fcsModel_U(),
  fcsModel_Y(),
  fcsModel_DW()
{
  // Currently there is no constructor body generated.
}

// Destructor
fcsModel::~fcsModel()
{
  // Currently there is no destructor body generated.
}

//
// File trailer for generated code.
//
// [EOF]
//
