//
// File: fcsModel.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.59
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Thu Sep 28 16:38:18 2023
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

// Named constants for Chart: '<S4>/Chart'
const uint8_T fcsModel_IN_ARM_MTRS{ 1U };

const uint8_T fcsModel_IN_INACTIVE{ 2U };

const uint8_T fcsModel_IN_INFLIGHT{ 3U };

//
// Output and update for atomic system:
//    '<S10>/pidWithDebug'
//    '<S51>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
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

  // Product: '<S45>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S45>/sample time'
  //
  //  About '<S45>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Product2 = rtu_pidParamBus->outputRateLimits[1] * 0.004;

  // Sum: '<S13>/Sum'
  rtb_Sum = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S13>/Discrete First Order Deriv Filter'
  // MATLAB Function: '<S44>/Compute Deriv Filter Numerator And Denominator'
  //  Call the main function
  // MATLAB Function 'Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator': '<S47>:1' 
  // '<S47>:1:4' [num, den] = computeFirstOrderDerivFilterNumAndDen_function(filterBandwidth_radps, sampleTime_s); 
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

  // DiscreteTransferFcn: '<S44>/Discrete Transfer Fcn'
  rtb_UnitDelay_i = rtb_Sum - localDW->den[1] *
    localDW->DiscreteTransferFcn_states;
  rtb_DiscreteTransferFcn_k = localDW->num[0] * rtb_UnitDelay_i + localDW->num[1]
    * localDW->DiscreteTransferFcn_states;

  // Update for DiscreteTransferFcn: '<S44>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states = rtb_UnitDelay_i;

  // End of Outputs for SubSystem: '<S13>/Discrete First Order Deriv Filter'

  // Product: '<S13>/Product'
  rtb_UnitDelay_i = rtb_DiscreteTransferFcn_k * rtu_pidParamBus->Kd;

  // Product: '<S13>/Product1'
  normalizer = rtb_Sum * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S13>/Discrete-Time Integrator'
  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = 0.0;
  }

  // Sum: '<S13>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rtb_DiscreteTransferFcn_k = ((rtu_feedForward + rtb_UnitDelay_i) + normalizer)
    + localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S46>/Switch2' incorporates:
  //   RelationalOperator: '<S46>/LowerRelop1'
  //   RelationalOperator: '<S46>/UpperRelop'
  //   Switch: '<S46>/Switch'

  if (rtb_DiscreteTransferFcn_k > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2 = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_DiscreteTransferFcn_k < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S46>/Switch'
    rtb_Switch2 = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2 = rtb_DiscreteTransferFcn_k;
  }

  // End of Switch: '<S46>/Switch2'

  // Sum: '<S45>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S45>/Delay Input2'
  //
  //  Block description for '<S45>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S45>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2 - localDW->DelayInput2_DSTATE;

  // Switch: '<S48>/Switch2' incorporates:
  //   RelationalOperator: '<S48>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Product2) {
    // Product: '<S45>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S45>/sample time'
    //
    //  About '<S45>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Product2 = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S48>/Switch' incorporates:
    //   RelationalOperator: '<S48>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Product2) {
      rtb_Product2 = rtb_UkYk1;
    }

    // End of Switch: '<S48>/Switch'
  }

  // End of Switch: '<S48>/Switch2'

  // Sum: '<S45>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S45>/Delay Input2'
  //
  //  Block description for '<S45>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S45>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Product2 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S13>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = normalizer;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_UnitDelay_i;

  // Update for DiscreteIntegrator: '<S13>/Discrete-Time Integrator' incorporates:
  //   Product: '<S13>/Product2'
  //   Product: '<S13>/Product3'
  //   Product: '<S13>/Product5'
  //   Sum: '<S13>/Sum2'
  //   Sum: '<S13>/Sum3'
  //   Sum: '<S13>/Sum4'
  //   Sum: '<S13>/Sum5'
  //   UnitDelay: '<S13>/Unit Delay'
  //   UnitDelay: '<S13>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum * rtu_pidParamBus->Ki) * 0.004;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S45>/Delay Input2'
  //
  //  Block description for '<S45>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S13>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2;

  // Update for UnitDelay: '<S13>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_DiscreteTransferFcn_k;
}

//
// Output and update for atomic system:
//    '<S29>/Compute Natural Frequency'
//    '<S30>/Compute Natural Frequency'
//    '<S14>/Compute Natural Frequency'
//    '<S15>/Compute Natural Frequency'
//    '<S70>/Compute Natural Frequency'
//    '<S71>/Compute Natural Frequency'
//    '<S55>/Compute Natural Frequency'
//    '<S56>/Compute Natural Frequency'
//
void fcsModel::fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps,
  real_T rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps)
{
  real_T tmp;

  //  call the main function
  // MATLAB Function 'Discrete Second Order Filter/Compute Natural Frequency': '<S37>:1' 
  // '<S37>:1:4' naturalFrequency_radps = computeSecondOrderSystemNaturalFrequency_function(bandwidth_radps, dampingRatio_nd); 
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
//    '<S10>/Signal Conditioning Block1'
//    '<S10>/Signal Conditioning Block'
//    '<S51>/Signal Conditioning Block1'
//    '<S51>/Signal Conditioning Block'
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

  // MATLAB Function: '<S29>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S30>/Compute Natural Frequency'
  //  call the main function
  // MATLAB Function 'Discrete Second Order Deriv Filter/Compute Numerator And Denominator': '<S38>:1' 
  // '<S38>:1:4' [rateNum, accelNum, den] = computeSecondOrderDerivFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
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

  // MATLAB Function: '<S30>/Compute Filter Numerator And Denominator'
  //  Call the main function
  // MATLAB Function 'Discrete Second Order Filter/Compute Filter Numerator And Denominator': '<S39>:1' 
  // '<S39>:1:4' [num, den] = computeSecondOrderFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
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

  // DiscreteTransferFcn: '<S30>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_Switch2 = (localDW->num[0] * localDW->DiscreteTransferFcn_tmp +
                 localDW->DiscreteTransferFcn_states[0] * localDW->num[1]) +
    localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S34>/Switch2' incorporates:
  //   RelationalOperator: '<S34>/LowerRelop1'
  //   RelationalOperator: '<S34>/UpperRelop'
  //   Switch: '<S34>/Switch'

  if (rtb_Switch2 > rtu_params->filteredInputLimits[1]) {
    rtb_Switch2 = rtu_params->filteredInputLimits[1];
  } else if (rtb_Switch2 < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S34>/Switch'
    rtb_Switch2 = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S34>/Switch2'

  // Sum: '<S31>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S31>/Delay Input2'
  //
  //  Block description for '<S31>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S31>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2 - localDW->DelayInput2_DSTATE;

  // Switch: '<S41>/Switch2' incorporates:
  //   Product: '<S31>/delta rise limit'
  //   SampleTimeMath: '<S31>/sample time'
  //
  //  About '<S31>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S41>/Switch2' incorporates:
  //   RelationalOperator: '<S41>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Switch2) {
    // Product: '<S31>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S31>/sample time'
    //
    //  About '<S31>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S41>/Switch' incorporates:
    //   RelationalOperator: '<S41>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Switch2) {
      // Switch: '<S41>/Switch2'
      rtb_Switch2 = rtb_UkYk1;
    }

    // End of Switch: '<S41>/Switch'
  }

  // End of Switch: '<S41>/Switch2'

  // Sum: '<S31>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S31>/Delay Input2'
  //
  //  Block description for '<S31>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S31>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S30>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S31>/Delay Input2'
  //
  //  Block description for '<S31>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_filteredInput;
}

//
// Function for Chart: '<S4>/Chart'
// function isTrue = checkRcCmds(cmds, params)
//
boolean_T fcsModel::fcsModel_checkRcCmds(const busRcInCmds
  *BusConversion_InsertedFor_Chart)
{
  boolean_T isTrue;

  // MATLAB Function 'checkRcCmds': '<S91>:7'
  // '<S91>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S91>:7:3' if(rcCmds.throttleCmd_nd <= pwmLowVal && ...
  // '<S91>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S91>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S91>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 987) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 987) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 987) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 987) {
          // '<S91>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S91>:7:8' else
          // '<S91>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S91>:7:8' else
        // '<S91>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S91>:7:8' else
      // '<S91>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S91>:7:8' else
    // '<S91>:7:9' isTrue = false;
    isTrue = false;
  }

  return isTrue;
}

// Model step function
void fcsModel::step()
{
  std::array<real_T, 9> conversionMatrix;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angRateC;
  std::array<busCtrlInputs, 3> rtb_ctrlInputsArray;
  busPidDebug rtb_BusCreator_o;
  real_T pCmd;
  real_T pCmd_tmp_tmp;
  real_T rCmd;
  real_T tCmd;
  real_T yCmd;
  boolean_T resetIntegrator;
  enumStateMachine state;

  // Chart: '<S4>/Chart' incorporates:
  //   Inport: '<Root>/rcCmdsIn'

  if (fcsModel_DW.temporalCounter_i1 < 8191U) {
    fcsModel_DW.temporalCounter_i1 = static_cast<uint16_T>
      (fcsModel_DW.temporalCounter_i1 + 1U);
  }

  // Gateway: rcInterpreter/Chart
  // During: rcInterpreter/Chart
  if (fcsModel_DW.is_active_c1_rcInterpreter == 0U) {
    // Entry: rcInterpreter/Chart
    fcsModel_DW.is_active_c1_rcInterpreter = 1U;

    // Entry Internal: rcInterpreter/Chart
    // Transition: '<S91>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S91>:1'
    // '<S91>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S91>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
      fcsModel_DW.durationCounter_1_d = 0;
    }

    // '<S91>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S91>:3'
      // '<S91>:10:1' sf_internal_predicateOutput = after(20, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 5000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S91>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S91>:1'
        // '<S91>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S91>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:1:4' resetIntegrator = true;

        // '<S91>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimits(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 987) {
        // Transition: '<S91>:12'
        fcsModel_DW.durationCounter_1_d = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S91>:11'
        // '<S91>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S91>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S91>:11:4' resetIntegrator = false;
      } else {
        // '<S91>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S91>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S91>:1'
      // '<S91>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 5 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 1250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S91>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S91>:3'
        // '<S91>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S91>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S91>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S91>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S91>:11'
      // '<S91>:17:1' sf_internal_predicateOutput = duration(rcCheckFlag == true, sec) >= 5; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1_d = 0;
      }

      if (fcsModel_DW.durationCounter_1_d >= 1250) {
        // Transition: '<S91>:17'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S91>:1'
        // '<S91>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S91>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S91>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S91>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
          fcsModel_DW.durationCounter_1_d = 0;
        }

        // '<S91>:11:4' resetIntegrator = false;
        resetIntegrator = false;
      }
      break;
    }
  }

  if (fcsModel_DW.rcCheckFlag) {
    fcsModel_DW.durationCounter_1++;
    fcsModel_DW.durationCounter_1_j++;
    fcsModel_DW.durationCounter_1_d++;
  } else {
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.durationCounter_1_j = 0;
    fcsModel_DW.durationCounter_1_d = 0;
  }

  // End of Chart: '<S4>/Chart'

  // MATLAB Function: '<S4>/Interpret RC In Cmds' incorporates:
  //   BusCreator generated from: '<S4>/Interpret RC In Cmds'
  //   Inport: '<Root>/rcCmdsIn'

  // Computes command and flight mode from the rc inputs
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S92>:1'
  // '<S92>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcCmds, rcParamsStruct); 
  // INTERPRETRCINPUTS_FUNCTION %Computes command and flight mode from the rc inputs 
  // 'interpretRcInputs_function:3' tCoeffs = rcParamsStruct.coeffs.throttle_nd; 
  // 'interpretRcInputs_function:4' tCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:5'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.throttleCmd_nd) ) ); 
  tCmd = std::fmin(2006.0, std::fmax(987.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.throttleCmd_nd)));

  // 'interpretRcInputs_function:6' rcOutCmds.throttleStick = -(tCoeffs(1)*tCmd^3 + tCoeffs(2)*tCmd^2 + ... 
  // 'interpretRcInputs_function:7'                           tCoeffs(3)*tCmd + tCoeffs(4))*rcParamsStruct.cmdLimits.zForce_N(2); 
  // 'interpretRcInputs_function:9' rCoeffs = rcParamsStruct.coeffs.roll_nd;
  // 'interpretRcInputs_function:10' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:11'        max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  rCmd = std::fmin(2006.0, std::fmax(987.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickXCmd_nd)));

  // 'interpretRcInputs_function:12' rcOutCmds.rollStick = (rCoeffs(1)*rCmd^3 + rCoeffs(2)*rCmd^2 + ... 
  // 'interpretRcInputs_function:13'                           rCoeffs(3)*rCmd + rCoeffs(4))*rcParamsStruct.cmdLimits.roll_rad(2); 
  // 'interpretRcInputs_function:15' pCoeffs = rcParamsStruct.coeffs.pitch_nd;
  //  Reverse the pitch cmd
  // 'interpretRcInputs_function:17' pCmd = double(rcInCmds.joystickYCmd_nd);
  // 'interpretRcInputs_function:18' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:19'        max( rcParamsStruct.pwmLimits(1),  pCmd) ); 
  pCmd = std::fmin(2006.0, std::fmax(987.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickYCmd_nd)));

  // 'interpretRcInputs_function:20' rcOutCmds.pitchStick = -(pCoeffs(1)*pCmd^3 + pCoeffs(2)*pCmd^2 + ... 
  // 'interpretRcInputs_function:21'                           pCoeffs(3)*pCmd + pCoeffs(4))*rcParamsStruct.cmdLimits.pitch_rad(2); 
  // 'interpretRcInputs_function:23' yCoeffs = rcParamsStruct.coeffs.yaw_nd;
  // 'interpretRcInputs_function:24' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:25'        max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  yCmd = std::fmin(2006.0, std::fmax(987.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickZCmd_nd)));

  // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
  //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
  //   Constant: '<S3>/Constant'
  //   Inport: '<Root>/stateEstimate'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'

  // 'interpretRcInputs_function:26' rcOutCmds.yawStick = (yCoeffs(1)*yCmd^3 + yCoeffs(2)*yCmd^2 + ... 
  // 'interpretRcInputs_function:27'                           yCoeffs(3)*yCmd + yCoeffs(4))*rcParamsStruct.cmdLimits.yawRate_radps(2); 
  // if(rcInCmds.rcSwitch1_nd < 1500)
  // 'interpretRcInputs_function:30' flightMode = enumFlightMode.STABILIZE;
  // end
  std::memcpy(&rtb_ctrlInputsArray[0],
              &fcsModel_ConstP.Constant_Value.attCtrlInputs.ctrlInputsArray[0],
              3U * sizeof(busCtrlInputs));

  // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S90>:1' 
  // '<S90>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
  //  This is a stop gap setup where we are only assuming that rate control
  //  is active and therefore not setting up attCtrlInputs for Euler angle
  //  control
  // '<S90>:1:6' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
  rtb_ctrlInputsArray[0].cmd = (((rCmd * rCmd * 1.0043823416568364E-5 +
    -2.2216571985540022E-9 * std::pow(rCmd, 3.0)) + -0.012613046004899127 * rCmd)
    + 3.8187541282968023) * 0.78539816339744828;

  // '<S90>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
  rtb_ctrlInputsArray[0].meas = fcsModel_U.stateEstimate.attitude_rad[0];

  // '<S90>:1:8' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
  rtb_ctrlInputsArray[1].cmd = -(((pCmd * pCmd * 1.0043823416568364E-5 +
    -2.2216571985540022E-9 * std::pow(pCmd, 3.0)) + -0.012613046004899127 * pCmd)
    + 3.8187541282968023) * 0.78539816339744828;

  // '<S90>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
  rtb_ctrlInputsArray[1].meas = fcsModel_U.stateEstimate.attitude_rad[1];

  // '<S90>:1:10' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
  rtb_ctrlInputsArray[2].cmd = (((yCmd * yCmd * 1.0043823416568364E-5 +
    -2.2216571985540022E-9 * std::pow(yCmd, 3.0)) + -0.012613046004899127 * yCmd)
    + 3.8187541282968023) * 3.1415926535897931;

  // '<S90>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
  rtb_ctrlInputsArray[2].meas = fcsModel_U.stateEstimate.attitude_rad[2];

  // '<S90>:1:13' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
  rtb_ctrlInputsArray[0].integratorReset = resetIntegrator;

  // '<S90>:1:14' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
  rtb_ctrlInputsArray[1].integratorReset = resetIntegrator;

  // '<S90>:1:15' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = resetIntegrator; 
  rtb_ctrlInputsArray[2].integratorReset = resetIntegrator;

  // Outputs for Iterator SubSystem: '<S9>/For Each Subsystem' incorporates:
  //   ForEach: '<S51>/For Each'

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S51>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S51>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[0],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[0].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[0],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[0].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S51>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, pCmd, rCmd, rtb_ctrlInputsArray[0].integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[0],
                        fcsModel_DW.CoreSubsys_p[0].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[0].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S51>/pidWithDebug'

  // Update for UnitDelay: '<S51>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[0].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S51>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug[0] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S51>/angRateCmds_radps'
  rtb_ImpAsg_InsertedFor_angRateC[0] = rCmd;

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S51>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S51>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[1],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[1].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[1],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[1].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S51>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, pCmd, rCmd, rtb_ctrlInputsArray[1].integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[1],
                        fcsModel_DW.CoreSubsys_p[1].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[1].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S51>/pidWithDebug'

  // Update for UnitDelay: '<S51>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[1].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S51>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug[1] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S51>/angRateCmds_radps'
  rtb_ImpAsg_InsertedFor_angRateC[1] = rCmd;

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S51>/ctrlInputs' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S51>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].cmd,
    &fcsModel_U.ctrlParams.attCtrlParams.cmdSignalConditioningParamsArray[2],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys_p[2].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].meas,
    &fcsModel_U.ctrlParams.attCtrlParams.measSignalConditioningParamsArray[2],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys_p[2].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S51>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, pCmd, rCmd, rtb_ctrlInputsArray[2].integratorReset,
                        &fcsModel_U.ctrlParams.attCtrlParams.ctrlParamsArray[2],
                        fcsModel_DW.CoreSubsys_p[2].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys_p[2].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S51>/pidWithDebug'

  // Update for UnitDelay: '<S51>/Unit Delay'
  fcsModel_DW.CoreSubsys_p[2].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S51>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug[2] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S51>/angRateCmds_radps'
  rtb_ImpAsg_InsertedFor_angRateC[2] = rCmd;

  // End of Outputs for SubSystem: '<S9>/For Each Subsystem'

  // MATLAB Function: '<S8>/EulerRates2BodyRates' incorporates:
  //   Inport: '<Root>/stateEstimate'
  //   SignalConversion generated from: '<S49>/ SFunction '
  //   Switch: '<S9>/Switch'

  // MATLAB Function 'EulerRates2BodyRates': '<S49>:1'
  // '<S49>:1:3' bodyRates_radps = eulerRates2bodyRates_function(taitBryanRates_radps,shipOrientation_rad); 
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
  rCmd = fcsModel_U.stateEstimate.attitude_rad[1];

  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  pCmd_tmp_tmp = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(pCmd_tmp_tmp - 1.5707963267948966) <= 0.004245395477824045) ||
      (std::abs(pCmd_tmp_tmp - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rCmd = -1.0;
      } else {
        rCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd *= pCmd_tmp_tmp - 0.004245395477824045;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rCmd = -1.0;
      } else {
        rCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd *= pCmd_tmp_tmp - 0.004245395477824045;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rCmd = -1.0;
      } else {
        rCmd = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rCmd *= pCmd_tmp_tmp + 0.004245395477824045;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  pCmd = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  yCmd = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  pCmd_tmp_tmp = std::cos(rCmd);
  conversionMatrix[0] = 1.0;
  conversionMatrix[3] = 0.0;
  conversionMatrix[6] = -std::sin(rCmd);
  conversionMatrix[1] = 0.0;
  conversionMatrix[4] = yCmd;
  conversionMatrix[7] = pCmd * pCmd_tmp_tmp;
  conversionMatrix[2] = 0.0;
  conversionMatrix[5] = -pCmd;
  conversionMatrix[8] = yCmd * pCmd_tmp_tmp;

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
  rCmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  pCmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  yCmd = rtb_ctrlInputsArray[2].cmd;
  for (int32_T ii{0}; ii < 3; ii++) {
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    pCmd_tmp_tmp = conversionMatrix[ii];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (pCmd_tmp_tmp <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      pCmd_tmp_tmp = 0.0;
    }

    conversionMatrix[ii] = pCmd_tmp_tmp;
    pCmd_tmp_tmp = conversionMatrix[ii + 3];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(pCmd_tmp_tmp) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      pCmd_tmp_tmp = 0.0;
    }

    conversionMatrix[ii + 3] = pCmd_tmp_tmp;
    pCmd_tmp_tmp = conversionMatrix[ii + 6];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(pCmd_tmp_tmp) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      pCmd_tmp_tmp = 0.0;
    }

    conversionMatrix[ii + 6] = pCmd_tmp_tmp;
    rtb_ImpAsg_InsertedFor_angRateC[ii] = (conversionMatrix[ii + 3] * pCmd +
      conversionMatrix[ii] * rCmd) + conversionMatrix[ii + 6] * yCmd;
  }

  // End of MATLAB Function: '<S8>/EulerRates2BodyRates'

  // BusCreator: '<S8>/Bus Creator' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'
  //   MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus'

  rtb_ctrlInputsArray[0].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[0].cmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_ctrlInputsArray[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_ctrlInputsArray[0].integratorReset = resetIntegrator;
  rtb_ctrlInputsArray[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator3' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'
  //   MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus'

  rtb_ctrlInputsArray[1].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[1].cmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  rtb_ctrlInputsArray[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_ctrlInputsArray[1].integratorReset = resetIntegrator;
  rtb_ctrlInputsArray[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator4' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'
  //   MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus'

  rtb_ctrlInputsArray[2].feedForwardCmd = 0.0;
  rtb_ctrlInputsArray[2].cmd = rtb_ImpAsg_InsertedFor_angRateC[2];
  rtb_ctrlInputsArray[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_ctrlInputsArray[2].integratorReset = resetIntegrator;
  rtb_ctrlInputsArray[2].trackingCtrlCmd = 0.0;

  // Outputs for Atomic SubSystem: '<S2>/Angular Rate Controller'
  // Outputs for Iterator SubSystem: '<S7>/For Each Subsystem' incorporates:
  //   ForEach: '<S10>/For Each'

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S10>/ctrlInputs' incorporates:
  //   BusCreator: '<S8>/Bus Creator1'
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S10>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[0],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[0].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[0].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[0],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[0].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, rCmd, pCmd, rtb_ctrlInputsArray[0].integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        0], fcsModel_DW.CoreSubsys[0].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[0].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S10>/pidWithDebug'

  // Update for UnitDelay: '<S10>/Unit Delay'
  fcsModel_DW.CoreSubsys[0].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S10>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug[0] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[0] = rCmd;

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S10>/ctrlInputs' incorporates:
  //   BusCreator: '<S8>/Bus Creator1'
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S10>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[1],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[1].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[1].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[1],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[1].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, rCmd, pCmd, rtb_ctrlInputsArray[1].integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        1], fcsModel_DW.CoreSubsys[1].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[1].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S10>/pidWithDebug'

  // Update for UnitDelay: '<S10>/Unit Delay'
  fcsModel_DW.CoreSubsys[1].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S10>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug[1] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[1] = rCmd;

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block'
  // ForEachSliceSelector generated from: '<S10>/ctrlInputs' incorporates:
  //   BusCreator: '<S8>/Bus Creator1'
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S10>/Unit Delay'

  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].cmd,
    &fcsModel_U.ctrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray[2],
    &rCmd, 0.004, &fcsModel_DW.CoreSubsys[2].SignalConditioningBlock);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block'

  // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block1'
  fcsMod_SignalConditioningBlock1(rtb_ctrlInputsArray[2].meas,
    &fcsModel_U.ctrlParams.angRateCtrlParams.measSignalConditioningParamsArray[2],
    &pCmd, 0.004, &fcsModel_DW.CoreSubsys[2].SignalConditioningBlock1);

  // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

  // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
  fcsModel_pidWithDebug(0.0, rCmd, pCmd, rtb_ctrlInputsArray[2].integratorReset,
                        &fcsModel_U.ctrlParams.angRateCtrlParams.ctrlParamsArray[
                        2], fcsModel_DW.CoreSubsys[2].UnitDelay_DSTATE, &rCmd,
                        &rtb_BusCreator_o, 0.004, &fcsModel_DW.CoreSubsys[2].
                        pidWithDebug);

  // End of Outputs for SubSystem: '<S10>/pidWithDebug'

  // Update for UnitDelay: '<S10>/Unit Delay'
  fcsModel_DW.CoreSubsys[2].UnitDelay_DSTATE = rCmd;

  // ForEachSliceAssignment generated from: '<S10>/pidDebug'
  fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug[2] = rtb_BusCreator_o;

  // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
  rtb_ImpAsg_InsertedFor_angAccel[2] = rCmd;

  // End of Outputs for SubSystem: '<S7>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'

  for (int32_T ii{0}; ii < 3; ii++) {
    rtb_ImpAsg_InsertedFor_angRateC[ii] = 0.0;
    rtb_ImpAsg_InsertedFor_angRateC[ii] += fcsModel_ConstP.Constant_Value_n[ii] *
      rtb_ImpAsg_InsertedFor_angAccel[0];
    rtb_ImpAsg_InsertedFor_angRateC[ii] += fcsModel_ConstP.Constant_Value_n[ii +
      3] * rtb_ImpAsg_InsertedFor_angAccel[1];
    rtb_ImpAsg_InsertedFor_angRateC[ii] += fcsModel_ConstP.Constant_Value_n[ii +
      6] * rtb_ImpAsg_InsertedFor_angAccel[2];
  }

  // End of Product: '<S2>/Matrix Multiply'

  // SignalConversion generated from: '<S1>/Matrix Multiply' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'

  pCmd_tmp_tmp = -(((tCmd * tCmd * 8.5458804152826638E-6 +
                     -1.9060550862738522E-9 * std::pow(tCmd, 3.0)) +
                    -0.011298232595401985 * tCmd) + 4.6650659714392448) * 46.0;
  rCmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  pCmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  yCmd = rtb_ImpAsg_InsertedFor_angRateC[2];

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);

  // BusCreator: '<Root>/Bus Creator' incorporates:
  //   Outport: '<Root>/fcsDebug'

  fcsModel_Y.fcsDebug.state = state;
  for (int32_T ii{0}; ii < 4; ii++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   Constant: '<S1>/Constant'

    tCmd = ((fcsModel_ConstP.Constant_Value_c[ii + 4] * rCmd +
             fcsModel_ConstP.Constant_Value_c[ii] * pCmd_tmp_tmp) +
            fcsModel_ConstP.Constant_Value_c[ii + 8] * pCmd) +
      fcsModel_ConstP.Constant_Value_c[ii + 12] * yCmd;

    // Saturate: '<S1>/Saturation'
    if (tCmd > 616850.27506808483) {
      tCmd = 616850.27506808483;
    } else if (tCmd < 0.0) {
      tCmd = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    tCmd = 9.5492965855137211 * std::sqrt(tCmd) - -0.92734095767679814 *
      fcsModel_DW.DiscreteTransferFcn_states[ii];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      // Outport: '<Root>/actuatorsCmds'
      fcsModel_Y.actuatorsCmds[ii] = -1.0;
    } else {
      // Outport: '<Root>/actuatorsCmds' incorporates:
      //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

      fcsModel_Y.actuatorsCmds[ii] = 0.036329521161600868 * tCmd +
        0.036329521161600868 * fcsModel_DW.DiscreteTransferFcn_states[ii];
    }

    // End of Switch: '<S1>/Switch'

    // Update for DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states[ii] = tCmd;
  }
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
