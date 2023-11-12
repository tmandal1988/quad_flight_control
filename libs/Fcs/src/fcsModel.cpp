//
// File: fcsModel.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.74
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Sun Nov  5 13:38:48 2023
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
#include <array>
#include <cstring>

// Named constants for Chart: '<S4>/Chart'
const uint8_T fcsModel_IN_ARM_MTRS{ 1U };

const uint8_T fcsModel_IN_INACTIVE{ 2U };

const uint8_T fcsModel_IN_INFLIGHT{ 3U };

const busOuterLoopToInnerLoop fcsModel_rtZbusOuterLoopToInnerLoop{
  {
    0.0                                // thrustCmd_N
  },                                   // outerLoopCmds

  {
    { {
        {
          0.0,                         // feedForwardCmd
          0.0,                         // cmd
          0.0,                         // meas
          false,                       // integratorReset
          0.0                          // trackingCtrlCmd
        }, {
          0.0,                         // feedForwardCmd
          0.0,                         // cmd
          0.0,                         // meas
          false,                       // integratorReset
          0.0                          // trackingCtrlCmd
        }, {
          0.0,                         // feedForwardCmd
          0.0,                         // cmd
          0.0,                         // meas
          false,                       // integratorReset
          0.0                          // trackingCtrlCmd
        } } }
    // ctrlInputsArray
  }                                    // attCtrlInputs
} ;                                    // busOuterLoopToInnerLoop ground

static void rate_scheduler(fcsModel::RT_MODEL_fcsModel_T *const fcsModel_M);

//
//         This function updates active task flag for each subrate.
//         The function is called at model base rate, hence the
//         generated code self-manages all its subrates.
//
static void rate_scheduler(fcsModel::RT_MODEL_fcsModel_T *const fcsModel_M)
{
  // Compute which subrates run during the next base time step.  Subrates
  //  are an integer multiple of the base rate counter.  Therefore, the subtask
  //  counter is reset when it reaches its limit (zero means run).

  (fcsModel_M->Timing.TaskCounters.TID[1])++;
  if ((fcsModel_M->Timing.TaskCounters.TID[1]) > 4) {// Sample time: [0.02s, 0.0s] 
    fcsModel_M->Timing.TaskCounters.TID[1] = 0;
  }
}

//
// Output and update for atomic system:
//    '<S13>/Discrete First Order Deriv Filter'
//    '<S57>/Discrete First Order Deriv Filter'
//    '<S103>/Discrete First Order Deriv Filter'
//
void fcsModel::f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
  rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
  rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW)
{
  real_T K;
  real_T normalizer;
  real_T num_tmp;

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
  K = 2.0 / rtp_sampleTime_s;

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
  normalizer = rtu_filterBandwidth_radps + K;

  // 'computeDiscreteTFNumAndDen_function:24' b0 = (B(1) + B(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:25' b1 = (B(1) - B(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:27' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:28' a1 = (A(1) - A(2)*K)/normalizer;
  // 'computeDiscreteTFNumAndDen_function:29' num = [b0, b1];
  num_tmp = rtu_filterBandwidth_radps * K;
  localDW->num[0] = num_tmp / normalizer;
  localDW->num[1] = (0.0 - num_tmp) / normalizer;

  // 'computeDiscreteTFNumAndDen_function:30' den = [a0, a1];
  localDW->den[0] = 1.0;
  localDW->den[1] = (rtu_filterBandwidth_radps - K) / normalizer;

  // DiscreteTransferFcn: '<S44>/Discrete Transfer Fcn'
  K = rtu_input - localDW->den[1] * localDW->DiscreteTransferFcn_states;
  *rty_filteredInputRate = localDW->num[0] * K + localDW->num[1] *
    localDW->DiscreteTransferFcn_states;

  // Update for DiscreteTransferFcn: '<S44>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states = K;
}

//
// Output and update for atomic system:
//    '<S10>/pidWithDebug'
//    '<S54>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
  *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd, busPidDebug
  *rty_pidDebug, real_T rtp_sampleTime_s, DW_pidWithDebug_fcsModel_T *localDW)
{
  real_T rtb_Product5;
  real_T rtb_Sum1_nj;
  real_T rtb_Sum_k;
  real_T rtb_Switch2_oz;
  real_T rtb_Switch2_p;
  real_T rtb_UkYk1_k;
  real_T rtb_UnitDelay_i;

  // Product: '<S45>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S45>/sample time'
  //
  //  About '<S45>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[1] * 0.004;

  // Sum: '<S13>/Sum'
  rtb_Sum_k = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S13>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_k,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S13>/Discrete First Order Deriv Filter'

  // Product: '<S13>/Product'
  rtb_Product5 *= rtu_pidParamBus->Kd;

  // Product: '<S13>/Product1'
  rtb_UnitDelay_i = rtb_Sum_k * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S13>/Discrete-Time Integrator'
  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = 0.0;
  }

  // Sum: '<S13>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rtb_Sum1_nj = ((rtu_feedForward + rtb_Product5) + rtb_UnitDelay_i) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S46>/Switch2' incorporates:
  //   RelationalOperator: '<S46>/LowerRelop1'
  //   RelationalOperator: '<S46>/UpperRelop'
  //   Switch: '<S46>/Switch'

  if (rtb_Sum1_nj > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_nj < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S46>/Switch'
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_oz = rtb_Sum1_nj;
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

  rtb_UkYk1_k = rtb_Switch2_oz - localDW->DelayInput2_DSTATE;

  // Switch: '<S48>/Switch2' incorporates:
  //   RelationalOperator: '<S48>/LowerRelop1'

  if (rtb_UkYk1_k <= rtb_Switch2_p) {
    // Product: '<S45>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S45>/sample time'
    //
    //  About '<S45>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S48>/Switch' incorporates:
    //   RelationalOperator: '<S48>/UpperRelop'

    if (rtb_UkYk1_k >= rtb_Switch2_p) {
      rtb_Switch2_p = rtb_UkYk1_k;
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

  *rty_ctrlCmd = rtb_Switch2_p + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S13>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_i;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5;

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
    rtu_pidParamBus->Kb) + rtb_Sum_k * rtu_pidParamBus->Ki) * 0.004;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S45>/Delay Input2'
  //
  //  Block description for '<S45>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S13>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_oz;

  // Update for UnitDelay: '<S13>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_nj;
}

//
// Output and update for atomic system:
//    '<S29>/Compute Natural Frequency'
//    '<S30>/Compute Natural Frequency'
//    '<S14>/Compute Natural Frequency'
//    '<S15>/Compute Natural Frequency'
//    '<S73>/Compute Natural Frequency'
//    '<S74>/Compute Natural Frequency'
//    '<S58>/Compute Natural Frequency'
//    '<S59>/Compute Natural Frequency'
//    '<S119>/Compute Natural Frequency'
//    '<S120>/Compute Natural Frequency'
//    ...
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
//    '<S29>/Compute Numerator And Denominator'
//    '<S14>/Compute Numerator And Denominator'
//    '<S73>/Compute Numerator And Denominator'
//    '<S58>/Compute Numerator And Denominator'
//    '<S119>/Compute Numerator And Denominator'
//    '<S104>/Compute Numerator And Denominator'
//
void fcsModel::ComputeNumeratorAndDenominator(real_T rtu_naturalFrequency_radps,
  real_T rtu_dampingRatio_nd, real_T rty_rateNum[3], real_T rty_accelNum[3],
  real_T rty_den[3], real_T rtp_sampleTime_s)
{
  real_T B_idx_1_tmp;
  real_T K;
  real_T normalizer_tmp;
  real_T normalizer_tmp_0;
  real_T normalizer_tmp_1;

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
  K = 2.0 / rtp_sampleTime_s;

  // 'computeSecondOrderDerivFilterNumAndDen_function:16' A0 = naturalFrequency_radps^2; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:17' A1 = 2*dampingRatio_nd*naturalFrequency_radps; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:18' A2 = 1;
  // 'computeSecondOrderDerivFilterNumAndDen_function:20' B0 = 0;
  // 'computeSecondOrderDerivFilterNumAndDen_function:22' B1 = naturalFrequency_radps^2; 
  // 'computeSecondOrderDerivFilterNumAndDen_function:23' B2 = 0;
  //  compute the rate transfer function numerator and the denominator
  // 'computeSecondOrderDerivFilterNumAndDen_function:26' [rateNum, den] = computeDiscreteTFNumAndDen_function([B0, B1, B2], [A0, A1, A2], K); 
  B_idx_1_tmp = rtu_naturalFrequency_radps * rtu_naturalFrequency_radps;

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
  normalizer_tmp = K * K;
  normalizer_tmp_0 = 2.0 * rtu_dampingRatio_nd * rtu_naturalFrequency_radps * K;
  normalizer_tmp_1 = (normalizer_tmp_0 + B_idx_1_tmp) + normalizer_tmp;

  // 'computeDiscreteTFNumAndDen_function:35' b0 = (B(1) + B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:36' b1 = (2*B(1) - 2*B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:37' b2 =  (B(1) - B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:39' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:40' a1 = (2*A(1) - 2*A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:41' a2 = (A(1) - A(2)*K + A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:42' num = [b0, b1, b2];
  K *= B_idx_1_tmp;
  rty_rateNum[0] = K / normalizer_tmp_1;
  rty_rateNum[1] = 0.0 / normalizer_tmp_1;
  rty_rateNum[2] = (0.0 - K) / normalizer_tmp_1;

  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
  rty_den[0] = 1.0;
  rty_den[1] = (2.0 * B_idx_1_tmp - normalizer_tmp * 2.0) / normalizer_tmp_1;
  rty_den[2] = ((B_idx_1_tmp - normalizer_tmp_0) + normalizer_tmp) /
    normalizer_tmp_1;

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
  K = normalizer_tmp * B_idx_1_tmp / normalizer_tmp_1;
  rty_accelNum[0] = K;
  rty_accelNum[1] = (0.0 - 2.0 * B_idx_1_tmp * normalizer_tmp) /
    normalizer_tmp_1;
  rty_accelNum[2] = K;

  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
}

//
// System initialize for atomic system:
//    '<S30>/Compute Filter Numerator And Denominator'
//    '<S15>/Compute Filter Numerator And Denominator'
//    '<S74>/Compute Filter Numerator And Denominator'
//    '<S59>/Compute Filter Numerator And Denominator'
//    '<S120>/Compute Filter Numerator And Denominator'
//    '<S105>/Compute Filter Numerator And Denominator'
//
void fcsModel::ComputeFilterNumeratorAndD_Init(real_T rty_num[3], real_T
  rty_den[3])
{
  rty_num[0] = 0.0;
  rty_den[0] = 0.0;
  rty_num[1] = 0.0;
  rty_den[1] = 0.0;
  rty_num[2] = 0.0;
  rty_den[2] = 0.0;
}

//
// Output and update for atomic system:
//    '<S30>/Compute Filter Numerator And Denominator'
//    '<S15>/Compute Filter Numerator And Denominator'
//    '<S74>/Compute Filter Numerator And Denominator'
//    '<S59>/Compute Filter Numerator And Denominator'
//    '<S120>/Compute Filter Numerator And Denominator'
//    '<S105>/Compute Filter Numerator And Denominator'
//
void fcsModel::ComputeFilterNumeratorAndDenomi(real_T rtu_naturalFrequency_radps,
  real_T rtu_dampingRatio_nd, real_T rty_num[3], real_T rty_den[3], real_T
  rtp_sampleTime_s)
{
  real_T B0;
  real_T B0_tmp;
  real_T B_idx_0;
  real_T K;
  real_T tmp;

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
  B0 = rtu_naturalFrequency_radps * rtu_naturalFrequency_radps;

  // 'computeSecondOrderFilterNumAndDen_function:15' B1 = 0;
  // 'computeSecondOrderFilterNumAndDen_function:16' B2 = 0;
  // 'computeSecondOrderFilterNumAndDen_function:18' A0 = B0;
  // 'computeSecondOrderFilterNumAndDen_function:19' A1 = 2*dampingRatio_nd*naturalFrequency_radps; 
  // 'computeSecondOrderFilterNumAndDen_function:20' A2 = 1;
  // 'computeSecondOrderFilterNumAndDen_function:22' K = 2/sampleTime_s;
  K = 2.0 / rtp_sampleTime_s;

  // 'computeSecondOrderFilterNumAndDen_function:24' [num, den] = computeDiscreteTFNumAndDen_function([B0, B1, B2], [A0, A1, A2], K); 
  B_idx_0 = B0;

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
  B0_tmp = K * K;
  K *= 2.0 * rtu_dampingRatio_nd * rtu_naturalFrequency_radps;
  B0 = (K + B0) + B0_tmp;

  // 'computeDiscreteTFNumAndDen_function:35' b0 = (B(1) + B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:36' b1 = (2*B(1) - 2*B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:37' b2 =  (B(1) - B(2)*K + B(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:39' a0 = 1;
  // 'computeDiscreteTFNumAndDen_function:40' a1 = (2*A(1) - 2*A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:41' a2 = (A(1) - A(2)*K + A(3)*K^2)/normalizer; 
  // 'computeDiscreteTFNumAndDen_function:42' num = [b0, b1, b2];
  tmp = B_idx_0 / B0;
  rty_num[0] = tmp;
  rty_num[1] = 2.0 * B_idx_0 / B0;
  rty_num[2] = tmp;

  // 'computeDiscreteTFNumAndDen_function:43' den = [a0, a1, a2];
  rty_den[0] = 1.0;
  rty_den[1] = (2.0 * B_idx_0 - B0_tmp * 2.0) / B0;
  rty_den[2] = ((B_idx_0 - K) + B0_tmp) / B0;
}

//
// System initialize for atomic system:
//    '<S10>/Signal Conditioning Block1'
//    '<S10>/Signal Conditioning Block'
//    '<S54>/Signal Conditioning Block1'
//    '<S54>/Signal Conditioning Block'
//
void fcsModel::f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S30>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S10>/Signal Conditioning Block1'
//    '<S10>/Signal Conditioning Block'
//    '<S54>/Signal Conditioning Block1'
//    '<S54>/Signal Conditioning Block'
//
void fcsModel::fcsMod_SignalConditioningBlock1(real_T rtu_input, const
  busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
  rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW)
{
  std::array<real_T, 3> rtb_accelNum;
  std::array<real_T, 3> rtb_den;
  std::array<real_T, 3> rtb_rateNum;
  real_T rtb_DiscreteTransferFcn_j;
  real_T rtb_Switch2_h;

  // MATLAB Function: '<S29>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S29>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S30>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S30>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S30>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_j = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S34>/Switch2' incorporates:
  //   RelationalOperator: '<S34>/LowerRelop1'
  //   RelationalOperator: '<S34>/UpperRelop'
  //   Switch: '<S34>/Switch'

  if (rtb_DiscreteTransferFcn_j > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_j < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S34>/Switch'
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[0];
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

  rtb_DiscreteTransferFcn_j -= localDW->DelayInput2_DSTATE;

  // Switch: '<S41>/Switch2' incorporates:
  //   Product: '<S31>/delta rise limit'
  //   SampleTimeMath: '<S31>/sample time'
  //
  //  About '<S31>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_h = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S41>/Switch2' incorporates:
  //   RelationalOperator: '<S41>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_j <= rtb_Switch2_h) {
    // Product: '<S31>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S31>/sample time'
    //
    //  About '<S31>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_h = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S41>/Switch' incorporates:
    //   RelationalOperator: '<S41>/UpperRelop'

    if (rtb_DiscreteTransferFcn_j >= rtb_Switch2_h) {
      // Switch: '<S41>/Switch2'
      rtb_Switch2_h = rtb_DiscreteTransferFcn_j;
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

  *rty_filteredInput = rtb_Switch2_h + localDW->DelayInput2_DSTATE;

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
// System initialize for atomic system:
//    '<S100>/Signal Conditioning Block1'
//    '<S100>/Signal Conditioning Block'
//
void fcsModel::SignalConditioningBlock1_j_Init(DW_SignalConditioningBlock1_k_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S120>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S100>/Signal Conditioning Block1'
//    '<S100>/Signal Conditioning Block'
//
void fcsModel::fcsM_SignalConditioningBlock1_d(real_T rtu_input, const
  busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
  rtp_sampleTime_s, DW_SignalConditioningBlock1_k_T *localDW)
{
  std::array<real_T, 3> rtb_accelNum;
  std::array<real_T, 3> rtb_den;
  std::array<real_T, 3> rtb_rateNum;
  real_T rtb_DiscreteTransferFcn_m;

  // MATLAB Function: '<S119>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, rty_filteredInput);

  // MATLAB Function: '<S119>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(*rty_filteredInput,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S120>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, rty_filteredInput);

  // MATLAB Function: '<S120>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(*rty_filteredInput,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S120>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_m = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S124>/Switch2' incorporates:
  //   RelationalOperator: '<S124>/LowerRelop1'
  //   RelationalOperator: '<S124>/UpperRelop'
  //   Switch: '<S124>/Switch'

  if (rtb_DiscreteTransferFcn_m > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_m = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_m < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S124>/Switch'
    rtb_DiscreteTransferFcn_m = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S124>/Switch2'

  // Sum: '<S121>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S121>/Delay Input2'
  //
  //  Block description for '<S121>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S121>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_m -= localDW->DelayInput2_DSTATE;

  // Product: '<S121>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S121>/sample time'
  //
  //  About '<S121>/sample time':
  //   y = K where K = ( w * Ts )

  *rty_filteredInput = rtu_params->filteredInputRateLimits[1] * 0.02;

  // Switch: '<S131>/Switch2' incorporates:
  //   RelationalOperator: '<S131>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_m <= *rty_filteredInput) {
    real_T rtb_deltafalllimit_i2;

    // Product: '<S121>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S121>/sample time'
    //
    //  About '<S121>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_deltafalllimit_i2 = rtu_params->filteredInputRateLimits[0] * 0.02;

    // Switch: '<S131>/Switch' incorporates:
    //   RelationalOperator: '<S131>/UpperRelop'

    if (rtb_DiscreteTransferFcn_m < rtb_deltafalllimit_i2) {
      *rty_filteredInput = rtb_deltafalllimit_i2;
    } else {
      *rty_filteredInput = rtb_DiscreteTransferFcn_m;
    }

    // End of Switch: '<S131>/Switch'
  }

  // End of Switch: '<S131>/Switch2'

  // Sum: '<S121>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S121>/Delay Input2'
  //
  //  Block description for '<S121>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S121>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput += localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S120>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S121>/Delay Input2'
  //
  //  Block description for '<S121>/Delay Input2':
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

  // MATLAB Function 'checkRcCmds': '<S139>:7'
  // '<S139>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S139>:7:3' if(rcCmds.throttleCmd_nd <= pwmLowVal && ...
  // '<S139>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S139>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S139>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 1000) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 1000) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 1000) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 1000) {
          // '<S139>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S139>:7:8' else
          // '<S139>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S139>:7:8' else
        // '<S139>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S139>:7:8' else
      // '<S139>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S139>:7:8' else
    // '<S139>:7:9' isTrue = false;
    isTrue = false;
  }

  return isTrue;
}

// Model step function
void fcsModel::step()
{
  // local scratch DWork variables
  int32_T ForEach_itr_h;
  int32_T ForEach_itr_p;
  int32_T ForEach_itr;
  std::array<real_T, 4> DiscreteTransferFcn_tmp;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angRateC;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtCmd_;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtMeas;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDebug;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_velCtrlO;
  std::array<real_T, 3> rtb_Switch_f;
  std::array<busCtrlInputs, 3> rtb_VectorConcatenate;
  std::array<real_T, 9> rtb_rFepTpNed;
  busCtrlInputs rtb_BusAssignment1;
  busCtrlInputs rtb_BusAssignment2;
  busCtrlInputs rtb_BusAssignment3_d;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_2;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_5;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_p;
  busPidDebug rtb_BusCreator_o;
  real_T c_psi;
  real_T rlim;
  real_T rtb_BusCreator_b_frcCmd_N;
  real_T rtb_BusCreator_b_velCtrlDebug_0;
  real_T rtb_BusCreator_b_velCtrlDebug_1;
  real_T rtb_BusCreator_b_velCtrlDebug_3;
  real_T rtb_BusCreator_b_velCtrlDebug_4;
  real_T rtb_BusCreator_b_velCtrlDebug_c;
  real_T rtb_BusCreator_b_velCtrlDebug_m;
  real_T rtb_Product3;
  real_T rtb_outDebug;
  real_T s_psi;
  real_T tlim;
  boolean_T resetIntegrator;
  enumFlightMode flightMode;
  enumStateMachine state;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // MATLAB Function: '<S93>/fepToNedRotationMatrix' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'Outer Loop Controller/Assemble Vel Ctrl Inputs/fepToNedRotationMatrix': '<S98>:1' 
    // '<S98>:1:2' c_psi = cos(yaw_rad);
    c_psi = std::cos(fcsModel_U.stateEstimate.attitude_rad[2]);

    // '<S98>:1:3' s_psi = sin(yaw_rad);
    s_psi = std::sin(fcsModel_U.stateEstimate.attitude_rad[2]);

    // '<S98>:1:5' rFepTpNed = [c_psi,  -s_psi, 0;
    // '<S98>:1:6' 			s_psi, c_psi, 0;
    // '<S98>:1:7' 			0, 0, 1];
    rtb_rFepTpNed[0] = c_psi;
    rtb_rFepTpNed[3] = -s_psi;
    rtb_rFepTpNed[6] = 0.0;
    rtb_rFepTpNed[1] = s_psi;
    rtb_rFepTpNed[4] = c_psi;
    rtb_rFepTpNed[7] = 0.0;
    rtb_rFepTpNed[2] = 0.0;
    rtb_rFepTpNed[5] = 0.0;
    rtb_rFepTpNed[8] = 1.0;
  }

  // Chart: '<S4>/Chart' incorporates:
  //   Inport: '<Root>/rcCmdsIn'

  if (fcsModel_DW.temporalCounter_i1 < 16383U) {
    fcsModel_DW.temporalCounter_i1 = static_cast<uint16_T>
      (fcsModel_DW.temporalCounter_i1 + 1U);
  }

  // Gateway: rcInterpreter/Chart
  // During: rcInterpreter/Chart
  if (fcsModel_DW.is_active_c1_rcInterpreter == 0U) {
    // Entry: rcInterpreter/Chart
    fcsModel_DW.is_active_c1_rcInterpreter = 1U;

    // Entry Internal: rcInterpreter/Chart
    // Transition: '<S139>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S139>:1'
    // '<S139>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S139>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
    }

    // '<S139>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S139>:3'
      // '<S139>:10:1' sf_internal_predicateOutput = after(60, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 15000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S139>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S139>:1'
        // '<S139>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S139>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S139>:1:4' resetIntegrator = true;

        // '<S139>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimits(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 1000) {
        // Transition: '<S139>:12'
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S139>:11'
        // '<S139>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S139>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S139>:11:4' resetIntegrator = false;
      } else {
        // '<S139>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S139>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S139>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S139>:1'
      // '<S139>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 5 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 1250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S139>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S139>:3'
        // '<S139>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S139>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S139>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S139>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S139>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S139>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S139>:11'
      // '<S139>:20:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd <= paramsStruct.pwmLimits(1); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1000) {
        // Transition: '<S139>:20'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S139>:3'
        // '<S139>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S139>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S139>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S139>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S139>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S139>:11:4' resetIntegrator = false;
        resetIntegrator = false;
      }
      break;
    }
  }

  if (fcsModel_DW.rcCheckFlag) {
    fcsModel_DW.durationCounter_1++;
    fcsModel_DW.durationCounter_1_j++;
  } else {
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.durationCounter_1_j = 0;
  }

  // End of Chart: '<S4>/Chart'

  // MATLAB Function: '<S4>/Interpret RC In Cmds' incorporates:
  //   BusCreator generated from: '<S4>/Interpret RC In Cmds'
  //   Inport: '<Root>/rcCmdsIn'

  // Computes command and flight mode from the rc inputs
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S140>:1'
  // '<S140>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcCmds, expo, rcParamsStruct); 
  // INTERPRETRCINPUTS_FUNCTION
  // Computes command and flight mode from the rc inputs
  // Select mode
  // 'interpretRcInputs_function:6' if(rcInCmds.rcSwitch1_nd < 1100)
  if (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1100) {
    // 'interpretRcInputs_function:7' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:8' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -60.0;

    // 'interpretRcInputs_function:9' rlim = rcParamsStruct.cmdLimits.rollRate_radps(2); 
    rlim = 1.0471975511965976;

    // 'interpretRcInputs_function:10' plim = rcParamsStruct.cmdLimits.pitchRate_radps(2); 
    s_psi = 1.0471975511965976;

    // 'interpretRcInputs_function:11' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    c_psi = 1.0471975511965976;
  } else if ((fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1100) &&
             (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1800)) {
    // 'interpretRcInputs_function:13' elseif (rcInCmds.rcSwitch1_nd >= 1100 && rcInCmds.rcSwitch1_nd < 1800) 
    // 'interpretRcInputs_function:14' flightMode = enumFlightMode.VEL_CONTROL;
    flightMode = enumFlightMode::VEL_CONTROL;

    // 'interpretRcInputs_function:15' tlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    tlim = -1.0;

    // 'interpretRcInputs_function:16' rlim = rcParamsStruct.cmdLimits.vy_mps(2); 
    rlim = 2.5;

    // 'interpretRcInputs_function:17' plim = rcParamsStruct.cmdLimits.vx_mps(2); 
    s_psi = 2.5;

    // 'interpretRcInputs_function:18' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    c_psi = 1.0471975511965976;
  } else {
    // 'interpretRcInputs_function:19' else
    // 'interpretRcInputs_function:20' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:21' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -60.0;

    // 'interpretRcInputs_function:22' rlim = rcParamsStruct.cmdLimits.rollRate_radps(2); 
    rlim = 1.0471975511965976;

    // 'interpretRcInputs_function:23' plim = rcParamsStruct.cmdLimits.pitchRate_radps(2); 
    s_psi = 1.0471975511965976;

    // 'interpretRcInputs_function:24' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    c_psi = 1.0471975511965976;
  }

  // 'interpretRcInputs_function:27' if expo
  // 'interpretRcInputs_function:60' else
  // 'interpretRcInputs_function:61' tCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:62'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.throttleCmd_nd) ) ); 
  // 'interpretRcInputs_function:63' tCmd_unitRange = -1 + tCmd/1000;
  // 'interpretRcInputs_function:65' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:66'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  // 'interpretRcInputs_function:67' rCmd_unitRange = -3 + rCmd/500;
  // 'interpretRcInputs_function:69' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:70'         max( rcParamsStruct.pwmLimits(1),  double(rcInCmds.joystickYCmd_nd) ) ); 
  // 'interpretRcInputs_function:71' pCmd_unitRange = -(-3 + pCmd/500);
  // 'interpretRcInputs_function:73' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:74'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  // 'interpretRcInputs_function:75' yCmd_unitRange = -3 + yCmd/500;
  // 'interpretRcInputs_function:77' rcOutCmds.throttleStick = tCmd_unitRange*tlim; 
  fcsModel_DW.rcOutCmds.throttleStick = (std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.throttleCmd_nd))) / 1000.0 + -1.0) *
    tlim;

  // 'interpretRcInputs_function:78' rcOutCmds.rollStick = rCmd_unitRange*rlim;
  fcsModel_DW.rcOutCmds.rollStick = (std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickXCmd_nd))) / 500.0 + -3.0) *
    rlim;

  // 'interpretRcInputs_function:79' rcOutCmds.pitchStick = pCmd_unitRange*plim; 
  fcsModel_DW.rcOutCmds.pitchStick = -(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickYCmd_nd))) / 500.0 + -3.0) *
    s_psi;

  // 'interpretRcInputs_function:80' rcOutCmds.yawStick = yCmd_unitRange*ylim;
  fcsModel_DW.rcOutCmds.yawStick = (std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickZCmd_nd))) / 500.0 + -3.0) *
    1.0471975511965976;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    busOuterLoopToInnerLoop rtb_outBus;

    // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
    //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
    //   Constant: '<S3>/Constant'
    //   Inport: '<Root>/stateEstimate'

    rtb_outBus = fcsModel_rtZbusOuterLoopToInnerLoop;

    // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S96>:1' 
    // '<S96>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
    rtb_outBus.outerLoopCmds.thrustCmd_N = fcsModel_DW.rcOutCmds.throttleStick;

    // '<S96>:1:3' outDebug = throttleCmd_N;
    rtb_outDebug = fcsModel_DW.rcOutCmds.throttleStick;

    //  This is a stop gap setup where we are only assuming that rate control
    //  is active and therefore not setting up attCtrlInputs for Euler angle
    //  control
    // '<S96>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].cmd =
      fcsModel_DW.rcOutCmds.rollStick;

    // '<S96>:1:8' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].meas =
      fcsModel_U.stateEstimate.attitude_rad[0];

    // '<S96>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].cmd =
      fcsModel_DW.rcOutCmds.pitchStick;

    // '<S96>:1:10' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].meas =
      fcsModel_U.stateEstimate.attitude_rad[1];

    // '<S96>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].cmd =
      fcsModel_DW.rcOutCmds.yawStick;

    // '<S96>:1:12' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].meas =
      fcsModel_U.stateEstimate.attitude_rad[2];

    // '<S96>:1:14' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].integratorReset =
      resetIntegrator;

    // '<S96>:1:15' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].integratorReset =
      resetIntegrator;

    // '<S96>:1:16' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].integratorReset =
      resetIntegrator;

    // Switch: '<S93>/Switch' incorporates:
    //   Constant: '<S93>/Constant'
    //   Constant: '<S97>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S97>/Compare'

    if (flightMode == enumFlightMode::VEL_CONTROL) {
      // Product: '<S93>/Matrix Multiply' incorporates:
      //   SignalConversion generated from: '<S93>/Vector Concatenate1'

      for (int32_T ii{0}; ii < 3; ii++) {
        rtb_Switch_f[ii] = 0.0;
        rtb_Switch_f[ii] += rtb_rFepTpNed[ii] * fcsModel_DW.rcOutCmds.pitchStick;
        rtb_Switch_f[ii] += rtb_rFepTpNed[ii + 3] *
          fcsModel_DW.rcOutCmds.rollStick;
        rtb_Switch_f[ii] += rtb_rFepTpNed[ii + 6] *
          fcsModel_DW.rcOutCmds.throttleStick;
      }

      // End of Product: '<S93>/Matrix Multiply'
    } else {
      rtb_Switch_f[0] = 0.0;
      rtb_Switch_f[1] = 0.0;
      rtb_Switch_f[2] = 0.0;
    }

    // End of Switch: '<S93>/Switch'

    // Concatenate: '<S93>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S93>/Bus Assignment' incorporates:
    //   Concatenate: '<S93>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[0].cmd = rtb_Switch_f[0];
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedVel_mps[0];
    rtb_VectorConcatenate[0].integratorReset = resetIntegrator;

    // Concatenate: '<S93>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S93>/Bus Assignment1' incorporates:
    //   Concatenate: '<S93>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[1].cmd = rtb_Switch_f[1];
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedVel_mps[1];
    rtb_VectorConcatenate[1].integratorReset = resetIntegrator;

    // Concatenate: '<S93>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S93>/Bus Assignment2' incorporates:
    //   Concatenate: '<S93>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[2].cmd = rtb_Switch_f[2];
    rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.nedVel_mps[2];
    rtb_VectorConcatenate[2].integratorReset = resetIntegrator;

    // Outputs for Iterator SubSystem: '<S95>/For Each Subsystem' incorporates:
    //   ForEach: '<S100>/For Each'

    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      real_T rtb_Product1;
      real_T rtb_Sum1_b;
      real_T rtb_Switch2;
      real_T tCmd_unitRange;

      // Outputs for Atomic SubSystem: '<S100>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S100>/ctrlInputs' incorporates:
      //   BusAssignment: '<S93>/Bus Assignment3'
      //   Concatenate: '<S99>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'

      fcsM_SignalConditioningBlock1_d(rtb_VectorConcatenate[ForEach_itr].cmd,
        &fcsModel_U.ctrlParams.velCtrlParams.cmdSignalConditioningParamsArray[ForEach_itr],
        &rlim, 0.02, &fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S100>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S100>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_d(rtb_VectorConcatenate[ForEach_itr].meas,
        &fcsModel_U.ctrlParams.velCtrlParams.measSignalConditioningParamsArray[ForEach_itr],
        &s_psi, 0.02, &fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S100>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S100>/pidWithDebug'
      // Product: '<S135>/delta rise limit' incorporates:
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'
      //   SampleTimeMath: '<S135>/sample time'
      //
      //  About '<S135>/sample time':
      //   y = K where K = ( w * Ts )

      rtb_Product1 =
        fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
        outputRateLimits[1] * 0.02;

      // Sum: '<S103>/Sum'
      tlim = rlim - s_psi;

      // Outputs for Atomic SubSystem: '<S103>/Discrete First Order Deriv Filter' 
      // ForEachSliceSelector generated from: '<S100>/pidParams' incorporates:
      //   Inport: '<Root>/ctrlParams'

      f_DiscreteFirstOrderDerivFilter(tlim,
        fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
        filterBandwidth_radps, &rtb_Product3, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteFirstOrderDerivFilter);

      // End of Outputs for SubSystem: '<S103>/Discrete First Order Deriv Filter' 

      // Product: '<S103>/Product' incorporates:
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'

      c_psi = rtb_Product3 *
        fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].Kd;

      // Product: '<S103>/Product1' incorporates:
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'

      tCmd_unitRange = tlim *
        fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].Kp;

      // DiscreteIntegrator: '<S103>/Discrete-Time Integrator' incorporates:
      //   BusAssignment: '<S93>/Bus Assignment3'
      //   Concatenate: '<S99>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S100>/ctrlInputs'

      if (rtb_VectorConcatenate[ForEach_itr].integratorReset ||
          (fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_PrevRese
           != 0)) {
        fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_DSTATE =
          0.0;
      }

      // Sum: '<S103>/Sum1' incorporates:
      //   DiscreteIntegrator: '<S103>/Discrete-Time Integrator'

      rtb_Sum1_b = (c_psi + tCmd_unitRange) +
        fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_DSTATE;

      // Switch: '<S136>/Switch2' incorporates:
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'
      //   RelationalOperator: '<S136>/LowerRelop1'
      //   RelationalOperator: '<S136>/UpperRelop'
      //   Switch: '<S136>/Switch'

      if (rtb_Sum1_b >
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
          outputLimits[1]) {
        rtb_Switch2 =
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
          outputLimits[1];
      } else if (rtb_Sum1_b <
                 fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr]
                 .outputLimits[0]) {
        // Switch: '<S136>/Switch'
        rtb_Switch2 =
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
          outputLimits[0];
      } else {
        rtb_Switch2 = rtb_Sum1_b;
      }

      // End of Switch: '<S136>/Switch2'

      // Sum: '<S135>/Difference Inputs1' incorporates:
      //   UnitDelay: '<S135>/Delay Input2'
      //
      //  Block description for '<S135>/Difference Inputs1':
      //
      //   Add in CPU
      //
      //  Block description for '<S135>/Delay Input2':
      //
      //   Store in Global RAM

      rtb_Product3 = rtb_Switch2 - fcsModel_DW.CoreSubsys_i[ForEach_itr].
        DelayInput2_DSTATE;

      // Switch: '<S138>/Switch2' incorporates:
      //   RelationalOperator: '<S138>/LowerRelop1'

      if (rtb_Product3 <= rtb_Product1) {
        // Product: '<S135>/delta fall limit' incorporates:
        //   ForEachSliceSelector generated from: '<S100>/pidParams'
        //   Inport: '<Root>/ctrlParams'
        //   SampleTimeMath: '<S135>/sample time'
        //
        //  About '<S135>/sample time':
        //   y = K where K = ( w * Ts )

        rtb_Product1 =
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].
          outputRateLimits[0] * 0.02;

        // Switch: '<S138>/Switch' incorporates:
        //   RelationalOperator: '<S138>/UpperRelop'

        if (rtb_Product3 >= rtb_Product1) {
          rtb_Product1 = rtb_Product3;
        }

        // End of Switch: '<S138>/Switch'
      }

      // End of Switch: '<S138>/Switch2'

      // Sum: '<S135>/Difference Inputs2' incorporates:
      //   UnitDelay: '<S135>/Delay Input2'
      //
      //  Block description for '<S135>/Difference Inputs2':
      //
      //   Add in CPU
      //
      //  Block description for '<S135>/Delay Input2':
      //
      //   Store in Global RAM

      rtb_Product1 += fcsModel_DW.CoreSubsys_i[ForEach_itr].DelayInput2_DSTATE;

      // ForEachSliceAssignment generated from: '<S100>/pidDebug' incorporates:
      //   BusCreator: '<S103>/Bus Creator'
      //   DiscreteIntegrator: '<S103>/Discrete-Time Integrator'

      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr].integralOutput =
        fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_DSTATE;

      // Product: '<S103>/Product3' incorporates:
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'

      rtb_Product3 = tlim *
        fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].Ki;

      // Update for DiscreteIntegrator: '<S103>/Discrete-Time Integrator' incorporates:
      //   BusAssignment: '<S93>/Bus Assignment3'
      //   Concatenate: '<S99>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S100>/pidParams'
      //   Inport: '<Root>/ctrlParams'
      //   Product: '<S103>/Product2'
      //   Product: '<S103>/Product5'
      //   Sum: '<S103>/Sum2'
      //   Sum: '<S103>/Sum3'
      //   Sum: '<S103>/Sum4'
      //   Sum: '<S103>/Sum5'
      //   UnitDelay: '<S100>/Unit Delay'
      //   UnitDelay: '<S103>/Unit Delay'
      //   UnitDelay: '<S103>/Unit Delay1'

      fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_DSTATE +=
        (((fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE -
           fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE_c) *
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].Kt +
          (fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE_c -
           fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay1_DSTATE) *
          fcsModel_U.ctrlParams.velCtrlParams.ctrlParamsArray[ForEach_itr].Kb) +
         rtb_Product3) * 0.02;
      fcsModel_DW.CoreSubsys_i[ForEach_itr].DiscreteTimeIntegrator_PrevRese =
        static_cast<int8_T>(rtb_VectorConcatenate[ForEach_itr].integratorReset);

      // Update for UnitDelay: '<S135>/Delay Input2'
      //
      //  Block description for '<S135>/Delay Input2':
      //
      //   Store in Global RAM

      fcsModel_DW.CoreSubsys_i[ForEach_itr].DelayInput2_DSTATE = rtb_Product1;

      // Update for UnitDelay: '<S103>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE_c = rtb_Switch2;

      // Update for UnitDelay: '<S103>/Unit Delay1'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay1_DSTATE = rtb_Sum1_b;

      // End of Outputs for SubSystem: '<S100>/pidWithDebug'

      // Update for UnitDelay: '<S100>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE = rtb_Product1;

      // ForEachSliceAssignment generated from: '<S100>/velCtrlOut '
      rtb_ImpAsg_InsertedFor_velCtrlO[ForEach_itr] = rtb_Product1;

      // Outputs for Atomic SubSystem: '<S100>/pidWithDebug'
      // ForEachSliceAssignment generated from: '<S100>/pidDebug' incorporates:
      //   BusCreator: '<S103>/Bus Creator'

      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr].output = rtb_Product1;
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr].proportionalOutput =
        tCmd_unitRange;
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr].derivativeOutput = c_psi;

      // End of Outputs for SubSystem: '<S100>/pidWithDebug'

      // ForEachSliceAssignment generated from: '<S100>/filtMeas'
      rtb_ImpAsg_InsertedFor_filtMeas[ForEach_itr] = s_psi;

      // ForEachSliceAssignment generated from: '<S100>/filtCmd'
      rtb_ImpAsg_InsertedFor_filtCmd_[ForEach_itr] = rlim;
    }

    // End of Outputs for SubSystem: '<S95>/For Each Subsystem'

    // Trigonometry: '<S99>/Sin' incorporates:
    //   Inport: '<Root>/stateEstimate'

    c_psi = std::sin(fcsModel_U.stateEstimate.attitude_rad[2]);

    // Trigonometry: '<S99>/Sin1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    s_psi = std::cos(fcsModel_U.stateEstimate.attitude_rad[2]);

    // BusAssignment: '<S99>/Bus Assignment1' incorporates:
    //   BusAssignment: '<S93>/Bus Assignment'
    //   Gain: '<S99>/Gain'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S99>/Product1'
    //   Product: '<S99>/Product2'
    //   Sum: '<S99>/Sum'

    std::memset(&rtb_BusAssignment1, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment1.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] * c_psi - s_psi
      * rtb_ImpAsg_InsertedFor_velCtrlO[1]) * -0.10197838058331635;
    rtb_BusAssignment1.meas = fcsModel_U.stateEstimate.attitude_rad[0];
    rtb_BusAssignment1.integratorReset = resetIntegrator;

    // BusAssignment: '<S99>/Bus Assignment2' incorporates:
    //   BusAssignment: '<S93>/Bus Assignment'
    //   Gain: '<S99>/Gain1'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S99>/Product3'
    //   Product: '<S99>/Product4'
    //   Sum: '<S99>/Sum1'

    std::memset(&rtb_BusAssignment2, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment2.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] * s_psi + c_psi
      * rtb_ImpAsg_InsertedFor_velCtrlO[1]) * -0.10197838058331635;
    rtb_BusAssignment2.meas = fcsModel_U.stateEstimate.attitude_rad[1];
    rtb_BusAssignment2.integratorReset = resetIntegrator;

    // BusAssignment: '<S99>/Bus Assignment3' incorporates:
    //   BusAssignment: '<S93>/Bus Assignment'
    //   Inport: '<Root>/stateEstimate'

    std::memset(&rtb_BusAssignment3_d, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment3_d.cmd = fcsModel_DW.rcOutCmds.yawStick;
    rtb_BusAssignment3_d.meas = fcsModel_U.stateEstimate.attitude_rad[2];
    rtb_BusAssignment3_d.integratorReset = resetIntegrator;

    // Switch: '<S3>/Switch' incorporates:
    //   Constant: '<S94>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S94>/Compare'

    if (flightMode != enumFlightMode::VEL_CONTROL) {
      // Switch: '<S3>/Switch'
      fcsModel_DW.Switch = rtb_outBus;
    } else {
      // Switch: '<S3>/Switch' incorporates:
      //   Concatenate: '<S99>/Vector Concatenate'
      //   Constant: '<S99>/Constant1'
      //   Constant: '<S99>/Constant5'
      //   Product: '<S99>/Product'
      //   Sum: '<S99>/Sum2'

      fcsModel_DW.Switch.outerLoopCmds.thrustCmd_N =
        (rtb_ImpAsg_InsertedFor_velCtrlO[2] + -9.806) * 1.771;
      fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[0] = rtb_BusAssignment1;
      fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[1] = rtb_BusAssignment2;
      fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment3_d;
    }

    // End of Switch: '<S3>/Switch'
  }

  // Outputs for Iterator SubSystem: '<S9>/For Each Subsystem' incorporates:
  //   ForEach: '<S54>/For Each'

  for (ForEach_itr_h = 0; ForEach_itr_h < 3; ForEach_itr_h++) {
    // Outputs for Atomic SubSystem: '<S54>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S54>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   UnitDelay: '<S54>/Unit Delay'

    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[ForEach_itr_h].cmd,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.cmdSignalConditioningParamsArray
       [ForEach_itr_h], &rlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
       SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S54>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S54>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[ForEach_itr_h].meas,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.measSignalConditioningParamsArray
       [ForEach_itr_h], &s_psi, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
       SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S54>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S54>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, rlim, s_psi,
                          fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[ForEach_itr_h]
                          .integratorReset,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.ctrlParamsArray
                          [ForEach_itr_h],
                          fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
                          UnitDelay_DSTATE, &c_psi, &rtb_BusCreator_o, 0.004,
                          &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].pidWithDebug);

    // End of Outputs for SubSystem: '<S54>/pidWithDebug'

    // Update for UnitDelay: '<S54>/Unit Delay'
    fcsModel_DW.CoreSubsys_p[ForEach_itr_h].UnitDelay_DSTATE = c_psi;

    // ForEachSliceAssignment generated from: '<S54>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[ForEach_itr_h] =
      rtb_BusCreator_o;

    // ForEachSliceAssignment generated from: '<S54>/angRateCmds_radps'
    rtb_ImpAsg_InsertedFor_angRateC[ForEach_itr_h] = c_psi;

    // ForEachSliceAssignment generated from: '<S54>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[ForEach_itr_h] =
      s_psi;

    // ForEachSliceAssignment generated from: '<S54>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[ForEach_itr_h] =
      rlim;
  }

  // End of Outputs for SubSystem: '<S9>/For Each Subsystem'

  // Switch: '<S8>/Switch'
  rtb_Switch_f[0] = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_Switch_f[1] = rtb_ImpAsg_InsertedFor_angRateC[1];

  // Switch: '<S9>/Switch' incorporates:
  //   Constant: '<S51>/Constant'
  //   Constant: '<S53>/Constant'
  //   Logic: '<S9>/Logical Operator'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   RelationalOperator: '<S51>/Compare'
  //   RelationalOperator: '<S53>/Compare'
  //   Switch: '<S8>/Switch'

  if ((flightMode == enumFlightMode::STABILIZE) || (flightMode == enumFlightMode::
       VEL_CONTROL)) {
    rtb_Switch_f[2] = fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[2].cmd;
  } else {
    rtb_Switch_f[2] = rtb_ImpAsg_InsertedFor_angRateC[2];
  }

  // End of Switch: '<S9>/Switch'

  // MATLAB Function: '<S8>/EulerRates2BodyRates' incorporates:
  //   Inport: '<Root>/stateEstimate'
  //   Switch: '<S8>/Switch'

  // MATLAB Function 'EulerRates2BodyRates': '<S50>:1'
  // '<S50>:1:3' bodyRates_radps = eulerRates2bodyRates_function(taitBryanRates_radps,shipOrientation_rad); 
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
  c_psi = fcsModel_U.stateEstimate.attitude_rad[1];

  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  s_psi = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(s_psi - 1.5707963267948966) <= 0.004245395477824045) || (std::
       abs(s_psi - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        c_psi = -1.0;
      } else {
        c_psi = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      c_psi *= s_psi - 0.004245395477824045;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        c_psi = -1.0;
      } else {
        c_psi = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      c_psi *= s_psi - 0.004245395477824045;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        c_psi = -1.0;
      } else {
        c_psi = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      c_psi *= s_psi + 0.004245395477824045;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  s_psi = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  rlim = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  tlim = std::cos(c_psi);
  rtb_rFepTpNed[0] = 1.0;
  rtb_rFepTpNed[3] = 0.0;
  rtb_rFepTpNed[6] = -std::sin(c_psi);
  rtb_rFepTpNed[1] = 0.0;
  rtb_rFepTpNed[4] = rlim;
  rtb_rFepTpNed[7] = s_psi * tlim;
  rtb_rFepTpNed[2] = 0.0;
  rtb_rFepTpNed[5] = -s_psi;
  rtb_rFepTpNed[8] = rlim * tlim;

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
  for (int32_T ii{0}; ii < 3; ii++) {
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    c_psi = rtb_rFepTpNed[ii];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (c_psi <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      c_psi = 0.0;
    }

    rtb_rFepTpNed[ii] = c_psi;
    s_psi = c_psi * rtb_Switch_f[0];
    c_psi = rtb_rFepTpNed[ii + 3];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(c_psi) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      c_psi = 0.0;
    }

    rtb_rFepTpNed[ii + 3] = c_psi;
    s_psi += c_psi * rtb_Switch_f[1];
    c_psi = rtb_rFepTpNed[ii + 6];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(c_psi) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      c_psi = 0.0;
    }

    rtb_rFepTpNed[ii + 6] = c_psi;
    rtb_ImpAsg_InsertedFor_angRateC[ii] = c_psi * rtb_Switch_f[2] + s_psi;
  }

  // End of MATLAB Function: '<S8>/EulerRates2BodyRates'

  // BusCreator: '<S8>/Bus Creator' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_VectorConcatenate[0].integratorReset =
    fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator3' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_VectorConcatenate[1].integratorReset =
    fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator4' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[2].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_angRateC[2];
  rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_VectorConcatenate[2].integratorReset =
    fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[2].trackingCtrlCmd = 0.0;

  // Outputs for Atomic SubSystem: '<S2>/Angular Rate Controller'
  // Outputs for Iterator SubSystem: '<S7>/For Each Subsystem' incorporates:
  //   ForEach: '<S10>/For Each'

  for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
    // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S10>/ctrlInputs' incorporates:
    //   BusCreator: '<S8>/Bus Creator1'
    //   Concatenate: '<S8>/Vector Concatenate'
    //   Inport: '<Root>/ctrlParams'
    //   UnitDelay: '<S10>/Unit Delay'

    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].cmd,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray
      [ForEach_itr_p], &c_psi, 0.004, &fcsModel_DW.CoreSubsys[ForEach_itr_p].
      SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].meas,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.measSignalConditioningParamsArray
      [ForEach_itr_p], &s_psi, 0.004, &fcsModel_DW.CoreSubsys[ForEach_itr_p].
      SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, c_psi, s_psi, rtb_VectorConcatenate[ForEach_itr_p]
                          .integratorReset,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.ctrlParamsArray
                          [ForEach_itr_p], fcsModel_DW.CoreSubsys[ForEach_itr_p]
                          .UnitDelay_DSTATE, &rlim, &rtb_BusCreator_o, 0.004,
                          &fcsModel_DW.CoreSubsys[ForEach_itr_p].pidWithDebug);

    // End of Outputs for SubSystem: '<S10>/pidWithDebug'

    // Update for UnitDelay: '<S10>/Unit Delay'
    fcsModel_DW.CoreSubsys[ForEach_itr_p].UnitDelay_DSTATE = rlim;

    // ForEachSliceAssignment generated from: '<S10>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[ForEach_itr_p]
      = rtb_BusCreator_o;

    // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
    rtb_ImpAsg_InsertedFor_angAccel[ForEach_itr_p] = rlim;

    // ForEachSliceAssignment generated from: '<S10>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[ForEach_itr_p] =
      s_psi;

    // ForEachSliceAssignment generated from: '<S10>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[ForEach_itr_p] =
      c_psi;
  }

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

  s_psi = rtb_ImpAsg_InsertedFor_angRateC[0];
  tlim = rtb_ImpAsg_InsertedFor_angRateC[1];
  rlim = rtb_ImpAsg_InsertedFor_angRateC[2];

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);
  for (int32_T ii{0}; ii < 4; ii++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   BusCreator: '<S2>/Bus Creator1'
    //   Constant: '<S1>/Constant'

    c_psi = ((fcsModel_ConstP.Constant_Value[ii + 4] * s_psi +
              fcsModel_ConstP.Constant_Value[ii] *
              fcsModel_DW.Switch.outerLoopCmds.thrustCmd_N) +
             fcsModel_ConstP.Constant_Value[ii + 8] * tlim) +
      fcsModel_ConstP.Constant_Value[ii + 12] * rlim;

    // Saturate: '<S1>/Saturation'
    if (c_psi > 616850.27506808483) {
      c_psi = 616850.27506808483;
    } else if (c_psi < 0.0) {
      c_psi = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    c_psi = 9.5492965855137211 * std::sqrt(c_psi) - -0.92734095767679814 *
      fcsModel_DW.DiscreteTransferFcn_states[ii];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      // Outport: '<Root>/actuatorsCmds'
      fcsModel_Y.actuatorsCmds[ii] = -1.0;
    } else {
      // Outport: '<Root>/actuatorsCmds' incorporates:
      //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

      fcsModel_Y.actuatorsCmds[ii] = 0.036329521161600868 * c_psi +
        0.036329521161600868 * fcsModel_DW.DiscreteTransferFcn_states[ii];
    }

    // End of Switch: '<S1>/Switch'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    DiscreteTransferFcn_tmp[ii] = c_psi;
  }

  // RateTransition: '<Root>/Rate Transition'
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_Y.fcsDebug.outerLoopCtrlDebug = fcsModel_DW.RateTransition_Buffer0;

    // BusCreator: '<S3>/Bus Creator' incorporates:
    //   BusCreator: '<S95>/Bus Creator'
    //   ForEachSliceAssignment generated from: '<S100>/filtCmd'
    //   RateTransition: '<Root>/Rate Transition'

    rtb_BusCreator_b_frcCmd_N = rtb_outDebug;
    rtb_BusCreator_b_velCtrlDebug_c = rtb_ImpAsg_InsertedFor_filtCmd_[0];
    rtb_BusCreator_b_velCtrlDebug_m = rtb_ImpAsg_InsertedFor_filtMeas[0];
    rtb_BusCreator_b_velCtrlDebug_p = rtb_ImpAsg_InsertedFor_pidDebug[0];
    rtb_BusCreator_b_velCtrlDebug_0 = rtb_ImpAsg_InsertedFor_filtCmd_[1];
    rtb_BusCreator_b_velCtrlDebug_1 = rtb_ImpAsg_InsertedFor_filtMeas[1];
    rtb_BusCreator_b_velCtrlDebug_2 = rtb_ImpAsg_InsertedFor_pidDebug[1];
    rtb_BusCreator_b_velCtrlDebug_3 = rtb_ImpAsg_InsertedFor_filtCmd_[2];
    rtb_BusCreator_b_velCtrlDebug_4 = rtb_ImpAsg_InsertedFor_filtMeas[2];
    rtb_BusCreator_b_velCtrlDebug_5 = rtb_ImpAsg_InsertedFor_pidDebug[2];
  }

  // End of RateTransition: '<Root>/Rate Transition'

  // BusCreator: '<Root>/Bus Creator' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   Outport: '<Root>/fcsDebug'

  fcsModel_Y.fcsDebug.allocDebug.thrustCmd_N =
    fcsModel_DW.Switch.outerLoopCmds.thrustCmd_N;
  fcsModel_Y.fcsDebug.allocDebug.xMomCmd_Nm = rtb_ImpAsg_InsertedFor_angRateC[0];
  fcsModel_Y.fcsDebug.allocDebug.yMomCmd_Nm = rtb_ImpAsg_InsertedFor_angRateC[1];
  fcsModel_Y.fcsDebug.allocDebug.zMomCmd_Nm = rtb_ImpAsg_InsertedFor_angRateC[2];
  fcsModel_Y.fcsDebug.state = state;
  fcsModel_Y.fcsDebug.flightMode = flightMode;

  // Update for DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
  fcsModel_DW.DiscreteTransferFcn_states[0] = DiscreteTransferFcn_tmp[0];
  fcsModel_DW.DiscreteTransferFcn_states[1] = DiscreteTransferFcn_tmp[1];
  fcsModel_DW.DiscreteTransferFcn_states[2] = DiscreteTransferFcn_tmp[2];
  fcsModel_DW.DiscreteTransferFcn_states[3] = DiscreteTransferFcn_tmp[3];

  // Update for RateTransition: '<Root>/Rate Transition' incorporates:
  //   BusCreator: '<S3>/Bus Creator'
  //
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_DW.RateTransition_Buffer0.frcCmd_N = rtb_BusCreator_b_frcCmd_N;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[0] =
      rtb_BusCreator_b_velCtrlDebug_c;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[0] =
      rtb_BusCreator_b_velCtrlDebug_m;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[0] =
      rtb_BusCreator_b_velCtrlDebug_p;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[1] =
      rtb_BusCreator_b_velCtrlDebug_0;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[1] =
      rtb_BusCreator_b_velCtrlDebug_1;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[1] =
      rtb_BusCreator_b_velCtrlDebug_2;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[2] =
      rtb_BusCreator_b_velCtrlDebug_3;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[2] =
      rtb_BusCreator_b_velCtrlDebug_4;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[2] =
      rtb_BusCreator_b_velCtrlDebug_5;
  }

  // End of Update for RateTransition: '<Root>/Rate Transition'
  rate_scheduler((&fcsModel_M));
}

// Model initialize function
void fcsModel::initialize()
{
  {
    // local scratch DWork variables
    int32_T ForEach_itr;
    int32_T ForEach_itr_h;
    int32_T ForEach_itr_p;

    // SystemInitialize for Iterator SubSystem: '<S95>/For Each Subsystem'
    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // SystemInitialize for Iterator SubSystem: '<S95>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S100>/Signal Conditioning Block' 
      SignalConditioningBlock1_j_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S100>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S100>/Signal Conditioning Block1' 
      SignalConditioningBlock1_j_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S100>/Signal Conditioning Block1' 
      // End of SystemInitialize for SubSystem: '<S95>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S95>/For Each Subsystem'
    // SystemInitialize for Iterator SubSystem: '<S9>/For Each Subsystem'
    for (ForEach_itr_h = 0; ForEach_itr_h < 3; ForEach_itr_h++) {
      // SystemInitialize for Iterator SubSystem: '<S9>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S54>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S54>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S54>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S54>/Signal Conditioning Block1' 
      // End of SystemInitialize for SubSystem: '<S9>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S9>/For Each Subsystem'
    // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
    // SystemInitialize for Iterator SubSystem: '<S7>/For Each Subsystem'
    for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
      // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
      // SystemInitialize for Iterator SubSystem: '<S7>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S10>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys[ForEach_itr_p].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S10>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S10>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys[ForEach_itr_p].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S10>/Signal Conditioning Block1' 
      // End of SystemInitialize for SubSystem: '<S7>/For Each Subsystem'
      // End of SystemInitialize for SubSystem: '<S2>/Angular Rate Controller'
    }

    // End of SystemInitialize for SubSystem: '<S7>/For Each Subsystem'
    // End of SystemInitialize for SubSystem: '<S2>/Angular Rate Controller'
  }
}

// Model terminate function
void fcsModel::terminate()
{
  // (no terminate code required)
}

// Constructor
fcsModel::fcsModel() :
  fcsModel_U(),
  fcsModel_Y(),
  fcsModel_DW(),
  fcsModel_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
fcsModel::~fcsModel()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
fcsModel::RT_MODEL_fcsModel_T * fcsModel::getRTM()
{
  return (&fcsModel_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
