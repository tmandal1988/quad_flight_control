//
// File: fcsModel.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.91
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Mon Nov 27 14:03:40 2023
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
//    '<S113>/Discrete First Order Deriv Filter'
//    '<S153>/Discrete First Order Deriv Filter'
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
  real_T rtb_Product5_o;
  real_T rtb_Sum1_b;
  real_T rtb_Sum_k;
  real_T rtb_Switch2;
  real_T rtb_Switch2_p;
  real_T rtb_UkYk1;
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
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_o, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S13>/Discrete First Order Deriv Filter'

  // Product: '<S13>/Product'
  rtb_Product5_o *= rtu_pidParamBus->Kd;

  // Product: '<S13>/Product1'
  rtb_UnitDelay_i = rtb_Sum_k * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S13>/Discrete-Time Integrator'
  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = 0.0;
  }

  // Sum: '<S13>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rtb_Sum1_b = ((rtu_feedForward + rtb_Product5_o) + rtb_UnitDelay_i) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S46>/Switch2' incorporates:
  //   RelationalOperator: '<S46>/LowerRelop1'
  //   RelationalOperator: '<S46>/UpperRelop'
  //   Switch: '<S46>/Switch'

  if (rtb_Sum1_b > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2 = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_b < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S46>/Switch'
    rtb_Switch2 = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2 = rtb_Sum1_b;
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

  if (rtb_UkYk1 <= rtb_Switch2_p) {
    // Product: '<S45>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S45>/sample time'
    //
    //  About '<S45>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S48>/Switch' incorporates:
    //   RelationalOperator: '<S48>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Switch2_p) {
      rtb_Switch2_p = rtb_UkYk1;
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
  rty_pidDebug->derivativeOutput = rtb_Product5_o;

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
  localDW->UnitDelay_DSTATE = rtb_Switch2;

  // Update for UnitDelay: '<S13>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_b;
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
//    '<S106>/Compute Natural Frequency'
//    '<S129>/Compute Natural Frequency'
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
//    '<S129>/Compute Numerator And Denominator'
//    '<S114>/Compute Numerator And Denominator'
//    '<S169>/Compute Numerator And Denominator'
//    '<S154>/Compute Numerator And Denominator'
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
//    '<S106>/Compute Filter Numerator And Denominator'
//    '<S130>/Compute Filter Numerator And Denominator'
//    '<S115>/Compute Filter Numerator And Denominator'
//    '<S170>/Compute Filter Numerator And Denominator'
//    '<S155>/Compute Filter Numerator And Denominator'
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
//    '<S106>/Compute Filter Numerator And Denominator'
//    '<S130>/Compute Filter Numerator And Denominator'
//    '<S115>/Compute Filter Numerator And Denominator'
//    '<S170>/Compute Filter Numerator And Denominator'
//    '<S155>/Compute Filter Numerator And Denominator'
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
  real_T rtb_Switch2;

  // MATLAB Function: '<S29>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S29>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S30>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S30>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2,
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

  rtb_Switch2 = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S41>/Switch2' incorporates:
  //   RelationalOperator: '<S41>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_j <= rtb_Switch2) {
    // Product: '<S31>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S31>/sample time'
    //
    //  About '<S31>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S41>/Switch' incorporates:
    //   RelationalOperator: '<S41>/UpperRelop'

    if (rtb_DiscreteTransferFcn_j >= rtb_Switch2) {
      // Switch: '<S41>/Switch2'
      rtb_Switch2 = rtb_DiscreteTransferFcn_j;
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
// Output and update for atomic system:
//    '<S101>/pidWithDebug'
//    '<S150>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
  *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd, busPidDebug
  *rty_pidDebug, real_T rtp_sampleTime_s, DW_pidWithDebug_fcsModel_i_T *localDW)
{
  real_T rtb_Product5_e;
  real_T rtb_Sum1_o;
  real_T rtb_Sum_b;
  real_T rtb_Switch2;
  real_T rtb_Switch2_n;
  real_T rtb_UkYk1;
  real_T rtb_UnitDelay_a;

  // Product: '<S145>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S145>/sample time'
  //
  //  About '<S145>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_pidParamBus->outputRateLimits[1] * 0.02;

  // Sum: '<S113>/Sum'
  rtb_Sum_b = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S113>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_b,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_e, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S113>/Discrete First Order Deriv Filter'

  // Product: '<S113>/Product'
  rtb_Product5_e *= rtu_pidParamBus->Kd;

  // Product: '<S113>/Product1'
  rtb_UnitDelay_a = rtb_Sum_b * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S113>/Discrete-Time Integrator'
  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = 0.0;
  }

  // Sum: '<S113>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S113>/Discrete-Time Integrator'

  rtb_Sum1_o = ((rtu_feedForward + rtb_Product5_e) + rtb_UnitDelay_a) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S146>/Switch2' incorporates:
  //   RelationalOperator: '<S146>/LowerRelop1'
  //   RelationalOperator: '<S146>/UpperRelop'
  //   Switch: '<S146>/Switch'

  if (rtb_Sum1_o > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_n = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_o < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S146>/Switch'
    rtb_Switch2_n = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_n = rtb_Sum1_o;
  }

  // End of Switch: '<S146>/Switch2'

  // Sum: '<S145>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S145>/Delay Input2'
  //
  //  Block description for '<S145>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S145>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2_n - localDW->DelayInput2_DSTATE;

  // Switch: '<S148>/Switch2' incorporates:
  //   RelationalOperator: '<S148>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Switch2) {
    // Product: '<S145>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S145>/sample time'
    //
    //  About '<S145>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_pidParamBus->outputRateLimits[0] * 0.02;

    // Switch: '<S148>/Switch' incorporates:
    //   RelationalOperator: '<S148>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Switch2) {
      rtb_Switch2 = rtb_UkYk1;
    }

    // End of Switch: '<S148>/Switch'
  }

  // End of Switch: '<S148>/Switch2'

  // Sum: '<S145>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S145>/Delay Input2'
  //
  //  Block description for '<S145>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S145>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S113>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S113>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_a;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_e;

  // Update for DiscreteIntegrator: '<S113>/Discrete-Time Integrator' incorporates:
  //   Product: '<S113>/Product2'
  //   Product: '<S113>/Product3'
  //   Product: '<S113>/Product5'
  //   Sum: '<S113>/Sum2'
  //   Sum: '<S113>/Sum3'
  //   Sum: '<S113>/Sum4'
  //   Sum: '<S113>/Sum5'
  //   UnitDelay: '<S113>/Unit Delay'
  //   UnitDelay: '<S113>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_b * rtu_pidParamBus->Ki) * 0.02;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S145>/Delay Input2'
  //
  //  Block description for '<S145>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S113>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_n;

  // Update for UnitDelay: '<S113>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_o;
}

//
// System initialize for atomic system:
//    '<S101>/Signal Conditioning Block1'
//    '<S101>/Signal Conditioning Block'
//    '<S150>/Signal Conditioning Block1'
//    '<S150>/Signal Conditioning Block'
//
void fcsModel::SignalConditioningBlock1_c_Init(DW_SignalConditioningBlock1_g_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S130>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S101>/Signal Conditioning Block1'
//    '<S101>/Signal Conditioning Block'
//    '<S150>/Signal Conditioning Block1'
//    '<S150>/Signal Conditioning Block'
//
void fcsModel::fcsM_SignalConditioningBlock1_f(real_T rtu_input, const
  busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
  rtp_sampleTime_s, DW_SignalConditioningBlock1_g_T *localDW)
{
  std::array<real_T, 3> rtb_accelNum;
  std::array<real_T, 3> rtb_den;
  std::array<real_T, 3> rtb_rateNum;
  real_T rtb_DiscreteTransferFcn_d;
  real_T rtb_Switch2;

  // MATLAB Function: '<S129>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S129>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S130>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S130>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S130>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_d = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S134>/Switch2' incorporates:
  //   RelationalOperator: '<S134>/LowerRelop1'
  //   RelationalOperator: '<S134>/UpperRelop'
  //   Switch: '<S134>/Switch'

  if (rtb_DiscreteTransferFcn_d > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_d < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S134>/Switch'
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S134>/Switch2'

  // Sum: '<S131>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S131>/Delay Input2'
  //
  //  Block description for '<S131>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S131>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_d -= localDW->DelayInput2_DSTATE;

  // Switch: '<S141>/Switch2' incorporates:
  //   Product: '<S131>/delta rise limit'
  //   SampleTimeMath: '<S131>/sample time'
  //
  //  About '<S131>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_params->filteredInputRateLimits[1] * 0.02;

  // Switch: '<S141>/Switch2' incorporates:
  //   RelationalOperator: '<S141>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_d <= rtb_Switch2) {
    // Product: '<S131>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S131>/sample time'
    //
    //  About '<S131>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_params->filteredInputRateLimits[0] * 0.02;

    // Switch: '<S141>/Switch' incorporates:
    //   RelationalOperator: '<S141>/UpperRelop'

    if (rtb_DiscreteTransferFcn_d >= rtb_Switch2) {
      // Switch: '<S141>/Switch2'
      rtb_Switch2 = rtb_DiscreteTransferFcn_d;
    }

    // End of Switch: '<S141>/Switch'
  }

  // End of Switch: '<S141>/Switch2'

  // Sum: '<S131>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S131>/Delay Input2'
  //
  //  Block description for '<S131>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S131>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S130>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S131>/Delay Input2'
  //
  //  Block description for '<S131>/Delay Input2':
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

  // MATLAB Function 'checkRcCmds': '<S189>:7'
  // '<S189>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S189>:7:3' if(rcCmds.throttleCmd_nd <= pwmLowVal && ...
  // '<S189>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S189>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S189>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 1000) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 1000) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 1000) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 1000) {
          // '<S189>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S189>:7:8' else
          // '<S189>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S189>:7:8' else
        // '<S189>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S189>:7:8' else
      // '<S189>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S189>:7:8' else
    // '<S189>:7:9' isTrue = false;
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
  int32_T ForEach_itr_i;
  std::array<real_T, 4> DiscreteTransferFcn_tmp;
  std::array<real_T, 9> conversionMatrix;
  std::array<real_T, 3> den;
  std::array<real_T, 3> num;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angRateC;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_cmd_at_i;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtCmd_;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtMeas;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_meas_at_;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_neVelCmd;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDeb_m;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDebug;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_velCtrlO;
  std::array<busCtrlInputs, 3> rtb_VectorConcatenate;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_2;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_5;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_p;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_2;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_5;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_p;
  busPidDebug rtb_BusCreator_o;
  real_T DiscreteTransferFcn;
  real_T plim;
  real_T rlim;
  real_T rtb_BusCreator_b_frcCmd_N;
  real_T rtb_BusCreator_b_posCtrlDebug_0;
  real_T rtb_BusCreator_b_posCtrlDebug_1;
  real_T rtb_BusCreator_b_posCtrlDebug_3;
  real_T rtb_BusCreator_b_posCtrlDebug_4;
  real_T rtb_BusCreator_b_posCtrlDebug_c;
  real_T rtb_BusCreator_b_posCtrlDebug_m;
  real_T rtb_BusCreator_b_velCtrlDebug_0;
  real_T rtb_BusCreator_b_velCtrlDebug_1;
  real_T rtb_BusCreator_b_velCtrlDebug_3;
  real_T rtb_BusCreator_b_velCtrlDebug_4;
  real_T rtb_BusCreator_b_velCtrlDebug_c;
  real_T rtb_BusCreator_b_velCtrlDebug_m;
  real_T rtb_frcCmd_N;
  real_T vzlim;
  real_T ylim;
  int32_T tlim;
  boolean_T resetIntegrator;
  boolean_T rtb_atCenter;
  enumFlightMode flightMode;
  enumStateMachine state;

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
    // Transition: '<S189>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S189>:1'
    // '<S189>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S189>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
    }

    // '<S189>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S189>:3'
      // '<S189>:10:1' sf_internal_predicateOutput = after(60, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 15000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S189>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S189>:1'
        // '<S189>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S189>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S189>:1:4' resetIntegrator = true;

        // '<S189>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimits(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 1000) {
        // Transition: '<S189>:12'
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S189>:11'
        // '<S189>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S189>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S189>:11:4' resetIntegrator = false;
      } else {
        // '<S189>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S189>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S189>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S189>:1'
      // '<S189>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 5 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 1250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S189>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S189>:3'
        // '<S189>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S189>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S189>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S189>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S189>:1:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S189>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S189>:11'
      // '<S189>:20:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd <= paramsStruct.pwmLimits(1); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1000) {
        // Transition: '<S189>:20'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S189>:3'
        // '<S189>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S189>:3:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S189>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S189>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S189>:11:3' rcCheckFlag = checkRcCmds(rcCmds, paramsStruct);
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S189>:11:4' resetIntegrator = false;
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
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S190>:1'
  // '<S190>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcCmds, expo, rcParamsStruct); 
  // INTERPRETRCINPUTS_FUNCTION
  // Computes command and flight mode from the rc inputs
  //  Used to directly set force commands in N that is fed into allocation in
  //  STABILIZE and ACRO flight modes
  // 'interpretRcInputs_function:7' rcOutCmds.throttleStick = 0;
  //  Used to command roll and pitch angles in STABILIZE flight mode
  // 'interpretRcInputs_function:10' rcOutCmds.rollStick = 0;
  // 'interpretRcInputs_function:11' rcOutCmds.pitchStick = 0;
  //  Used to control yaw rate in all flight modes
  // 'interpretRcInputs_function:14' rcOutCmds.yawStick = 0;
  //  Not used currently
  // 'interpretRcInputs_function:17' rcOutCmds.vxStick_mps = 0;
  // 'interpretRcInputs_function:18' rcOutCmds.vyStick_mps = 0;
  //  Used to command vertical velocity in ALT_CONTROL flight mode. In this
  //  mode rcOutCmds.throttleStick is ignored.
  // 'interpretRcInputs_function:22' rcOutCmds.vzStick_mps = 0;
  // Select mode and set the command limits
  // 'interpretRcInputs_function:25' if(rcInCmds.rcSwitch1_nd < 1100)
  if (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1100) {
    // 'interpretRcInputs_function:26' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:27' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -60;

    // 'interpretRcInputs_function:28' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    // 'interpretRcInputs_function:29' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    // 'interpretRcInputs_function:30' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:31' vxlim = 0;
    // 'interpretRcInputs_function:32' vylim = 0;
    // 'interpretRcInputs_function:33' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    vzlim = -2.0;
  } else if ((fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1100) &&
             (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1800)) {
    // 'interpretRcInputs_function:35' elseif (rcInCmds.rcSwitch1_nd >= 1100 && rcInCmds.rcSwitch1_nd < 1800) 
    // 'interpretRcInputs_function:36' flightMode = enumFlightMode.ALT_CONTROL;
    flightMode = enumFlightMode::ALT_CONTROL;

    // 'interpretRcInputs_function:37' vxlim = 0;
    // 'interpretRcInputs_function:38' vylim = 0;
    // 'interpretRcInputs_function:39' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    vzlim = -2.0;

    // 'interpretRcInputs_function:40' tlim = 0;
    tlim = 0;

    // 'interpretRcInputs_function:41' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    // 'interpretRcInputs_function:42' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    // 'interpretRcInputs_function:43' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else {
    // 'interpretRcInputs_function:44' else
    // 'interpretRcInputs_function:45' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:46' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -60;

    // 'interpretRcInputs_function:47' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    // 'interpretRcInputs_function:48' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    // 'interpretRcInputs_function:49' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:50' vxlim = 0;
    // 'interpretRcInputs_function:51' vylim = 0;
    // 'interpretRcInputs_function:52' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    vzlim = -2.0;
  }

  // 'interpretRcInputs_function:55' tCmd = min( rcParamsStruct.pwmLimitsThrottle(2), ... 
  // 'interpretRcInputs_function:56'         max( rcParamsStruct.pwmLimitsThrottle(1), double(rcInCmds.throttleCmd_nd) ) ); 
  plim = std::fmin(1882.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.throttleCmd_nd)));

  // 'interpretRcInputs_function:58' if (flightMode == enumFlightMode.STABILIZE || flightMode == enumFlightMode.ACRO) 
  if (flightMode == enumFlightMode::STABILIZE) {
    //  In stabilize mode throttle stick starts at 0
    // 'interpretRcInputs_function:60' tCmd_unitRange = -1 + tCmd/1000;
    ylim = plim / 1000.0 + -1.0;
  } else {
    // 'interpretRcInputs_function:61' else
    //  In other modes throttle is used with symmetry around middle stick
    //  position. NOT USED CURRENTLY. Defined so that tCmd_unitRange is
    //  defined in all execution path
    // 'interpretRcInputs_function:65' tCmd_unitRange = -1 + tCmd/500;
    ylim = plim / 500.0 + -1.0;
  }

  //  Use throttle stick to set Vz commands. Always set but only used
  //  in ALT_CONTROL flight mode. Has different slopes about the center point
  //  as the center point is not always at 1500 which is PWM center
  // 'interpretRcInputs_function:72' if ((tCmd <= rcParamsStruct.pwmThrottleMidHigh) && ... 
  // 'interpretRcInputs_function:73'         tCmd >= rcParamsStruct.pwmThrottleMidLow) 
  if ((plim <= 1470.0) && (plim >= 1260.0)) {
    // 'interpretRcInputs_function:74' vzCmd_unitRange  = 0;
    plim = 0.0;
  } else if (plim < 1260.0) {
    // 'interpretRcInputs_function:75' elseif (tCmd < rcParamsStruct.pwmThrottleMidLow) 
    // 'interpretRcInputs_function:76' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
    // 'interpretRcInputs_function:77'         rcParamsStruct.pwmToCmdThrottleIncptLow; 
    plim = 0.0038461538461538464 * plim + -4.8461538461538458;

    // ;
  } else {
    // 'interpretRcInputs_function:78' else
    // 'interpretRcInputs_function:79' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeHigh*tCmd + ...; 
    // 'interpretRcInputs_function:80'         rcParamsStruct.pwmToCmdThrottleIncptHigh; 
    plim = 0.0024271844660194173 * plim + -3.5679611650485437;

    // ;
  }

  //  Set roll, pitch and yaw stick
  // 'interpretRcInputs_function:84' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:85'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  // 'interpretRcInputs_function:86' rCmd_unitRange = -3 + rCmd/500;
  // 'interpretRcInputs_function:88' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:89'         max( rcParamsStruct.pwmLimits(1),  double(rcInCmds.joystickYCmd_nd) ) ); 
  //  Reverse the pitch cmd
  // 'interpretRcInputs_function:92' pCmd_unitRange = -(-3 + pCmd/500);
  // 'interpretRcInputs_function:94' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:95'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  // 'interpretRcInputs_function:96' yCmd_unitRange = -3 + yCmd/500;
  // 'interpretRcInputs_function:99' if expo
  // 'interpretRcInputs_function:122' else
  //  Usually expo is set in the Tx hence simply use a linear map here
  // 'interpretRcInputs_function:124' rcOutCmds.throttleStick = tCmd_unitRange*tlim; 
  fcsModel_DW.rcOutCmds.throttleStick = ylim * static_cast<real_T>(tlim);

  // 'interpretRcInputs_function:125' rcOutCmds.rollStick = rCmd_unitRange*rlim; 
  fcsModel_DW.rcOutCmds.rollStick = (std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickXCmd_nd))) / 500.0 + -3.0) *
    0.78539816339744828;

  // 'interpretRcInputs_function:126' rcOutCmds.pitchStick = pCmd_unitRange*plim; 
  fcsModel_DW.rcOutCmds.pitchStick = -(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickYCmd_nd))) / 500.0 + -3.0) *
    0.78539816339744828;

  // 'interpretRcInputs_function:127' rcOutCmds.yawStick = yCmd_unitRange*ylim;
  fcsModel_DW.rcOutCmds.yawStick = (std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickZCmd_nd))) / 500.0 + -3.0) *
    1.0471975511965976;

  // 'interpretRcInputs_function:128' rcOutCmds.vzStick_mps = vzCmd_unitRange*vzlim; 
  fcsModel_DW.rcOutCmds.vzStick_mps = plim * -2.0;

  // 'interpretRcInputs_function:129' rcOutCmds.vxStick_mps = pCmd_unitRange*vxlim; 
  fcsModel_DW.rcOutCmds.vxStick_mps = 0.0;

  // 'interpretRcInputs_function:130' rcOutCmds.vyStick_mps = rCmd_unitRange*vylim; 
  fcsModel_DW.rcOutCmds.vyStick_mps = 0.0;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    busOuterLoopToInnerLoop rtb_outBus;

    // Concatenate: '<S100>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S100>/Bus Assignment1' incorporates:
    //   Concatenate: '<S100>/Vector Concatenate'
    //   Constant: '<S100>/Constant1'

    rtb_VectorConcatenate[0].cmd = 0.0;
    rtb_VectorConcatenate[0].meas = 0.0;
    rtb_VectorConcatenate[0].integratorReset = resetIntegrator;

    // Concatenate: '<S100>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S100>/Bus Assignment2' incorporates:
    //   Concatenate: '<S100>/Vector Concatenate'
    //   Constant: '<S100>/Constant3'

    rtb_VectorConcatenate[1].cmd = 0.0;
    rtb_VectorConcatenate[1].meas = 0.0;
    rtb_VectorConcatenate[1].integratorReset = resetIntegrator;

    // Outputs for Atomic SubSystem: '<S100>/baroAltEstimate'
    // MATLAB Function: '<S106>/Compute Natural Frequency' incorporates:
    //   BusCreator: '<S104>/Bus Creator'
    //   Constant: '<S104>/Constant'

    fcsMode_ComputeNaturalFrequency(157.07963267948966, 1.0, &rtb_frcCmd_N);

    // MATLAB Function: '<S106>/Compute Filter Numerator And Denominator' incorporates:
    //   BusCreator: '<S104>/Bus Creator'
    //   Constant: '<S104>/Constant'

    ComputeFilterNumeratorAndDenomi(rtb_frcCmd_N, 1.0, &num[0], &den[0], 0.02);

    // DiscreteTransferFcn: '<S106>/Discrete Transfer Fcn' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = (fcsModel_U.stateEstimate.pressure_mbar -
                    fcsModel_DW.DiscreteTransferFcn_states_e[0] * den[1]) -
      fcsModel_DW.DiscreteTransferFcn_states_e[1] * den[2];

    // MATLAB Function: '<S104>/computeAltFromBaro' incorporates:
    //   DiscreteTransferFcn: '<S106>/Discrete Transfer Fcn'

    //  Call the scripted function
    // MATLAB Function 'baroAltEstimate/computeAltFromBaro': '<S107>:1'
    // '<S107>:1:4' [agl_baro_m, base_alt_point_m, alt_m] = computeAltFromBaro_function(pressure_mbar, state, params); 
    // 'computeAltFromBaro_function:3' if isempty(pers_base_alt_point_m)
    // Convert baro measurement in millibar to Pascal
    // 'computeAltFromBaro_function:8' p_pa = pressure_mbar*100.0;
    // Use barometric formula to compute altitue
    // 'computeAltFromBaro_function:11' alt_m = 44330.0 * (1 - ( p_pa/params.mean_sea_level_pressure_pa )^(1.0/5.255) ); 
    DiscreteTransferFcn = (1.0 - std::pow(((num[0] * rtb_frcCmd_N +
      fcsModel_DW.DiscreteTransferFcn_states_e[0] * num[1]) +
      fcsModel_DW.DiscreteTransferFcn_states_e[1] * num[2]) * 100.0 / 101325.0,
      0.19029495718363465)) * 44330.0;

    // Use alt at MTR_ARMED state as the base_alt_point_m
    // 'computeAltFromBaro_function:16' if (state ==  enumStateMachine.MTR_ARMED) 
    if (state == enumStateMachine::MTR_ARMED) {
      // 'computeAltFromBaro_function:17' pers_base_alt_point_m = alt_m;
      fcsModel_DW.pers_base_alt_point_m = DiscreteTransferFcn;
    }

    // 'computeAltFromBaro_function:20' base_alt_point_m = pers_base_alt_point_m; 
    // 'computeAltFromBaro_function:22' agl_baro_m = alt_m - base_alt_point_m;
    plim = DiscreteTransferFcn - fcsModel_DW.pers_base_alt_point_m;

    // End of MATLAB Function: '<S104>/computeAltFromBaro'

    // Update for DiscreteTransferFcn: '<S106>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states_e[1] =
      fcsModel_DW.DiscreteTransferFcn_states_e[0];
    fcsModel_DW.DiscreteTransferFcn_states_e[0] = rtb_frcCmd_N;

    // End of Outputs for SubSystem: '<S100>/baroAltEstimate'

    // Outputs for Atomic SubSystem: '<S100>/holdOutputAtCenter'
    // MATLAB Function: '<S105>/holdOutputAtCenter'
    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S110>:1'
    // '<S110>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
    // HOLDOUTPUTATCENTER_FUNCTION holds the output constant at last input if the 
    // trigger value is within user defined delta from the center
    // 'holdOutputAtCenter_function:5' if isempty(last_input)
    // 'holdOutputAtCenter_function:9' if(trigger <= (params.center + params.posDeltaFromCenter) && ... 
    // 'holdOutputAtCenter_function:10'         trigger >=(params.center - params.negDeltaFromCenter)) 
    if ((fcsModel_DW.rcOutCmds.vzStick_mps <= 0.05) &&
        (fcsModel_DW.rcOutCmds.vzStick_mps >= -0.05)) {
      // 'holdOutputAtCenter_function:11' atCenter = true;
      rtb_atCenter = true;
    } else {
      // 'holdOutputAtCenter_function:12' else
      // 'holdOutputAtCenter_function:13' atCenter = false;
      rtb_atCenter = false;

      // 'holdOutputAtCenter_function:14' last_input = input;
      fcsModel_DW.last_input = plim;
    }

    // Gain: '<S100>/Gain' incorporates:
    //   MATLAB Function: '<S105>/holdOutputAtCenter'

    // 'holdOutputAtCenter_function:17' output = last_input;
    DiscreteTransferFcn = -fcsModel_DW.last_input;

    // End of Outputs for SubSystem: '<S100>/holdOutputAtCenter'

    // Gain: '<S100>/Gain1'
    plim = -plim;

    // Concatenate: '<S100>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // Switch: '<S100>/Switch' incorporates:
    //   RelationalOperator: '<S102>/Compare'

    if (rtb_atCenter) {
      // BusAssignment: '<S100>/Bus Assignment3' incorporates:
      //   Concatenate: '<S100>/Vector Concatenate'
      //   Constant: '<S100>/Constant5'

      rtb_VectorConcatenate[2].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S100>/Bus Assignment3' incorporates:
      //   Concatenate: '<S100>/Vector Concatenate'

      rtb_VectorConcatenate[2].feedForwardCmd =
        fcsModel_DW.rcOutCmds.vzStick_mps;
    }

    // End of Switch: '<S100>/Switch'

    // Outputs for Atomic SubSystem: '<S100>/holdOutputAtCenter'
    // BusAssignment: '<S100>/Bus Assignment3' incorporates:
    //   Concatenate: '<S100>/Vector Concatenate'
    //   Constant: '<S103>/Constant'
    //   Gain: '<S100>/Gain'
    //   Logic: '<S100>/Logical Operator'
    //   MATLAB Function: '<S105>/holdOutputAtCenter'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S103>/Compare'

    rtb_VectorConcatenate[2].cmd = -fcsModel_DW.last_input;

    // End of Outputs for SubSystem: '<S100>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].meas = plim;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::ALT_CONTROL));

    // Outputs for Iterator SubSystem: '<S97>/NED Position Control' incorporates:
    //   ForEach: '<S101>/For Each'

    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // Outputs for Atomic SubSystem: '<S101>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S101>/ctrlInputs' incorporates:
      //   BusAssignment: '<S100>/Bus Assignment'
      //   Concatenate: '<S149>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S101>/Unit Delay'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_frcCmd_N, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S101>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S101>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr_i], &plim, 0.02, &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S101>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S101>/pidWithDebug'
      fcsModel_pidWithDebug_j(rtb_VectorConcatenate[ForEach_itr_i].
        feedForwardCmd, rtb_frcCmd_N, plim, rtb_VectorConcatenate[ForEach_itr_i]
        .integratorReset,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.ctrlParamsArray[
        ForEach_itr_i], fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE,
        &DiscreteTransferFcn, &rtb_BusCreator_o, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].pidWithDebug);

      // End of Outputs for SubSystem: '<S101>/pidWithDebug'

      // Update for UnitDelay: '<S101>/Unit Delay'
      fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE =
        DiscreteTransferFcn;

      // ForEachSliceAssignment generated from: '<S101>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDeb_m[ForEach_itr_i] = rtb_BusCreator_o;

      // ForEachSliceAssignment generated from: '<S101>/neVelCmd_mps'
      rtb_ImpAsg_InsertedFor_neVelCmd[ForEach_itr_i] = DiscreteTransferFcn;

      // ForEachSliceAssignment generated from: '<S101>/meas'
      rtb_ImpAsg_InsertedFor_meas_at_[ForEach_itr_i] = plim;

      // ForEachSliceAssignment generated from: '<S101>/cmd'
      rtb_ImpAsg_InsertedFor_cmd_at_i[ForEach_itr_i] = rtb_frcCmd_N;
    }

    // End of Outputs for SubSystem: '<S97>/NED Position Control'

    // Logic: '<S96>/Logical Operator' incorporates:
    //   Constant: '<S99>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S99>/Compare'

    rtb_atCenter = (resetIntegrator || (flightMode != enumFlightMode::
      ALT_CONTROL));

    // Concatenate: '<S96>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S96>/Bus Assignment' incorporates:
    //   Concatenate: '<S96>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[0];
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedVel_mps[0];
    rtb_VectorConcatenate[0].integratorReset = rtb_atCenter;

    // Concatenate: '<S96>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S96>/Bus Assignment1' incorporates:
    //   Concatenate: '<S96>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[1];
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedVel_mps[1];
    rtb_VectorConcatenate[1].integratorReset = rtb_atCenter;

    // Concatenate: '<S96>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S96>/Bus Assignment2' incorporates:
    //   Concatenate: '<S96>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.nedVel_mps[2];
    rtb_VectorConcatenate[2].integratorReset = rtb_atCenter;

    // Outputs for Iterator SubSystem: '<S98>/For Each Subsystem' incorporates:
    //   ForEach: '<S150>/For Each'

    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // Outputs for Atomic SubSystem: '<S150>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S150>/ctrlInputs' incorporates:
      //   BusAssignment: '<S96>/Bus Assignment3'
      //   Concatenate: '<S149>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S150>/Unit Delay'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr], &rtb_frcCmd_N, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S150>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S150>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr], &plim, 0.02, &fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S150>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S150>/pidWithDebug'
      fcsModel_pidWithDebug_j(0.0, rtb_frcCmd_N, plim,
        rtb_VectorConcatenate[ForEach_itr].integratorReset,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ctrlParamsArray[
        ForEach_itr], fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE,
        &DiscreteTransferFcn, &rtb_BusCreator_o, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].pidWithDebug);

      // End of Outputs for SubSystem: '<S150>/pidWithDebug'

      // Update for UnitDelay: '<S150>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE =
        DiscreteTransferFcn;

      // ForEachSliceAssignment generated from: '<S150>/velCtrlOut '
      rtb_ImpAsg_InsertedFor_velCtrlO[ForEach_itr] = DiscreteTransferFcn;

      // ForEachSliceAssignment generated from: '<S150>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr] = rtb_BusCreator_o;

      // ForEachSliceAssignment generated from: '<S150>/filtMeas'
      rtb_ImpAsg_InsertedFor_filtMeas[ForEach_itr] = plim;

      // ForEachSliceAssignment generated from: '<S150>/filtCmd'
      rtb_ImpAsg_InsertedFor_filtCmd_[ForEach_itr] = rtb_frcCmd_N;
    }

    // End of Outputs for SubSystem: '<S98>/For Each Subsystem'

    // Product: '<S149>/Product' incorporates:
    //   Constant: '<S149>/Constant1'
    //   Constant: '<S149>/Constant5'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S149>/Divide'
    //   Product: '<S149>/Product5'
    //   Sum: '<S149>/Sum2'
    //   Trigonometry: '<S149>/Sin2'
    //   Trigonometry: '<S149>/Sin3'

    rtb_frcCmd_N = 1.0 / (std::cos(fcsModel_U.stateEstimate.attitude_rad[0]) *
                          std::cos(fcsModel_U.stateEstimate.attitude_rad[1])) *
      (rtb_ImpAsg_InsertedFor_velCtrlO[2] + -9.806) * 1.771;

    // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
    //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
    //   Constant: '<S3>/Constant'
    //   Inport: '<Root>/stateEstimate'

    rtb_outBus = fcsModel_rtZbusOuterLoopToInnerLoop;

    // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S95>:1' 
    // '<S95>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
    rtb_outBus.outerLoopCmds.thrustCmd_N = fcsModel_DW.rcOutCmds.throttleStick;

    // '<S95>:1:3' outDebug = throttleCmd_N;
    DiscreteTransferFcn = fcsModel_DW.rcOutCmds.throttleStick;

    //  This is a stop gap setup where we are only assuming that rate control
    //  is active and therefore not setting up attCtrlInputs for Euler angle
    //  control
    // '<S95>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].cmd =
      fcsModel_DW.rcOutCmds.rollStick;

    // '<S95>:1:8' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].meas =
      fcsModel_U.stateEstimate.attitude_rad[0];

    // '<S95>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].cmd =
      fcsModel_DW.rcOutCmds.pitchStick;

    // '<S95>:1:10' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].meas =
      fcsModel_U.stateEstimate.attitude_rad[1];

    // '<S95>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].cmd =
      fcsModel_DW.rcOutCmds.yawStick;

    // '<S95>:1:12' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].meas =
      fcsModel_U.stateEstimate.attitude_rad[2];

    // '<S95>:1:14' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[0].integratorReset =
      resetIntegrator;

    // '<S95>:1:15' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[1].integratorReset =
      resetIntegrator;

    // '<S95>:1:16' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = resetIntegrator; 
    rtb_outBus.attCtrlInputs.ctrlInputsArray[2].integratorReset =
      resetIntegrator;

    // RelationalOperator: '<S93>/Compare' incorporates:
    //   Constant: '<S93>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_atCenter = (flightMode != enumFlightMode::ALT_CONTROL);

    // Switch: '<S3>/Switch'
    if (rtb_atCenter) {
      // Switch: '<S3>/Switch'
      fcsModel_DW.Switch = rtb_outBus;
    } else {
      // Switch: '<S3>/Switch' incorporates:
      //   BusAssignment: '<S149>/Bus Assignment'
      //   BusAssignment: '<S3>/Bus Assignment'

      fcsModel_DW.Switch.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      fcsModel_DW.Switch.attCtrlInputs = rtb_outBus.attCtrlInputs;
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
       [ForEach_itr_h], &vzlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
       SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S54>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S54>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[ForEach_itr_h].meas,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.measSignalConditioningParamsArray
       [ForEach_itr_h], &ylim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
       SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S54>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S54>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, vzlim, ylim,
                          fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[ForEach_itr_h]
                          .integratorReset,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.ctrlParamsArray
                          [ForEach_itr_h],
                          fcsModel_DW.CoreSubsys_p[ForEach_itr_h].
                          UnitDelay_DSTATE, &plim, &rtb_BusCreator_o, 0.004,
                          &fcsModel_DW.CoreSubsys_p[ForEach_itr_h].pidWithDebug);

    // End of Outputs for SubSystem: '<S54>/pidWithDebug'

    // Update for UnitDelay: '<S54>/Unit Delay'
    fcsModel_DW.CoreSubsys_p[ForEach_itr_h].UnitDelay_DSTATE = plim;

    // ForEachSliceAssignment generated from: '<S54>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[ForEach_itr_h] =
      rtb_BusCreator_o;

    // ForEachSliceAssignment generated from: '<S54>/angRateCmds_radps'
    rtb_ImpAsg_InsertedFor_angRateC[ForEach_itr_h] = plim;

    // ForEachSliceAssignment generated from: '<S54>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[ForEach_itr_h] =
      ylim;

    // ForEachSliceAssignment generated from: '<S54>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[ForEach_itr_h] =
      vzlim;
  }

  // End of Outputs for SubSystem: '<S9>/For Each Subsystem'

  // Switch: '<S8>/Switch'
  num[0] = rtb_ImpAsg_InsertedFor_angRateC[0];
  num[1] = rtb_ImpAsg_InsertedFor_angRateC[1];

  // Switch: '<S9>/Switch' incorporates:
  //   Constant: '<S51>/Constant'
  //   Constant: '<S53>/Constant'
  //   Logic: '<S9>/Logical Operator'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   RelationalOperator: '<S51>/Compare'
  //   RelationalOperator: '<S53>/Compare'
  //   Switch: '<S8>/Switch'

  if ((flightMode == enumFlightMode::STABILIZE) || (flightMode == enumFlightMode::
       ALT_CONTROL)) {
    num[2] = fcsModel_DW.Switch.attCtrlInputs.ctrlInputsArray[2].cmd;
  } else {
    num[2] = rtb_ImpAsg_InsertedFor_angRateC[2];
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
  plim = fcsModel_U.stateEstimate.attitude_rad[1];

  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  ylim = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(ylim - 1.5707963267948966) <= 0.004245395477824045) || (std::abs
       (ylim - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      plim *= ylim - 0.004245395477824045;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      plim *= ylim - 0.004245395477824045;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      plim *= ylim + 0.004245395477824045;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  ylim = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  vzlim = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  rlim = std::cos(plim);
  conversionMatrix[0] = 1.0;
  conversionMatrix[3] = 0.0;
  conversionMatrix[6] = -std::sin(plim);
  conversionMatrix[1] = 0.0;
  conversionMatrix[4] = vzlim;
  conversionMatrix[7] = ylim * rlim;
  conversionMatrix[2] = 0.0;
  conversionMatrix[5] = -ylim;
  conversionMatrix[8] = vzlim * rlim;

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
  for (tlim = 0; tlim < 3; tlim++) {
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    plim = conversionMatrix[tlim];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (plim <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      plim = 0.0;
    }

    conversionMatrix[tlim] = plim;
    ylim = plim * num[0];
    plim = conversionMatrix[tlim + 3];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(plim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      plim = 0.0;
    }

    conversionMatrix[tlim + 3] = plim;
    ylim += plim * num[1];
    plim = conversionMatrix[tlim + 6];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(plim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      plim = 0.0;
    }

    conversionMatrix[tlim + 6] = plim;
    rtb_ImpAsg_InsertedFor_angRateC[tlim] = plim * num[2] + ylim;
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
      [ForEach_itr_p], &plim, 0.004, &fcsModel_DW.CoreSubsys[ForEach_itr_p].
      SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S10>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].meas,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.measSignalConditioningParamsArray
      [ForEach_itr_p], &ylim, 0.004, &fcsModel_DW.CoreSubsys[ForEach_itr_p].
      SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, plim, ylim, rtb_VectorConcatenate[ForEach_itr_p].
                          integratorReset,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.ctrlParamsArray
                          [ForEach_itr_p], fcsModel_DW.CoreSubsys[ForEach_itr_p]
                          .UnitDelay_DSTATE, &vzlim, &rtb_BusCreator_o, 0.004,
                          &fcsModel_DW.CoreSubsys[ForEach_itr_p].pidWithDebug);

    // End of Outputs for SubSystem: '<S10>/pidWithDebug'

    // Update for UnitDelay: '<S10>/Unit Delay'
    fcsModel_DW.CoreSubsys[ForEach_itr_p].UnitDelay_DSTATE = vzlim;

    // ForEachSliceAssignment generated from: '<S10>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[ForEach_itr_p]
      = rtb_BusCreator_o;

    // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
    rtb_ImpAsg_InsertedFor_angAccel[ForEach_itr_p] = vzlim;

    // ForEachSliceAssignment generated from: '<S10>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[ForEach_itr_p] =
      ylim;

    // ForEachSliceAssignment generated from: '<S10>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[ForEach_itr_p] =
      plim;
  }

  // End of Outputs for SubSystem: '<S7>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'

  for (tlim = 0; tlim < 3; tlim++) {
    rtb_ImpAsg_InsertedFor_angRateC[tlim] = 0.0;
    rtb_ImpAsg_InsertedFor_angRateC[tlim] +=
      fcsModel_ConstP.Constant_Value_n[tlim] * rtb_ImpAsg_InsertedFor_angAccel[0];
    rtb_ImpAsg_InsertedFor_angRateC[tlim] +=
      fcsModel_ConstP.Constant_Value_n[tlim + 3] *
      rtb_ImpAsg_InsertedFor_angAccel[1];
    rtb_ImpAsg_InsertedFor_angRateC[tlim] +=
      fcsModel_ConstP.Constant_Value_n[tlim + 6] *
      rtb_ImpAsg_InsertedFor_angAccel[2];
  }

  // End of Product: '<S2>/Matrix Multiply'

  // SignalConversion generated from: '<S1>/Matrix Multiply' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'

  ylim = rtb_ImpAsg_InsertedFor_angRateC[0];
  vzlim = rtb_ImpAsg_InsertedFor_angRateC[1];
  rlim = rtb_ImpAsg_InsertedFor_angRateC[2];

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);
  for (tlim = 0; tlim < 4; tlim++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   BusCreator: '<S2>/Bus Creator1'
    //   Constant: '<S1>/Constant'

    plim = ((fcsModel_ConstP.Constant_Value_c[tlim + 4] * ylim +
             fcsModel_ConstP.Constant_Value_c[tlim] *
             fcsModel_DW.Switch.outerLoopCmds.thrustCmd_N) +
            fcsModel_ConstP.Constant_Value_c[tlim + 8] * vzlim) +
      fcsModel_ConstP.Constant_Value_c[tlim + 12] * rlim;

    // Saturate: '<S1>/Saturation'
    if (plim > 616850.27506808483) {
      plim = 616850.27506808483;
    } else if (plim < 0.0) {
      plim = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    plim = 9.5492965855137211 * std::sqrt(plim) - -0.92734095767679814 *
      fcsModel_DW.DiscreteTransferFcn_states[tlim];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      // Outport: '<Root>/actuatorsCmds'
      fcsModel_Y.actuatorsCmds[tlim] = -1.0;
    } else {
      // Outport: '<Root>/actuatorsCmds' incorporates:
      //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

      fcsModel_Y.actuatorsCmds[tlim] = 0.036329521161600868 * plim +
        0.036329521161600868 * fcsModel_DW.DiscreteTransferFcn_states[tlim];
    }

    // End of Switch: '<S1>/Switch'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    DiscreteTransferFcn_tmp[tlim] = plim;
  }

  // RateTransition: '<Root>/Rate Transition'
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_Y.fcsDebug.outerLoopCtrlDebug = fcsModel_DW.RateTransition_Buffer0;

    // Switch: '<S3>/Switch1' incorporates:
    //   RateTransition: '<Root>/Rate Transition'

    if (rtb_atCenter) {
      // BusCreator: '<S3>/Bus Creator'
      rtb_BusCreator_b_frcCmd_N = DiscreteTransferFcn;
    } else {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S149>/Bus Assignment'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    }

    // End of Switch: '<S3>/Switch1'

    // BusCreator: '<S3>/Bus Creator' incorporates:
    //   BusCreator: '<S97>/Bus Creator'
    //   BusCreator: '<S98>/Bus Creator'
    //   ForEachSliceAssignment generated from: '<S101>/cmd'
    //   ForEachSliceAssignment generated from: '<S150>/filtCmd'

    rtb_BusCreator_b_velCtrlDebug_c = rtb_ImpAsg_InsertedFor_filtCmd_[0];
    rtb_BusCreator_b_velCtrlDebug_m = rtb_ImpAsg_InsertedFor_filtMeas[0];
    rtb_BusCreator_b_velCtrlDebug_p = rtb_ImpAsg_InsertedFor_pidDebug[0];
    rtb_BusCreator_b_posCtrlDebug_c = rtb_ImpAsg_InsertedFor_cmd_at_i[0];
    rtb_BusCreator_b_posCtrlDebug_m = rtb_ImpAsg_InsertedFor_meas_at_[0];
    rtb_BusCreator_b_posCtrlDebug_p = rtb_ImpAsg_InsertedFor_pidDeb_m[0];
    rtb_BusCreator_b_velCtrlDebug_0 = rtb_ImpAsg_InsertedFor_filtCmd_[1];
    rtb_BusCreator_b_velCtrlDebug_1 = rtb_ImpAsg_InsertedFor_filtMeas[1];
    rtb_BusCreator_b_velCtrlDebug_2 = rtb_ImpAsg_InsertedFor_pidDebug[1];
    rtb_BusCreator_b_posCtrlDebug_0 = rtb_ImpAsg_InsertedFor_cmd_at_i[1];
    rtb_BusCreator_b_posCtrlDebug_1 = rtb_ImpAsg_InsertedFor_meas_at_[1];
    rtb_BusCreator_b_posCtrlDebug_2 = rtb_ImpAsg_InsertedFor_pidDeb_m[1];
    rtb_BusCreator_b_velCtrlDebug_3 = rtb_ImpAsg_InsertedFor_filtCmd_[2];
    rtb_BusCreator_b_velCtrlDebug_4 = rtb_ImpAsg_InsertedFor_filtMeas[2];
    rtb_BusCreator_b_velCtrlDebug_5 = rtb_ImpAsg_InsertedFor_pidDebug[2];
    rtb_BusCreator_b_posCtrlDebug_3 = rtb_ImpAsg_InsertedFor_cmd_at_i[2];
    rtb_BusCreator_b_posCtrlDebug_4 = rtb_ImpAsg_InsertedFor_meas_at_[2];
    rtb_BusCreator_b_posCtrlDebug_5 = rtb_ImpAsg_InsertedFor_pidDeb_m[2];
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
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.cmd[0] =
      rtb_BusCreator_b_posCtrlDebug_c;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.meas[0] =
      rtb_BusCreator_b_posCtrlDebug_m;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.pidDebug[0] =
      rtb_BusCreator_b_posCtrlDebug_p;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[1] =
      rtb_BusCreator_b_velCtrlDebug_0;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[1] =
      rtb_BusCreator_b_velCtrlDebug_1;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[1] =
      rtb_BusCreator_b_velCtrlDebug_2;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.cmd[1] =
      rtb_BusCreator_b_posCtrlDebug_0;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.meas[1] =
      rtb_BusCreator_b_posCtrlDebug_1;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.pidDebug[1] =
      rtb_BusCreator_b_posCtrlDebug_2;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[2] =
      rtb_BusCreator_b_velCtrlDebug_3;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[2] =
      rtb_BusCreator_b_velCtrlDebug_4;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[2] =
      rtb_BusCreator_b_velCtrlDebug_5;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.cmd[2] =
      rtb_BusCreator_b_posCtrlDebug_3;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.meas[2] =
      rtb_BusCreator_b_posCtrlDebug_4;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.pidDebug[2] =
      rtb_BusCreator_b_posCtrlDebug_5;
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
    int32_T ForEach_itr_i;
    int32_T ForEach_itr_h;
    int32_T ForEach_itr_p;
    std::array<real_T, 3> den;
    std::array<real_T, 3> num;

    // SystemInitialize for Atomic SubSystem: '<S100>/baroAltEstimate'
    // SystemInitialize for MATLAB Function: '<S106>/Compute Filter Numerator And Denominator' 
    ComputeFilterNumeratorAndD_Init(&num[0], &den[0]);

    // End of SystemInitialize for SubSystem: '<S100>/baroAltEstimate'
    // 'computeAltFromBaro_function:4' pers_base_alt_point_m = 0;
    // 'holdOutputAtCenter_function:6' last_input = 0;
    // SystemInitialize for Iterator SubSystem: '<S97>/NED Position Control'
    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // SystemInitialize for Iterator SubSystem: '<S97>/NED Position Control'
      // SystemInitialize for Atomic SubSystem: '<S101>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S101>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S101>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S101>/Signal Conditioning Block1' 
      // End of SystemInitialize for SubSystem: '<S97>/NED Position Control'
    }

    // End of SystemInitialize for SubSystem: '<S97>/NED Position Control'
    // SystemInitialize for Iterator SubSystem: '<S98>/For Each Subsystem'
    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // SystemInitialize for Iterator SubSystem: '<S98>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S150>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S150>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S150>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S150>/Signal Conditioning Block1' 
      // End of SystemInitialize for SubSystem: '<S98>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S98>/For Each Subsystem'
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
