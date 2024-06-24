//
// File: fcsModel.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.103
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Fri Jun 14 19:31:46 2024
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
#include "fcsModel_private.h"

// Named constants for Chart: '<S4>/Chart'
const uint8_T fcsModel_IN_ARM_MTRS{ 1U };

const uint8_T fcsModel_IN_INACTIVE{ 2U };

const uint8_T fcsModel_IN_INFLIGHT{ 3U };

const busXyBodyAccelCtrIDebug fcsModel_rtZbusXyBodyAccelCtrIDebug{
  { {
      0.0, 0.0 } }
  ,                                    // cmd

  { {
      0.0, 0.0 } }
  ,                                    // meas

  { {
      {
        0.0,                           // output
        0.0,                           // proportionalOutput
        0.0,                           // integralOutput
        0.0                            // derivativeOutput
      }, {
        0.0,                           // output
        0.0,                           // proportionalOutput
        0.0,                           // integralOutput
        0.0                            // derivativeOutput
      } } }
  // pidDebug
} ;                                    // busXyBodyAccelCtrIDebug ground

static void rate_scheduler(fcsModel::RT_MODEL_fcsModel_T *const fcsModel_M);
const busFcsDebug fcsModel_rtZbusFcsDebug{ { 0.0,// thrustCmd_N
    0.0,                               // xMomCmd_Nm
    0.0,                               // yMomCmd_Nm
    0.0                                // zMomCmd_Nm
  },                                   // allocDebug

  { { { { 0.0, 0.0, 0.0 } },           // cmd

      { { 0.0, 0.0, 0.0 } },           // meas

      { { { 0.0,                       // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          } } }                        // pidDebug
    },                                 // angRateCtrlDebug

    { { { 0.0, 0.0, 0.0 } },           // cmd

      { { 0.0, 0.0, 0.0 } },           // meas

      { { { 0.0,                       // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          } } }                        // pidDebug
    }                                  // attCtrlDebug
  },                                   // innerLoopCtrlDebug

  { 0.0,                               // frcCmd_N
    { { { 0.0, 0.0, 0.0 } },           // cmd

      { { 0.0, 0.0, 0.0 } },           // meas

      { { 0.0, 0.0, 0.0 } },           // velCtrlFf

      { { { 0.0,                       // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          } } }                        // pidDebug
    },                                 // velCtrlDebug

    { { { 0.0, 0.0, 0.0 } },           // cmd

      { { 0.0, 0.0, 0.0 } },           // meas

      { { { 0.0,                       // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          } } }                        // pidDebug
    },                                 // posCtrlDebug

    { 0.0,                             // cmd
      0.0,                             // meas

      { 0.0,                           // output
        0.0,                           // proportionalOutput
        0.0,                           // integralOutput
        0.0                            // derivativeOutput
      }                                // pidDebug
    },                                 // zAccelCtrlDebug

    { { { 0.0, 0.0 } },                // cmd

      { { 0.0, 0.0 } },                // meas

      { { { 0.0,                       // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          }, { 0.0,                    // output
            0.0,                       // proportionalOutput
            0.0,                       // integralOutput
            0.0                        // derivativeOutput
          } } }                        // pidDebug
    }                                  // xyBodyAccelCtrlDebug
  },                                   // outerLoopCtrlDebug
  enumStateMachine::INACTIVE,          // state
  enumFlightMode::STABILIZE            // flightMode
};

uint32_T plook_bincpag(real_T u, const real_T bp[], uint32_T maxIndex, real_T
  *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  // Prelookup - Index and Fraction
  // Index Search method: 'binary'
  // Use previous index: 'on'
  // Use last breakpoint for index at or above upper limit: 'on'
  // Remove protection against out-of-range input in generated code: 'on'

  if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d_prevIdx(u, bp, *prevIndex, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex;
    *fraction = 0.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp1d_la(uint32_T bpIndex, real_T frac, const real_T table[], uint32_T
                  maxIndex)
{
  real_T y;

  // Column-major Interpolation 1-D
  // Interpolation method: 'Linear point-slope'
  // Use last breakpoint for index at or above upper limit: 'on'
  // Overflow mode: 'wrapping'

  if (bpIndex == maxIndex) {
    y = table[bpIndex];
  } else {
    real_T yL_0d0;
    yL_0d0 = table[bpIndex];
    y = (table[bpIndex + 1U] - yL_0d0) * frac + yL_0d0;
  }

  return y;
}

uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T found;
  uint32_T iLeft;
  uint32_T iRght;

  // Binary Search using Previous Index
  bpIndex = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIndex]) {
      iRght = bpIndex - 1U;
      bpIndex = ((bpIndex + iLeft) - 1U) >> 1U;
    } else if (u < bp[bpIndex + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIndex + 1U;
      bpIndex = ((bpIndex + iRght) + 1U) >> 1U;
    }
  }

  return bpIndex;
}

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
//    '<S14>/Discrete First Order Deriv Filter'
//    '<S58>/Discrete First Order Deriv Filter'
//    '<S120>/Discrete First Order Deriv Filter'
//    '<S174>/Discrete First Order Deriv Filter'
//    '<S213>/Discrete First Order Deriv Filter'
//
void fcsModel::f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
  rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
  rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW)
{
  real_T K;
  real_T normalizer;
  real_T num_tmp;

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

  // DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn'
  K = rtu_input - localDW->den[1] * localDW->DiscreteTransferFcn_states;
  *rty_filteredInputRate = localDW->num[0] * K + localDW->num[1] *
    localDW->DiscreteTransferFcn_states;

  // Update for DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states = K;
}

//
// System initialize for atomic system:
//    '<S11>/pidWithDebug'
//    '<S52>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_Init(DW_pidWithDebug_fcsModel_T *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S14>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S11>/pidWithDebug'
//    '<S52>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc, const
  busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
  busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s, DW_pidWithDebug_fcsModel_T
  *localDW)
{
  real_T rtb_Product5_o;
  real_T rtb_Sum1_b;
  real_T rtb_Sum_k;
  real_T rtb_Switch2_oz;
  real_T rtb_Switch2_p;
  real_T rtb_UkYk1_k;
  real_T rtb_UnitDelay_i;

  // Product: '<S46>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S46>/sample time'
  //
  //  About '<S46>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[1] * 0.004;

  // Sum: '<S14>/Sum'
  rtb_Sum_k = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S14>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_k,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_o, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S14>/Discrete First Order Deriv Filter'

  // Product: '<S14>/Product'
  rtb_Product5_o *= rtu_pidParamBus->Kd;

  // Product: '<S14>/Product1'
  rtb_UnitDelay_i = rtb_Sum_k * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S14>/Discrete-Time Integrator'
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S14>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S14>/Discrete-Time Integrator'

  rtb_Sum1_b = ((rtu_feedForward + rtb_Product5_o) + rtb_UnitDelay_i) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S47>/Switch2' incorporates:
  //   RelationalOperator: '<S47>/LowerRelop1'
  //   RelationalOperator: '<S47>/UpperRelop'
  //   Switch: '<S47>/Switch'

  if (rtb_Sum1_b > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_b < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S47>/Switch'
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_oz = rtb_Sum1_b;
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

  rtb_UkYk1_k = rtb_Switch2_oz - localDW->DelayInput2_DSTATE;

  // Switch: '<S49>/Switch2' incorporates:
  //   RelationalOperator: '<S49>/LowerRelop1'

  if (rtb_UkYk1_k <= rtb_Switch2_p) {
    // Product: '<S46>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S46>/sample time'
    //
    //  About '<S46>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S49>/Switch' incorporates:
    //   RelationalOperator: '<S49>/UpperRelop'

    if (rtb_UkYk1_k >= rtb_Switch2_p) {
      rtb_Switch2_p = rtb_UkYk1_k;
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

  *rty_ctrlCmd = rtb_Switch2_p + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S14>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S14>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_i;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_o;

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

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_k * rtu_pidParamBus->Ki) * 0.004;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S46>/Delay Input2'
  //
  //  Block description for '<S46>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S14>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_oz;

  // Update for UnitDelay: '<S14>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_b;
}

//
// Output and update for atomic system:
//    '<S30>/Compute Natural Frequency'
//    '<S31>/Compute Natural Frequency'
//    '<S15>/Compute Natural Frequency'
//    '<S16>/Compute Natural Frequency'
//    '<S74>/Compute Natural Frequency'
//    '<S75>/Compute Natural Frequency'
//    '<S59>/Compute Natural Frequency'
//    '<S60>/Compute Natural Frequency'
//    '<S136>/Compute Natural Frequency'
//    '<S137>/Compute Natural Frequency'
//    ...
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
//    '<S30>/Compute Numerator And Denominator'
//    '<S15>/Compute Numerator And Denominator'
//    '<S74>/Compute Numerator And Denominator'
//    '<S59>/Compute Numerator And Denominator'
//    '<S136>/Compute Numerator And Denominator'
//    '<S121>/Compute Numerator And Denominator'
//    '<S175>/Compute Numerator And Denominator'
//    '<S190>/Compute Numerator And Denominator'
//    '<S244>/Compute Numerator And Denominator'
//    '<S229>/Compute Numerator And Denominator'
//    '<S214>/Compute Numerator And Denominator'
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
//    '<S31>/Compute Filter Numerator And Denominator'
//    '<S16>/Compute Filter Numerator And Denominator'
//    '<S75>/Compute Filter Numerator And Denominator'
//    '<S60>/Compute Filter Numerator And Denominator'
//    '<S137>/Compute Filter Numerator And Denominator'
//    '<S122>/Compute Filter Numerator And Denominator'
//    '<S176>/Compute Filter Numerator And Denominator'
//    '<S191>/Compute Filter Numerator And Denominator'
//    '<S245>/Compute Filter Numerator And Denominator'
//    '<S230>/Compute Filter Numerator And Denominator'
//    '<S215>/Compute Filter Numerator And Denominator'
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
//    '<S31>/Compute Filter Numerator And Denominator'
//    '<S16>/Compute Filter Numerator And Denominator'
//    '<S75>/Compute Filter Numerator And Denominator'
//    '<S60>/Compute Filter Numerator And Denominator'
//    '<S137>/Compute Filter Numerator And Denominator'
//    '<S122>/Compute Filter Numerator And Denominator'
//    '<S176>/Compute Filter Numerator And Denominator'
//    '<S191>/Compute Filter Numerator And Denominator'
//    '<S245>/Compute Filter Numerator And Denominator'
//    '<S230>/Compute Filter Numerator And Denominator'
//    '<S215>/Compute Filter Numerator And Denominator'
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
//    '<S11>/Signal Conditioning Block1'
//    '<S11>/Signal Conditioning Block'
//    '<S52>/Signal Conditioning Block1'
//    '<S52>/Signal Conditioning Block'
//
void fcsModel::f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S31>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
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
  std::array<real_T, 3> rtb_accelNum_b;
  std::array<real_T, 3> rtb_den;
  std::array<real_T, 3> rtb_rateNum;
  real_T rtb_DiscreteTransferFcn_j;
  real_T rtb_Switch2_h;

  // MATLAB Function: '<S30>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S30>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum_b[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S31>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S31>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S31>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_j = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S35>/Switch2' incorporates:
  //   RelationalOperator: '<S35>/LowerRelop1'
  //   RelationalOperator: '<S35>/UpperRelop'
  //   Switch: '<S35>/Switch'

  if (rtb_DiscreteTransferFcn_j > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_j < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S35>/Switch'
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[0];
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

  rtb_DiscreteTransferFcn_j -= localDW->DelayInput2_DSTATE;

  // Switch: '<S42>/Switch2' incorporates:
  //   Product: '<S32>/delta rise limit'
  //   SampleTimeMath: '<S32>/sample time'
  //
  //  About '<S32>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_h = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S42>/Switch2' incorporates:
  //   RelationalOperator: '<S42>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_j <= rtb_Switch2_h) {
    // Product: '<S32>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S32>/sample time'
    //
    //  About '<S32>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_h = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S42>/Switch' incorporates:
    //   RelationalOperator: '<S42>/UpperRelop'

    if (rtb_DiscreteTransferFcn_j >= rtb_Switch2_h) {
      // Switch: '<S42>/Switch2'
      rtb_Switch2_h = rtb_DiscreteTransferFcn_j;
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

  *rty_filteredInput = rtb_Switch2_h + localDW->DelayInput2_DSTATE;

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
// Output and update for atomic system:
//    '<S103>/holdOutputAtCenter1'
//    '<S103>/holdOutputAtCenter2'
//
void fcsModel::fcsModel_holdOutputAtCenter1(real_T rtu_input, real_T rtu_trigger,
  real_T *rty_output, boolean_T *rty_atCenter, DW_holdOutputAtCenter1_fcsMod_T
  *localDW)
{
  // MATLAB Function: '<S113>/holdOutputAtCenter'
  // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S116>:1'
  // '<S116>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
  // HOLDOUTPUTATCENTER_FUNCTION holds the output constant at last input if the
  // trigger value is within user defined delta from the center
  // 'holdOutputAtCenter_function:5' if isempty(last_input)
  // 'holdOutputAtCenter_function:9' if(trigger <= (params.center + params.posDeltaFromCenter) && ... 
  // 'holdOutputAtCenter_function:10'         trigger >=(params.center - params.negDeltaFromCenter)) 
  if ((rtu_trigger <= 0.05) && (rtu_trigger >= -0.05)) {
    // 'holdOutputAtCenter_function:11' atCenter = true;
    *rty_atCenter = true;
  } else {
    // 'holdOutputAtCenter_function:12' else
    // 'holdOutputAtCenter_function:13' atCenter = false;
    *rty_atCenter = false;

    // 'holdOutputAtCenter_function:14' last_input = input;
    localDW->last_input = rtu_input;
  }

  // 'holdOutputAtCenter_function:17' output = last_input;
  *rty_output = localDW->last_input;

  // End of MATLAB Function: '<S113>/holdOutputAtCenter'
}

//
// System initialize for atomic system:
//    '<S104>/pidWithDebug'
//    '<S171>/pidWithDebug'
//    '<S157>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_m_Init(DW_pidWithDebug_fcsModel_i_T
  *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S120>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S104>/pidWithDebug'
//    '<S171>/pidWithDebug'
//    '<S157>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc, const
  busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
  busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
  DW_pidWithDebug_fcsModel_i_T *localDW)
{
  real_T rtb_Product5_e;
  real_T rtb_Sum1_o;
  real_T rtb_Sum_c4;
  real_T rtb_Switch2_b;
  real_T rtb_Switch2_d3;
  real_T rtb_UkYk1_h;
  real_T rtb_UnitDelay_a;

  // Product: '<S152>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S152>/sample time'
  //
  //  About '<S152>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_d3 = rtu_pidParamBus->outputRateLimits[1] * 0.02;

  // Sum: '<S120>/Sum'
  rtb_Sum_c4 = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S120>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_c4,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_e, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S120>/Discrete First Order Deriv Filter'

  // Product: '<S120>/Product'
  rtb_Product5_e *= rtu_pidParamBus->Kd;

  // Product: '<S120>/Product1'
  rtb_UnitDelay_a = rtb_Sum_c4 * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S120>/Discrete-Time Integrator'
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S120>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S120>/Discrete-Time Integrator'

  rtb_Sum1_o = ((rtu_feedForward + rtb_Product5_e) + rtb_UnitDelay_a) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S153>/Switch2' incorporates:
  //   RelationalOperator: '<S153>/LowerRelop1'
  //   RelationalOperator: '<S153>/UpperRelop'
  //   Switch: '<S153>/Switch'

  if (rtb_Sum1_o > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_b = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_o < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S153>/Switch'
    rtb_Switch2_b = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_b = rtb_Sum1_o;
  }

  // End of Switch: '<S153>/Switch2'

  // Sum: '<S152>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S152>/Delay Input2'
  //
  //  Block description for '<S152>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S152>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1_h = rtb_Switch2_b - localDW->DelayInput2_DSTATE;

  // Switch: '<S155>/Switch2' incorporates:
  //   RelationalOperator: '<S155>/LowerRelop1'

  if (rtb_UkYk1_h <= rtb_Switch2_d3) {
    // Product: '<S152>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S152>/sample time'
    //
    //  About '<S152>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_d3 = rtu_pidParamBus->outputRateLimits[0] * 0.02;

    // Switch: '<S155>/Switch' incorporates:
    //   RelationalOperator: '<S155>/UpperRelop'

    if (rtb_UkYk1_h >= rtb_Switch2_d3) {
      rtb_Switch2_d3 = rtb_UkYk1_h;
    }

    // End of Switch: '<S155>/Switch'
  }

  // End of Switch: '<S155>/Switch2'

  // Sum: '<S152>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S152>/Delay Input2'
  //
  //  Block description for '<S152>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S152>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Switch2_d3 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S120>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S120>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_a;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_e;

  // Update for DiscreteIntegrator: '<S120>/Discrete-Time Integrator' incorporates:
  //   Product: '<S120>/Product2'
  //   Product: '<S120>/Product3'
  //   Product: '<S120>/Product5'
  //   Sum: '<S120>/Sum2'
  //   Sum: '<S120>/Sum3'
  //   Sum: '<S120>/Sum4'
  //   Sum: '<S120>/Sum5'
  //   UnitDelay: '<S120>/Unit Delay'
  //   UnitDelay: '<S120>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_c4 * rtu_pidParamBus->Ki) * 0.02;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S152>/Delay Input2'
  //
  //  Block description for '<S152>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S120>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_b;

  // Update for UnitDelay: '<S120>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_o;
}

//
// System initialize for atomic system:
//    '<S104>/Signal Conditioning Block1'
//    '<S104>/Signal Conditioning Block'
//    '<S171>/Signal Conditioning Block'
//    '<S157>/Signal Conditioning Block2'
//    '<S157>/Signal Conditioning Block1'
//    '<S157>/Signal Conditioning Block'
//
void fcsModel::SignalConditioningBlock1_c_Init(DW_SignalConditioningBlock1_g_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S137>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S104>/Signal Conditioning Block1'
//    '<S104>/Signal Conditioning Block'
//    '<S171>/Signal Conditioning Block'
//    '<S157>/Signal Conditioning Block2'
//    '<S157>/Signal Conditioning Block1'
//    '<S157>/Signal Conditioning Block'
//
void fcsModel::fcsM_SignalConditioningBlock1_f(real_T rtu_input, const
  busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
  rtp_sampleTime_s, DW_SignalConditioningBlock1_g_T *localDW)
{
  std::array<real_T, 3> rtb_accelNum_k;
  std::array<real_T, 3> rtb_den;
  std::array<real_T, 3> rtb_rateNum;
  real_T rtb_DiscreteTransferFcn_d;
  real_T rtb_Switch2_n;

  // MATLAB Function: '<S136>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_n);

  // MATLAB Function: '<S136>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2_n,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum_k[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S137>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_n);

  // MATLAB Function: '<S137>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2_n,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S137>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_d = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S141>/Switch2' incorporates:
  //   RelationalOperator: '<S141>/LowerRelop1'
  //   RelationalOperator: '<S141>/UpperRelop'
  //   Switch: '<S141>/Switch'

  if (rtb_DiscreteTransferFcn_d > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_d < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S141>/Switch'
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S141>/Switch2'

  // Sum: '<S138>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S138>/Delay Input2'
  //
  //  Block description for '<S138>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S138>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_d -= localDW->DelayInput2_DSTATE;

  // Switch: '<S148>/Switch2' incorporates:
  //   Product: '<S138>/delta rise limit'
  //   SampleTimeMath: '<S138>/sample time'
  //
  //  About '<S138>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_n = rtu_params->filteredInputRateLimits[1] * 0.02;

  // Switch: '<S148>/Switch2' incorporates:
  //   RelationalOperator: '<S148>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_d <= rtb_Switch2_n) {
    // Product: '<S138>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S138>/sample time'
    //
    //  About '<S138>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_n = rtu_params->filteredInputRateLimits[0] * 0.02;

    // Switch: '<S148>/Switch' incorporates:
    //   RelationalOperator: '<S148>/UpperRelop'

    if (rtb_DiscreteTransferFcn_d >= rtb_Switch2_n) {
      // Switch: '<S148>/Switch2'
      rtb_Switch2_n = rtb_DiscreteTransferFcn_d;
    }

    // End of Switch: '<S148>/Switch'
  }

  // End of Switch: '<S148>/Switch2'

  // Sum: '<S138>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S138>/Delay Input2'
  //
  //  Block description for '<S138>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S138>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2_n + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S137>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S138>/Delay Input2'
  //
  //  Block description for '<S138>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_filteredInput;
}

//
// Function for Chart: '<S4>/Chart'
// function isTrue = checkRcCmds
//
boolean_T fcsModel::fcsModel_checkRcCmds(const busRcInCmds
  *BusConversion_InsertedFor_Chart)
{
  boolean_T isTrue;

  // MATLAB Function 'checkRcCmds': '<S264>:7'
  // '<S264>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S264>:7:3' if(rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle(1) && ... 
  // '<S264>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S264>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S264>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 1000) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 1000) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 1000) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 1000) {
          // '<S264>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S264>:7:8' else
          // '<S264>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S264>:7:8' else
        // '<S264>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S264>:7:8' else
      // '<S264>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S264>:7:8' else
    // '<S264>:7:9' isTrue = false;
    isTrue = false;
  }

  return isTrue;
}

// Model step function
void fcsModel::step()
{
  // local scratch DWork variables
  int32_T ForEach_itr_l;
  int32_T ForEach_itr_p;
  int32_T ForEach_itr_g;
  int32_T ForEach_itr;
  int32_T ForEach_itr_i;
  std::array<real_T, 4> DiscreteTransferFcn_tmp;
  std::array<real_T, 4> rtb_DiscreteTransferFcn_e;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angRateC;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_cmd_at_i;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtCmd_;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtMeas;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_meas_at_;
  std::array<real_T, 4> rtb_ImpAsg_InsertedFor_mtrPwmCm;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_neVelCmd;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDeb_m;
  std::array<busPidDebug, 3> rtb_ImpAsg_InsertedFor_pidDebug;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_velCtrlF;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_velCtrlO;
  std::array<real_T, 3> rtb_MatrixMultiply;
  std::array<real_T, 9> rtb_Transpose;
  std::array<busCtrlInputs, 3> rtb_VectorConcatenate;
  std::array<real_T, 3> rtb_VectorConcatenate1;
  busCtrlInputs rtb_BusAssignment1_a;
  busCtrlInputs rtb_BusAssignment2;
  busCtrlInputs rtb_BusAssignment4;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_2;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_5;
  busPidDebug rtb_BusCreator_b_posCtrlDebug_p;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_3;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_7;
  busPidDebug rtb_BusCreator_b_velCtrlDebug_p;
  busPidDebug rtb_BusCreator_em;
  busPidDebug rtb_BusCreator_f;
  busXyBodyAccelCtrIDebug rtb_BusCreator_b_xyBodyAccelCtr;
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
  real_T rtb_BusCreator_b_velCtrlDebug_2;
  real_T rtb_BusCreator_b_velCtrlDebug_4;
  real_T rtb_BusCreator_b_velCtrlDebug_5;
  real_T rtb_BusCreator_b_velCtrlDebug_6;
  real_T rtb_BusCreator_b_velCtrlDebug_c;
  real_T rtb_BusCreator_b_velCtrlDebug_m;
  real_T rtb_BusCreator_b_velCtrlDebug_v;
  real_T rtb_BusCreator_b_zAccelCtrlDe_0;
  real_T rtb_BusCreator_b_zAccelCtrlDe_1;
  real_T rtb_BusCreator_b_zAccelCtrlDe_2;
  real_T rtb_BusCreator_b_zAccelCtrlDe_3;
  real_T rtb_BusCreator_b_zAccelCtrlDe_4;
  real_T rtb_BusCreator_b_zAccelCtrlDebu;
  real_T rtb_DifferenceInputs2_a;
  real_T rtb_DifferenceInputs2_p;
  real_T rtb_DiscreteTransferFcn;
  real_T rtb_DiscreteTransferFcn_k;
  real_T rtb_frcCmd_N;
  real_T vxCmd_unitRange;
  real_T yCmd;
  real_T ylim;
  int32_T pCmd;
  int32_T rCmd;
  int32_T tCmd;
  boolean_T resetIntegrator;
  boolean_T rtb_Compare_d;
  boolean_T rtb_Compare_mi;
  enumFlightMode flightMode;
  enumStateMachine state;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // Math: '<S103>/Transpose' incorporates:
    //   Inport: '<Root>/stateEstimate'

    tCmd = 0;
    for (rCmd = 0; rCmd < 3; rCmd++) {
      rtb_Transpose[tCmd] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd];
      rtb_Transpose[tCmd + 1] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd + 3];
      rtb_Transpose[tCmd + 2] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd + 6];
      tCmd += 3;
    }

    // End of Math: '<S103>/Transpose'
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
    // Transition: '<S264>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S264>:1'
    // '<S264>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S264>:1:3' rcCheckFlag = checkRcCmds;
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
    }

    // '<S264>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S264>:3'
      // '<S264>:10:1' sf_internal_predicateOutput = after(60, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 15000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S264>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S264>:1'
        // '<S264>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S264>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S264>:1:4' resetIntegrator = true;

        // '<S264>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimitsThrottle(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 1000) {
        // Transition: '<S264>:12'
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S264>:11'
        // '<S264>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S264>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S264>:11:4' resetIntegrator = false;
      } else {
        // '<S264>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S264>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S264>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S264>:1'
      // '<S264>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 1 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S264>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S264>:3'
        // '<S264>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S264>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S264>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S264>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S264>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S264>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S264>:11'
      // '<S264>:20:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle (1); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1000) {
        // Transition: '<S264>:20'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S264>:3'
        // '<S264>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S264>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S264>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S264>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S264>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S264>:11:4' resetIntegrator = false;
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
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S265>:1'
  // '<S265>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcCmds, rcParamsStruct); 
  // INTERPRETRCINPUTS_FUNCTION
  // Computes command and flight mode from the rc inputs
  //  In Alt hold mode or Pos hold mode throttle bottom position is
  //  mapped to max descent rate which then decreases to zero when throttle is
  //  mid-stick on the RC transmitter and from there as throttle is increased
  //  throttle top position is mapped to max asent rate. This works well when
  //  vehicle is already in the air. For special case when we start the Alt
  //  hold mode or Pos hold from takeoff having bottom half of throttle mapped
  //  to descent rate windsup the climbrate controller integrator as the
  //  vehicle is on the ground with zero vertical velocity and vehicle is being
  //  commanded to descent. Also having the bottom half of the throttle mapped
  //  to descent rate means pilot has to move the stick to half throttle before
  //  vehicle is commanded to ascent during takeoff. A simple solution to fix this special 
  //  case is to have a flag that is not set when throttle is at bottom and set
  //  when the throttle reaches the mid stick. This way the flag can be used to
  //  set entire throttle stick range to be climb rate when we are taking off
  //  and once we reach mid throttle the throttle range will be split into
  //  equal parts descent and ascent rate
  // 'interpretRcInputs_function:23' if isempty(throttle_is_up)
  //  Used to directly set force commands in N that is fed into allocation in
  //  STABILIZE and ACRO flight modes
  // 'interpretRcInputs_function:29' rcOutCmds.throttleStick = 0;
  //  Used to command roll and pitch angles in STABILIZE flight mode
  // 'interpretRcInputs_function:32' rcOutCmds.rollStick = 0;
  // 'interpretRcInputs_function:33' rcOutCmds.pitchStick = 0;
  //  Used to control yaw rate in all flight modes
  // 'interpretRcInputs_function:36' rcOutCmds.yawStick = 0;
  //  Not used currently
  // 'interpretRcInputs_function:39' rcOutCmds.vxStick_mps = 0;
  // 'interpretRcInputs_function:40' rcOutCmds.vyStick_mps = 0;
  //  Used to command vertical velocity in ALT_CONTROL flight mode. In this
  //  mode rcOutCmds.throttleStick is ignored.
  // 'interpretRcInputs_function:44' rcOutCmds.vzStick_mps = 0;
  // Select mode and set the command limits
  // 'interpretRcInputs_function:47' if(rcInCmds.rcSwitch1_nd < 1100)
  if (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1100) {
    // 'interpretRcInputs_function:48' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:49' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:50' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:51' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:52' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:53' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
  } else if ((fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1100) &&
             (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1700)) {
    // 'interpretRcInputs_function:55' elseif (rcInCmds.rcSwitch1_nd >= 1100 && rcInCmds.rcSwitch1_nd < 1700) 
    // 'interpretRcInputs_function:56' flightMode = enumFlightMode.ALT_CONTROL;
    flightMode = enumFlightMode::ALT_CONTROL;

    // 'interpretRcInputs_function:57' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:58' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:59' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:60' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:61' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else if (fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1700) {
    // 'interpretRcInputs_function:62' elseif (rcInCmds.rcSwitch1_nd >= 1700)
    // 'interpretRcInputs_function:63' flightMode = enumFlightMode.POS_CONTROL;
    flightMode = enumFlightMode::POS_CONTROL;

    // 'interpretRcInputs_function:64' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:65' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:66' rlim = 0;
    rlim = 0.0;

    // 'interpretRcInputs_function:67' plim = 0;
    plim = 0.0;

    // 'interpretRcInputs_function:68' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else {
    // 'interpretRcInputs_function:69' else
    // 'interpretRcInputs_function:70' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:71' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:72' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:73' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:74' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:75' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
  }

  // 'interpretRcInputs_function:78' if double(rcInCmds.throttleCmd_nd) < rcParamsStruct.pwmLimits(1) 
  if (fcsModel_U.rcCmdsIn.throttleCmd_nd < 1000) {
    // This means either we haven't take off yet or we landed after a flight
    // and might take of again so set the throttle_is_up to false
    // 'interpretRcInputs_function:81' throttle_is_up = false;
    fcsModel_DW.throttle_is_up = false;
  } else if ((fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1450) &&
             (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 1300)) {
    // 'interpretRcInputs_function:82' elseif ((double(rcInCmds.throttleCmd_nd) <= rcParamsStruct.pwmThrottleMidHigh) && ... 
    // 'interpretRcInputs_function:83'         double(rcInCmds.throttleCmd_nd) >= rcParamsStruct.pwmThrottleMidLow) 
    // Since throttle stick on RC transmitter is half way up this probably
    // means the vehicle is flying so set throttle_is_up flag to true
    // 'interpretRcInputs_function:86' throttle_is_up = true;
    fcsModel_DW.throttle_is_up = true;
  } else {
    // 'interpretRcInputs_function:87' else
    // Otherwise preserve the previous flag
  }

  // 'interpretRcInputs_function:91' tCmd = min( rcParamsStruct.pwmLimitsThrottle(2), ... 
  // 'interpretRcInputs_function:92'         max( rcParamsStruct.pwmLimitsThrottle(1), double(rcInCmds.throttleCmd_nd) ) ); 
  tCmd = static_cast<int32_T>(std::fmin(1900.0, std::fmax(1000.0, static_cast<
    real_T>(fcsModel_U.rcCmdsIn.throttleCmd_nd))));

  //  In stabilize mode throttle stick starts at 0
  // 'interpretRcInputs_function:95' tCmd_unitRange = (rcParamsStruct.throttleUnitRangeMapCoeff.a*tCmd + ... 
  // 'interpretRcInputs_function:96'     rcParamsStruct.throttleUnitRangeMapCoeff.c)/rcParamsStruct.throttleUnitRangeMapCoeff.b; 
  //  Use throttle stick to set Vz commands. Always set but only used
  //  in ALT_CONTROL flight mode. Has different slopes about the center point
  //  as the center point is not always at 1500 which is PWM center
  // 'interpretRcInputs_function:101' if ((tCmd <= rcParamsStruct.pwmThrottleMidHigh) && ... 
  // 'interpretRcInputs_function:102'         tCmd >= rcParamsStruct.pwmThrottleMidLow) 
  if ((tCmd <= 1450) && (tCmd >= 1300)) {
    // 'interpretRcInputs_function:103' rcOutCmds.vzStick_mps = 0;
    fcsModel_DW.rcOutCmds.vzStick_mps = 0.0;
  } else if (tCmd < 1300) {
    // 'interpretRcInputs_function:104' elseif (tCmd < rcParamsStruct.pwmThrottleMidLow) 
    // 'interpretRcInputs_function:105' if (throttle_is_up)
    if (fcsModel_DW.throttle_is_up) {
      //  This means we are already flying, make lower half of throttle
      //  stick map to descent rates
      // 'interpretRcInputs_function:108' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:109'         rcParamsStruct.pwmToCmdThrottleIncptLow; 
      // ;
      // 'interpretRcInputs_function:110' rcOutCmds.vzStick_mps = (rcParamsStruct.vzLowRangeMapCoeff.a*vzCmd_unitRange + ... 
      // 'interpretRcInputs_function:111'                         rcParamsStruct.vzLowRangeMapCoeff.c)/ ... 
      // 'interpretRcInputs_function:112'                         rcParamsStruct.vzLowRangeMapCoeff.b; 
      fcsModel_DW.rcOutCmds.vzStick_mps = -(0.0033333333333333335 * static_cast<
        real_T>(tCmd) + -4.333333333333333);
    } else {
      // 'interpretRcInputs_function:113' else
      //  This means we haven't taken off yet or we landed and might take
      //  off again. Map the lower half of the throttle stick also to
      //  ascent rate as we want to climb up from the ground and only go
      //  upto 3/4 of the max ascent rate at mid point. Lowest position is
      //  throttle maps to -1 so add +1 make lowest throttle point 0 and
      //  increase from there on
      // 'interpretRcInputs_function:120' vzCmd_unitRange = ((rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:121'             rcParamsStruct.pwmToCmdThrottleIncptLow) + 1)*0.75; 
      // ;
      // 'interpretRcInputs_function:122' rcOutCmds.vzStick_mps = (rcParamsStruct.vzHighRangeMapCoeff.a*vzCmd_unitRange + ... 
      // 'interpretRcInputs_function:123'                         rcParamsStruct.vzHighRangeMapCoeff.c)/ ... 
      // 'interpretRcInputs_function:124'                         rcParamsStruct.vzHighRangeMapCoeff.b; 
      fcsModel_DW.rcOutCmds.vzStick_mps = ((0.0033333333333333335 * static_cast<
        real_T>(tCmd) + -4.333333333333333) + 1.0) * 0.75 * -1.5;
    }
  } else {
    // 'interpretRcInputs_function:126' else
    // 'interpretRcInputs_function:127' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeHigh*tCmd + ...; 
    // 'interpretRcInputs_function:128'         rcParamsStruct.pwmToCmdThrottleIncptHigh; 
    // ;
    // 'interpretRcInputs_function:129' rcOutCmds.vzStick_mps = (rcParamsStruct.vzHighRangeMapCoeff.a*vzCmd_unitRange + ... 
    // 'interpretRcInputs_function:130'                         rcParamsStruct.vzHighRangeMapCoeff.c)/ ... 
    // 'interpretRcInputs_function:131'                         rcParamsStruct.vzHighRangeMapCoeff.b; 
    fcsModel_DW.rcOutCmds.vzStick_mps = (0.0022222222222222222 *
      static_cast<real_T>(tCmd) + -3.2222222222222223) * -1.5;
  }

  //  Set roll, pitch and yaw stick
  // 'interpretRcInputs_function:135' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:136'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  rCmd = static_cast<int32_T>(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickXCmd_nd))));

  //  Use roll stick to set FEP Vy to be used for POS control mode
  // 'interpretRcInputs_function:139' if ((rCmd <= rcParamsStruct.pwmRollStickMidHigh) && ... 
  // 'interpretRcInputs_function:140'         rCmd >= rcParamsStruct.pwmRollStickMidLow) 
  if ((rCmd <= 1650) && (rCmd >= 1350)) {
    // 'interpretRcInputs_function:141' vyCmd_unitRange  = 0;
    ylim = 0.0;
  } else {
    // 'interpretRcInputs_function:142' else
    // 'interpretRcInputs_function:143' vyCmd_unitRange = -3 + rCmd/500;
    ylim = static_cast<real_T>(rCmd) / 500.0 + -3.0;
  }

  // 'interpretRcInputs_function:145' rCmd_unitRange = -3 + rCmd/500;
  // 'interpretRcInputs_function:148' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:149'         max( rcParamsStruct.pwmLimits(1),  double(rcInCmds.joystickYCmd_nd) ) ); 
  pCmd = static_cast<int32_T>(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickYCmd_nd))));

  //  Use pitch stick to set FEP Vx to be used for POS control mode
  // 'interpretRcInputs_function:152' if ((pCmd <= rcParamsStruct.pwmPitchStickMidHigh) && ... 
  // 'interpretRcInputs_function:153'         pCmd >= rcParamsStruct.pwmPitchStickMidLow) 
  if ((pCmd <= 1650) && (pCmd >= 1350)) {
    // 'interpretRcInputs_function:154' vxCmd_unitRange  = 0;
    vxCmd_unitRange = 0.0;
  } else {
    // 'interpretRcInputs_function:155' else
    // 'interpretRcInputs_function:156' vxCmd_unitRange = -3 + pCmd/500;
    vxCmd_unitRange = static_cast<real_T>(pCmd) / 500.0 + -3.0;
  }

  //  Reverse the pitch cmd
  // 'interpretRcInputs_function:160' pCmd_unitRange = -(-3 + pCmd/500);
  // 'interpretRcInputs_function:163' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:164'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  yCmd = std::fmin(2000.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickZCmd_nd)));

  //  Use yaw stick to also pick a Yaw angle in ALT or POS control mode
  // 'interpretRcInputs_function:167' if (flightMode == enumFlightMode.STABILIZE) 
  if (flightMode == enumFlightMode::STABILIZE) {
    // 'interpretRcInputs_function:168' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;

    // 'interpretRcInputs_function:169' else
    // 'interpretRcInputs_function:170' if ((yCmd <= rcParamsStruct.pwmYawStickMidHigh) && ... 
    // 'interpretRcInputs_function:171'         yCmd >= rcParamsStruct.pwmYawStickMidLow) 
  } else if ((yCmd <= 1600.0) && (yCmd >= 1400.0)) {
    // 'interpretRcInputs_function:172' yCmd_unitRange  = 0;
    yCmd = 0.0;
  } else {
    // 'interpretRcInputs_function:173' else
    // 'interpretRcInputs_function:174' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;
  }

  //  Usually expo is set in the Tx hence simply use a linear map here
  // 'interpretRcInputs_function:179' rcOutCmds.throttleStick = (rcParamsStruct.fCmdRangeMapCoeff.a*tCmd_unitRange + ... 
  // 'interpretRcInputs_function:180'                         rcParamsStruct.fCmdRangeMapCoeff.c)/ ... 
  // 'interpretRcInputs_function:181'                         rcParamsStruct.fCmdRangeMapCoeff.b; 
  fcsModel_DW.rcOutCmds.throttleStick = (static_cast<real_T>(tCmd) + -1000.0) /
    900.0 * -20.0 + -15.0;

  // 'interpretRcInputs_function:182' rcOutCmds.rollStick = rCmd_unitRange*rlim; 
  fcsModel_DW.rcOutCmds.rollStick = (static_cast<real_T>(rCmd) / 500.0 + -3.0) *
    rlim;

  // 'interpretRcInputs_function:183' rcOutCmds.pitchStick = pCmd_unitRange*plim; 
  fcsModel_DW.rcOutCmds.pitchStick = -(static_cast<real_T>(pCmd) / 500.0 + -3.0)
    * plim;

  // 'interpretRcInputs_function:184' rcOutCmds.yawStick = yCmd_unitRange*ylim;
  fcsModel_DW.rcOutCmds.yawStick = yCmd * 1.0471975511965976;

  // 'interpretRcInputs_function:185' rcOutCmds.vxStick_mps = vxCmd_unitRange*vxlim; 
  fcsModel_DW.rcOutCmds.vxStick_mps = vxCmd_unitRange * 2.5;

  // 'interpretRcInputs_function:186' rcOutCmds.vyStick_mps = vyCmd_unitRange*vylim; 
  fcsModel_DW.rcOutCmds.vyStick_mps = ylim * 2.5;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // SignalConversion generated from: '<S103>/Vector Concatenate1'
    rtb_VectorConcatenate1[0] = fcsModel_DW.rcOutCmds.vxStick_mps;

    // SignalConversion generated from: '<S103>/Vector Concatenate1'
    rtb_VectorConcatenate1[1] = fcsModel_DW.rcOutCmds.vyStick_mps;

    // Product: '<S103>/Matrix Multiply' incorporates:
    //   Math: '<S103>/Transpose'

    for (tCmd = 0; tCmd < 3; tCmd++) {
      rtb_MatrixMultiply[tCmd] = 0.0;
      rtb_MatrixMultiply[tCmd] += rtb_Transpose[tCmd] * rtb_VectorConcatenate1[0];
      rtb_MatrixMultiply[tCmd] += rtb_Transpose[tCmd + 3] *
        rtb_VectorConcatenate1[1];
    }

    // End of Product: '<S103>/Matrix Multiply'

    // Outputs for Atomic SubSystem: '<S103>/holdOutputAtCenter1'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[0],
      rtb_MatrixMultiply[0], &rtb_frcCmd_N, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter1);

    // End of Outputs for SubSystem: '<S103>/holdOutputAtCenter1'

    // Concatenate: '<S103>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // Switch: '<S103>/Switch1' incorporates:
    //   RelationalOperator: '<S107>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S103>/Bus Assignment1' incorporates:
      //   Concatenate: '<S103>/Vector Concatenate'
      //   Constant: '<S103>/Constant1'

      rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S103>/Bus Assignment1' incorporates:
      //   Concatenate: '<S103>/Vector Concatenate'

      rtb_VectorConcatenate[0].feedForwardCmd = rtb_MatrixMultiply[0];
    }

    // End of Switch: '<S103>/Switch1'

    // BusAssignment: '<S103>/Bus Assignment1' incorporates:
    //   Concatenate: '<S103>/Vector Concatenate'
    //   Constant: '<S110>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S103>/Logical Operator2'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S110>/Compare'

    rtb_VectorConcatenate[0].cmd = rtb_frcCmd_N;
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedPos_m[0];
    rtb_VectorConcatenate[0].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S103>/holdOutputAtCenter2'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[1],
      rtb_MatrixMultiply[1], &rtb_frcCmd_N, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter2);

    // End of Outputs for SubSystem: '<S103>/holdOutputAtCenter2'

    // Concatenate: '<S103>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // Switch: '<S103>/Switch2' incorporates:
    //   RelationalOperator: '<S108>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S103>/Bus Assignment2' incorporates:
      //   Concatenate: '<S103>/Vector Concatenate'
      //   Constant: '<S103>/Constant3'

      rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S103>/Bus Assignment2' incorporates:
      //   Concatenate: '<S103>/Vector Concatenate'

      rtb_VectorConcatenate[1].feedForwardCmd = rtb_MatrixMultiply[1];
    }

    // End of Switch: '<S103>/Switch2'

    // BusAssignment: '<S103>/Bus Assignment2' incorporates:
    //   Concatenate: '<S103>/Vector Concatenate'
    //   Constant: '<S111>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S103>/Logical Operator3'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S111>/Compare'

    rtb_VectorConcatenate[1].cmd = rtb_frcCmd_N;
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedPos_m[1];
    rtb_VectorConcatenate[1].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S103>/holdOutputAtCenter'
    // MATLAB Function: '<S112>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S115>:1'
    // '<S115>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
    // HOLDOUTPUTATCENTER_FUNCTION holds the output constant at last input if the 
    // trigger value is within user defined delta from the center
    // 'holdOutputAtCenter_function:5' if isempty(last_input)
    // 'holdOutputAtCenter_function:9' if(trigger <= (params.center + params.posDeltaFromCenter) && ... 
    // 'holdOutputAtCenter_function:10'         trigger >=(params.center - params.negDeltaFromCenter)) 
    if ((fcsModel_DW.rcOutCmds.vzStick_mps <= 0.05) &&
        (fcsModel_DW.rcOutCmds.vzStick_mps >= -0.05)) {
      // 'holdOutputAtCenter_function:11' atCenter = true;
      rtb_Compare_d = true;
    } else {
      // 'holdOutputAtCenter_function:12' else
      // 'holdOutputAtCenter_function:13' atCenter = false;
      rtb_Compare_d = false;

      // 'holdOutputAtCenter_function:14' last_input = input;
      fcsModel_DW.last_input_c = fcsModel_U.stateEstimate.aglEst_m;
    }

    // End of Outputs for SubSystem: '<S103>/holdOutputAtCenter'

    // Switch: '<S103>/Switch' incorporates:
    //   Constant: '<S103>/Constant5'
    //   RelationalOperator: '<S105>/Compare'

    // 'holdOutputAtCenter_function:17' output = last_input;
    if (rtb_Compare_d) {
      rtb_frcCmd_N = 0.0;
    } else {
      rtb_frcCmd_N = fcsModel_DW.rcOutCmds.vzStick_mps;
    }

    // End of Switch: '<S103>/Switch'

    // Outputs for Atomic SubSystem: '<S103>/holdOutputAtCenter'
    // Gain: '<S103>/Gain' incorporates:
    //   MATLAB Function: '<S112>/holdOutputAtCenter'

    rtb_DifferenceInputs2_a = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S103>/holdOutputAtCenter'

    // Gain: '<S103>/Gain1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_DifferenceInputs2_p = -fcsModel_U.stateEstimate.aglEst_m;

    // Concatenate: '<S103>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S103>/Bus Assignment3' incorporates:
    //   Concatenate: '<S103>/Vector Concatenate'
    //   Constant: '<S106>/Constant'
    //   Constant: '<S109>/Constant'
    //   Gain: '<S103>/Gain'
    //   Gain: '<S103>/Gain1'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S103>/Logical Operator'
    //   Logic: '<S103>/Logical Operator1'
    //   MATLAB Function: '<S112>/holdOutputAtCenter'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S106>/Compare'
    //   RelationalOperator: '<S109>/Compare'

    rtb_VectorConcatenate[2].feedForwardCmd = rtb_frcCmd_N;

    // Outputs for Atomic SubSystem: '<S103>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].cmd = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S103>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.aglEst_m;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && (flightMode != enumFlightMode::
      POS_CONTROL)));

    // Outputs for Iterator SubSystem: '<S99>/NED Position Control' incorporates:
    //   ForEach: '<S104>/For Each'

    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // Outputs for Atomic SubSystem: '<S104>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S104>/ctrlInputs' incorporates:
      //   BusAssignment: '<S103>/Bus Assignment'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S104>/Unit Delay'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_frcCmd_N, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S104>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S104>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_DifferenceInputs2_p, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S104>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S104>/pidWithDebug'
      fcsModel_pidWithDebug_j(rtb_VectorConcatenate[ForEach_itr_i].
        feedForwardCmd, rtb_frcCmd_N, rtb_DifferenceInputs2_p,
        rtb_VectorConcatenate[ForEach_itr_i].integratorReset, 0.0,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.ctrlParamsArray[
        ForEach_itr_i], fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE,
        &rtb_DifferenceInputs2_a, &rtb_BusCreator_em, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].pidWithDebug);

      // End of Outputs for SubSystem: '<S104>/pidWithDebug'

      // Update for UnitDelay: '<S104>/Unit Delay'
      fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE =
        rtb_DifferenceInputs2_a;

      // ForEachSliceAssignment generated from: '<S104>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDeb_m[ForEach_itr_i] = rtb_BusCreator_em;

      // ForEachSliceAssignment generated from: '<S104>/neVelCmd_mps'
      rtb_ImpAsg_InsertedFor_neVelCmd[ForEach_itr_i] = rtb_DifferenceInputs2_a;

      // ForEachSliceAssignment generated from: '<S104>/meas'
      rtb_ImpAsg_InsertedFor_meas_at_[ForEach_itr_i] = rtb_DifferenceInputs2_p;

      // ForEachSliceAssignment generated from: '<S104>/cmd'
      rtb_ImpAsg_InsertedFor_cmd_at_i[ForEach_itr_i] = rtb_frcCmd_N;
    }

    // End of Outputs for SubSystem: '<S99>/NED Position Control'

    // RelationalOperator: '<S102>/Compare' incorporates:
    //   Constant: '<S102>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_mi = (flightMode != enumFlightMode::POS_CONTROL);

    // Logic: '<S98>/Logical Operator2'
    rtb_Compare_d = (resetIntegrator || rtb_Compare_mi);

    // Concatenate: '<S98>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S98>/Bus Assignment' incorporates:
    //   Concatenate: '<S98>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[0];
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedVel_mps[0];
    rtb_VectorConcatenate[0].integratorReset = rtb_Compare_d;

    // Concatenate: '<S98>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S98>/Bus Assignment1' incorporates:
    //   Concatenate: '<S98>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[1];
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedVel_mps[1];
    rtb_VectorConcatenate[1].integratorReset = rtb_Compare_d;

    // Gain: '<S98>/Gain' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_DifferenceInputs2_p = -fcsModel_U.stateEstimate.climbRateEst_mps;

    // Concatenate: '<S98>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S98>/Bus Assignment2' incorporates:
    //   Concatenate: '<S98>/Vector Concatenate'
    //   Constant: '<S101>/Constant'
    //   Gain: '<S98>/Gain'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S98>/Logical Operator'
    //   Logic: '<S98>/Logical Operator1'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S101>/Compare'

    rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.climbRateEst_mps;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && rtb_Compare_mi));

    // SignalConversion generated from: '<S100>/For Each Subsystem'
    rtb_VectorConcatenate1[0] = 0.0;
    rtb_VectorConcatenate1[1] = 0.0;
    rtb_VectorConcatenate1[2] = 0.0;

    // Outputs for Iterator SubSystem: '<S100>/For Each Subsystem' incorporates:
    //   ForEach: '<S157>/For Each'

    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // Outputs for Atomic SubSystem: '<S157>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S157>/ctrlInputs' incorporates:
      //   BusAssignment: '<S98>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr], &rtb_DifferenceInputs2_p, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S157>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S157>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr], &rtb_DifferenceInputs2_a, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S157>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S157>/Signal Conditioning Block2'
      // ForEachSliceSelector generated from: '<S157>/nedAccel_mps2' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   Inport: '<Root>/stateEstimate'

      fcsM_SignalConditioningBlock1_f
        (fcsModel_U.stateEstimate.nedAccel_mps2[ForEach_itr],
         &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelSignalConditioningParamsArray
         [ForEach_itr], &rtb_frcCmd_N, 0.02,
         &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock2);

      // End of Outputs for SubSystem: '<S157>/Signal Conditioning Block2'

      // Product: '<S157>/Product' incorporates:
      //   ForEachSliceSelector generated from: '<S157>/accelFbGain'
      //   Inport: '<Root>/ctrlParams'

      rlim =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelFbGainsArray
        [ForEach_itr] * rtb_frcCmd_N;

      // Product: '<S157>/Product1' incorporates:
      //   BusAssignment: '<S98>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S157>/ctrlInputs'
      //   ForEachSliceSelector generated from: '<S157>/ffGain'
      //   Inport: '<Root>/ctrlParams'

      rtb_frcCmd_N =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ffGainsArray[ForEach_itr]
        * rtb_VectorConcatenate[ForEach_itr].cmd;

      // Outputs for Atomic SubSystem: '<S157>/pidWithDebug'
      // Sum: '<S157>/Sum' incorporates:
      //   BusAssignment: '<S98>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S157>/ctrlInputs'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S157>/Unit Delay'

      fcsModel_pidWithDebug_j(rtb_frcCmd_N - rlim, rtb_DifferenceInputs2_p,
        rtb_DifferenceInputs2_a, rtb_VectorConcatenate[ForEach_itr].
        integratorReset, rtb_VectorConcatenate1[ForEach_itr],
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ctrlParamsArray[
        ForEach_itr], fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE,
        &rtb_frcCmd_N, &rtb_BusCreator_em, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].pidWithDebug);

      // End of Outputs for SubSystem: '<S157>/pidWithDebug'

      // Update for UnitDelay: '<S157>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S157>/velCtrlOut '
      rtb_ImpAsg_InsertedFor_velCtrlO[ForEach_itr] = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S157>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr] = rtb_BusCreator_em;

      // ForEachSliceAssignment generated from: '<S157>/velCtrlFf'
      rtb_ImpAsg_InsertedFor_velCtrlF[ForEach_itr] = rlim;

      // ForEachSliceAssignment generated from: '<S157>/filtMeas'
      rtb_ImpAsg_InsertedFor_filtMeas[ForEach_itr] = rtb_DifferenceInputs2_a;

      // ForEachSliceAssignment generated from: '<S157>/filtCmd'
      rtb_ImpAsg_InsertedFor_filtCmd_[ForEach_itr] = rtb_DifferenceInputs2_p;
    }

    // End of Outputs for SubSystem: '<S100>/For Each Subsystem'

    // Trigonometry: '<S169>/Sin' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_DifferenceInputs2_p = std::sin(fcsModel_U.stateEstimate.attitude_rad[2]);

    // Trigonometry: '<S169>/Sin1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = std::cos(fcsModel_U.stateEstimate.attitude_rad[2]);

    // BusAssignment: '<S156>/Bus Assignment1' incorporates:
    //   Constant: '<S169>/Constant4'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S169>/Product1'
    //   Product: '<S169>/Product2'
    //   Product: '<S169>/Product6'
    //   Sum: '<S169>/Sum'

    std::memset(&rtb_BusAssignment1_a, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment1_a.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] *
      rtb_DifferenceInputs2_p - rtb_frcCmd_N * rtb_ImpAsg_InsertedFor_velCtrlO[1])
      * -0.10197838058331635;
    rtb_BusAssignment1_a.meas = fcsModel_U.stateEstimate.attitude_rad[0];
    rtb_BusAssignment1_a.integratorReset = resetIntegrator;

    // BusAssignment: '<S156>/Bus Assignment2' incorporates:
    //   Constant: '<S169>/Constant4'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S169>/Product3'
    //   Product: '<S169>/Product4'
    //   Product: '<S169>/Product7'
    //   Sum: '<S169>/Sum1'

    std::memset(&rtb_BusAssignment2, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment2.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] * rtb_frcCmd_N
      + rtb_DifferenceInputs2_p * rtb_ImpAsg_InsertedFor_velCtrlO[1]) *
      -0.10197838058331635;
    rtb_BusAssignment2.meas = fcsModel_U.stateEstimate.attitude_rad[1];
    rtb_BusAssignment2.integratorReset = resetIntegrator;

    // Logic: '<S156>/Logical Operator' incorporates:
    //   Constant: '<S159>/Constant'
    //   Constant: '<S160>/Constant'
    //   Logic: '<S156>/Logical Operator1'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S159>/Compare'
    //   RelationalOperator: '<S160>/Compare'

    rtb_Compare_d = (resetIntegrator || ((flightMode != enumFlightMode::
      ALT_CONTROL) && (flightMode != enumFlightMode::POS_CONTROL)));

    // Outputs for Atomic SubSystem: '<S156>/holdOutputAtCenter'
    // MATLAB Function: '<S162>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S167>:1'
    // '<S167>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
    // HOLDOUTPUTATCENTER_FUNCTION holds the output constant at last input if the 
    // trigger value is within user defined delta from the center
    // 'holdOutputAtCenter_function:5' if isempty(last_input)
    // 'holdOutputAtCenter_function:9' if(trigger <= (params.center + params.posDeltaFromCenter) && ... 
    // 'holdOutputAtCenter_function:10'         trigger >=(params.center - params.negDeltaFromCenter)) 
    if ((fcsModel_DW.rcOutCmds.yawStick <= 0.017453292519943295) &&
        (fcsModel_DW.rcOutCmds.yawStick >= -0.017453292519943295)) {
      // 'holdOutputAtCenter_function:11' atCenter = true;
      rtb_Compare_mi = true;
    } else {
      // 'holdOutputAtCenter_function:12' else
      // 'holdOutputAtCenter_function:13' atCenter = false;
      rtb_Compare_mi = false;

      // 'holdOutputAtCenter_function:14' last_input = input;
      fcsModel_DW.last_input = fcsModel_U.stateEstimate.attitude_rad[2];
    }

    // End of Outputs for SubSystem: '<S156>/holdOutputAtCenter'

    // BusAssignment: '<S156>/Bus Assignment4'
    // 'holdOutputAtCenter_function:17' output = last_input;
    std::memset(&rtb_BusAssignment4, 0, sizeof(busCtrlInputs));

    // Switch: '<S156>/Switch' incorporates:
    //   RelationalOperator: '<S158>/Compare'

    if (rtb_Compare_mi) {
      // BusAssignment: '<S156>/Bus Assignment4' incorporates:
      //   Constant: '<S156>/Constant5'

      rtb_BusAssignment4.feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S156>/Bus Assignment4'
      rtb_BusAssignment4.feedForwardCmd = fcsModel_DW.rcOutCmds.yawStick;
    }

    // End of Switch: '<S156>/Switch'

    // Outputs for Atomic SubSystem: '<S156>/holdOutputAtCenter'
    // BusAssignment: '<S156>/Bus Assignment4' incorporates:
    //   Inport: '<Root>/stateEstimate'
    //   MATLAB Function: '<S162>/holdOutputAtCenter'

    rtb_BusAssignment4.cmd = fcsModel_DW.last_input;

    // End of Outputs for SubSystem: '<S156>/holdOutputAtCenter'
    rtb_BusAssignment4.meas = fcsModel_U.stateEstimate.attitude_rad[2];
    rtb_BusAssignment4.integratorReset = rtb_Compare_d;

    // Trigonometry: '<S156>/Sin2' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);

    // Product: '<S156>/Product5' incorporates:
    //   Inport: '<Root>/stateEstimate'
    //   Trigonometry: '<S156>/Sin3'

    rlim = rtb_frcCmd_N * std::cos(fcsModel_U.stateEstimate.attitude_rad[1]);

    // BusAssignment: '<S156>/Bus Assignment3' incorporates:
    //   Inport: '<Root>/stateEstimate'

    plim = fcsModel_U.stateEstimate.nedAccel_mps2[2];

    // Outputs for Atomic SubSystem: '<S171>/Signal Conditioning Block1'
    // MATLAB Function: '<S190>/Compute Natural Frequency' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMode_ComputeNaturalFrequency
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.filterBandwidth_radps,
       fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
       &rtb_frcCmd_N);

    // MATLAB Function: '<S190>/Compute Numerator And Denominator' incorporates:
    //   Inport: '<Root>/ctrlParams'

    ComputeNumeratorAndDenominator(rtb_frcCmd_N,
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
      &rtb_VectorConcatenate1[0], &rtb_MatrixMultiply[0],
      &rtb_ImpAsg_InsertedFor_neVelCmd[0], 0.02);

    // DiscreteTransferFcn: '<S190>/Discrete Transfer Fcn' incorporates:
    //   BusAssignment: '<S156>/Bus Assignment3'
    //   Inport: '<Root>/stateEstimate'

    rtb_DifferenceInputs2_a = (fcsModel_U.stateEstimate.nedAccel_mps2[2] -
      fcsModel_DW.DiscreteTransferFcn_states_n[0] *
      rtb_ImpAsg_InsertedFor_neVelCmd[1]) -
      fcsModel_DW.DiscreteTransferFcn_states_n[1] *
      rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_DiscreteTransferFcn = (rtb_VectorConcatenate1[0] *
      rtb_DifferenceInputs2_a + fcsModel_DW.DiscreteTransferFcn_states_n[0] *
      rtb_VectorConcatenate1[1]) + fcsModel_DW.DiscreteTransferFcn_states_n[1] *
      rtb_VectorConcatenate1[2];

    // Switch: '<S196>/Switch2' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   RelationalOperator: '<S196>/LowerRelop1'
    //   RelationalOperator: '<S196>/UpperRelop'
    //   Switch: '<S196>/Switch'

    if (rtb_DiscreteTransferFcn >
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [1]) {
      rtb_DiscreteTransferFcn =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [1];
    } else if (rtb_DiscreteTransferFcn <
               fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
               [0]) {
      // Switch: '<S196>/Switch'
      rtb_DiscreteTransferFcn =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [0];
    }

    // End of Switch: '<S196>/Switch2'

    // Sum: '<S193>/Difference Inputs1' incorporates:
    //   UnitDelay: '<S193>/Delay Input2'
    //
    //  Block description for '<S193>/Difference Inputs1':
    //
    //   Add in CPU
    //
    //  Block description for '<S193>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_frcCmd_N = rtb_DiscreteTransferFcn - fcsModel_DW.DelayInput2_DSTATE;

    // Product: '<S193>/delta rise limit' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   SampleTimeMath: '<S193>/sample time'
    //
    //  About '<S193>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_DifferenceInputs2_p =
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputAccelLimits
      [1] * 0.02;

    // Switch: '<S203>/Switch2' incorporates:
    //   RelationalOperator: '<S203>/LowerRelop1'

    if (rtb_frcCmd_N <= rtb_DifferenceInputs2_p) {
      // Product: '<S193>/delta fall limit' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   SampleTimeMath: '<S193>/sample time'
      //
      //  About '<S193>/sample time':
      //   y = K where K = ( w * Ts )

      rtb_DifferenceInputs2_p =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputAccelLimits
        [0] * 0.02;

      // Switch: '<S203>/Switch' incorporates:
      //   RelationalOperator: '<S203>/UpperRelop'

      if (rtb_frcCmd_N >= rtb_DifferenceInputs2_p) {
        rtb_DifferenceInputs2_p = rtb_frcCmd_N;
      }

      // End of Switch: '<S203>/Switch'
    }

    // End of Switch: '<S203>/Switch2'

    // Sum: '<S193>/Difference Inputs2' incorporates:
    //   UnitDelay: '<S193>/Delay Input2'
    //
    //  Block description for '<S193>/Difference Inputs2':
    //
    //   Add in CPU
    //
    //  Block description for '<S193>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DiscreteTransferFcn = rtb_DifferenceInputs2_p +
      fcsModel_DW.DelayInput2_DSTATE;

    // MATLAB Function: '<S191>/Compute Natural Frequency' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMode_ComputeNaturalFrequency
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.filterBandwidth_radps,
       fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
       &rtb_frcCmd_N);

    // MATLAB Function: '<S191>/Compute Filter Numerator And Denominator' incorporates:
    //   Inport: '<Root>/ctrlParams'

    ComputeFilterNumeratorAndDenomi(rtb_frcCmd_N,
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
      &rtb_VectorConcatenate1[0], &rtb_MatrixMultiply[0], 0.02);

    // DiscreteTransferFcn: '<S191>/Discrete Transfer Fcn' incorporates:
    //   BusAssignment: '<S156>/Bus Assignment3'
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = (fcsModel_U.stateEstimate.nedAccel_mps2[2] -
                    fcsModel_DW.DiscreteTransferFcn_states_nh[0] *
                    rtb_MatrixMultiply[1]) -
      fcsModel_DW.DiscreteTransferFcn_states_nh[1] * rtb_MatrixMultiply[2];
    rtb_DiscreteTransferFcn_k = (rtb_VectorConcatenate1[0] * rtb_frcCmd_N +
      fcsModel_DW.DiscreteTransferFcn_states_nh[0] * rtb_VectorConcatenate1[1])
      + fcsModel_DW.DiscreteTransferFcn_states_nh[1] * rtb_VectorConcatenate1[2];

    // Switch: '<S195>/Switch2' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   RelationalOperator: '<S195>/LowerRelop1'
    //   RelationalOperator: '<S195>/UpperRelop'
    //   Switch: '<S195>/Switch'

    if (rtb_DiscreteTransferFcn_k >
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [1]) {
      rtb_DiscreteTransferFcn_k =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [1];
    } else if (rtb_DiscreteTransferFcn_k <
               fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
               [0]) {
      // Switch: '<S195>/Switch'
      rtb_DiscreteTransferFcn_k =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [0];
    }

    // End of Switch: '<S195>/Switch2'

    // Sum: '<S192>/Difference Inputs1' incorporates:
    //   UnitDelay: '<S192>/Delay Input2'
    //
    //  Block description for '<S192>/Difference Inputs1':
    //
    //   Add in CPU
    //
    //  Block description for '<S192>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DiscreteTransferFcn_k -= fcsModel_DW.DelayInput2_DSTATE_e;

    // Product: '<S192>/delta rise limit' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   SampleTimeMath: '<S192>/sample time'
    //
    //  About '<S192>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_DifferenceInputs2_p =
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
      [1] * 0.02;

    // Switch: '<S202>/Switch2' incorporates:
    //   RelationalOperator: '<S202>/LowerRelop1'

    if (rtb_DiscreteTransferFcn_k <= rtb_DifferenceInputs2_p) {
      // Product: '<S192>/delta fall limit' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   SampleTimeMath: '<S192>/sample time'
      //
      //  About '<S192>/sample time':
      //   y = K where K = ( w * Ts )

      rtb_DifferenceInputs2_p =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [0] * 0.02;

      // Switch: '<S202>/Switch' incorporates:
      //   RelationalOperator: '<S202>/UpperRelop'

      if (rtb_DiscreteTransferFcn_k >= rtb_DifferenceInputs2_p) {
        rtb_DifferenceInputs2_p = rtb_DiscreteTransferFcn_k;
      }

      // End of Switch: '<S202>/Switch'
    }

    // End of Switch: '<S202>/Switch2'

    // Sum: '<S192>/Difference Inputs2' incorporates:
    //   UnitDelay: '<S192>/Delay Input2'
    //
    //  Block description for '<S192>/Difference Inputs2':
    //
    //   Add in CPU
    //
    //  Block description for '<S192>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DiscreteTransferFcn_k = rtb_DifferenceInputs2_p +
      fcsModel_DW.DelayInput2_DSTATE_e;

    // Update for DiscreteTransferFcn: '<S190>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states_n[1] =
      fcsModel_DW.DiscreteTransferFcn_states_n[0];
    fcsModel_DW.DiscreteTransferFcn_states_n[0] = rtb_DifferenceInputs2_a;

    // Update for UnitDelay: '<S193>/Delay Input2'
    //
    //  Block description for '<S193>/Delay Input2':
    //
    //   Store in Global RAM

    fcsModel_DW.DelayInput2_DSTATE = rtb_DiscreteTransferFcn;

    // Update for DiscreteTransferFcn: '<S191>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states_nh[1] =
      fcsModel_DW.DiscreteTransferFcn_states_nh[0];
    fcsModel_DW.DiscreteTransferFcn_states_nh[0] = rtb_frcCmd_N;

    // Update for UnitDelay: '<S192>/Delay Input2'
    //
    //  Block description for '<S192>/Delay Input2':
    //
    //   Store in Global RAM

    fcsModel_DW.DelayInput2_DSTATE_e = rtb_DiscreteTransferFcn_k;

    // End of Outputs for SubSystem: '<S171>/Signal Conditioning Block1'

    // Product: '<S171>/Product' incorporates:
    //   Gain: '<S171>/Gain'
    //   Inport: '<Root>/ctrlParams'

    rtb_DifferenceInputs2_p = -rtb_DiscreteTransferFcn *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ctrlParams.Kd;

    // Outputs for Atomic SubSystem: '<S171>/Signal Conditioning Block'
    // BusAssignment: '<S156>/Bus Assignment3' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsM_SignalConditioningBlock1_f(rtb_ImpAsg_InsertedFor_velCtrlO[2],
      &fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.cmdSignalConditioningParams,
      &rtb_DifferenceInputs2_a, 0.02, &fcsModel_DW.SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S171>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S171>/pidWithDebug'
    // Sum: '<S171>/Sum' incorporates:
    //   BusAssignment: '<S156>/Bus Assignment3'
    //   Inport: '<Root>/ctrlParams'
    //   Product: '<S171>/Product1'
    //   UnitDelay: '<S171>/Unit Delay'

    fcsModel_pidWithDebug_j(rtb_ImpAsg_InsertedFor_velCtrlO[2] *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ffGain +
      rtb_DifferenceInputs2_p, rtb_DifferenceInputs2_a,
      rtb_DiscreteTransferFcn_k, rtb_Compare_d, 0.0,
      &fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ctrlParams,
      fcsModel_DW.UnitDelay_DSTATE, &rtb_DifferenceInputs2_p, &rtb_BusCreator_em,
      0.02, &fcsModel_DW.pidWithDebug);

    // End of Outputs for SubSystem: '<S171>/pidWithDebug'

    // Product: '<S156>/Divide1' incorporates:
    //   Constant: '<S168>/g'
    //   Inport: '<Root>/ctrlParams'
    //   Product: '<S168>/Divide3'
    //   Sum: '<S156>/Sum2'

    rtb_frcCmd_N =
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.baseMass_kg *
       -9.806 + rtb_DifferenceInputs2_p) * (1.0 / rlim);

    // RelationalOperator: '<S95>/Compare' incorporates:
    //   Constant: '<S95>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_d = (flightMode == enumFlightMode::POS_CONTROL);

    // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
    //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
    //   Constant: '<S3>/Constant'
    //   Inport: '<Root>/stateEstimate'

    std::memcpy(&rtb_VectorConcatenate[0],
                &fcsModel_ConstP.pooled3.attCtrlInputs.ctrlInputsArray[0], 3U *
                sizeof(busCtrlInputs));

    // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S97>:1' 
    // '<S97>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
    // '<S97>:1:3' outDebug = throttleCmd_N;
    rtb_DiscreteTransferFcn = fcsModel_DW.rcOutCmds.throttleStick;

    //  This is a stop gap setup where we are only assuming that rate control
    //  is active and therefore not setting up attCtrlInputs for Euler angle
    //  control
    // '<S97>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
    rtb_VectorConcatenate[0].cmd = fcsModel_DW.rcOutCmds.rollStick;

    // '<S97>:1:8' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.attitude_rad[0];

    // '<S97>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
    rtb_VectorConcatenate[1].cmd = fcsModel_DW.rcOutCmds.pitchStick;

    // '<S97>:1:10' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.attitude_rad[1];

    // '<S97>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
    rtb_VectorConcatenate[2].cmd = fcsModel_DW.rcOutCmds.yawStick;

    // '<S97>:1:12' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
    rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.attitude_rad[2];

    // '<S97>:1:14' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[0].integratorReset = resetIntegrator;

    // '<S97>:1:15' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[1].integratorReset = resetIntegrator;

    // '<S97>:1:16' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = true; 
    rtb_VectorConcatenate[2].integratorReset = true;

    // RelationalOperator: '<S94>/Compare' incorporates:
    //   Constant: '<S94>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_mi = (flightMode != enumFlightMode::ALT_CONTROL);

    // Switch: '<S3>/Switch2' incorporates:
    //   Switch: '<S3>/Switch'

    if (rtb_Compare_d) {
      // Switch: '<S3>/Switch2' incorporates:
      //   BusAssignment: '<S156>/Bus Assignment'
      //   Concatenate: '<S156>/Vector Concatenate'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0] =
        rtb_BusAssignment1_a;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[1] = rtb_BusAssignment2;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment4;
    } else if (rtb_Compare_mi) {
      // Switch: '<S3>/Switch2' incorporates:
      //   MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus'
      //   Switch: '<S3>/Switch'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N =
        fcsModel_DW.rcOutCmds.throttleStick;
      std::memcpy(&fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0],
                  &rtb_VectorConcatenate[0], 3U * sizeof(busCtrlInputs));
    } else {
      // Switch: '<S3>/Switch2' incorporates:
      //   BusAssignment: '<S156>/Bus Assignment'
      //   Concatenate: '<S156>/Vector Concatenate'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Switch: '<S3>/Switch'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      std::memcpy(&fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0],
                  &rtb_VectorConcatenate[0], sizeof(busCtrlInputs) << 1U);
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment4;
    }

    // End of Switch: '<S3>/Switch2'
  }

  // Outputs for Iterator SubSystem: '<S10>/Attitude Control' incorporates:
  //   ForEach: '<S52>/For Each'

  for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
    // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S52>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].cmd,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.cmdSignalConditioningParamsArray
       [ForEach_itr_l], &rlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S52>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].meas,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.measSignalConditioningParamsArray
       [ForEach_itr_l], &plim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S52>/Signal Conditioning Block1'

    // MATLAB Function: '<S52>/pickAttitudeCmdAndMeas' incorporates:
    //   Constant: '<S10>/Constant'
    //   ForEachSliceSelector generated from: '<S52>/index'

    vxCmd_unitRange = rlim;
    yCmd = plim;

    //  Passes cmd and meas as it is for roll and pitch channel
    //  but for yaw channel computes shortest angular distance between cmd Yaw
    //  and meas Yaw and overwrites Yaw cmd with that error and sets the meas Yaw to 
    //  zero for PID block
    // MATLAB Function 'Attitude Controller/Attitude Control/pickAttitudeCmdAndMeas': '<S57>:1' 
    // '<S57>:1:6' if index == cast(3, 'uint8')
    if (fcsModel_ConstP.Constant_Value_e[ForEach_itr_l] == 3) {
      // '<S57>:1:7' diff = mod(( cmd - meas + pi ), 2*pi) - pi;
      vxCmd_unitRange = (rlim - plim) + 3.1415926535897931;
      if (vxCmd_unitRange == 0.0) {
        ylim = 0.0;
      } else {
        ylim = std::fmod(vxCmd_unitRange, 6.2831853071795862);
        resetIntegrator = (ylim == 0.0);
        if (!resetIntegrator) {
          yCmd = std::abs(vxCmd_unitRange / 6.2831853071795862);
          resetIntegrator = (std::abs(yCmd - std::floor(yCmd + 0.5)) <=
                             2.2204460492503131E-16 * yCmd);
        }

        if (resetIntegrator) {
          ylim = 0.0;
        } else if (vxCmd_unitRange < 0.0) {
          ylim += 6.2831853071795862;
        }
      }

      vxCmd_unitRange = ylim - 3.1415926535897931;

      // '<S57>:1:8' if diff < -pi
      if (ylim - 3.1415926535897931 < -3.1415926535897931) {
        // '<S57>:1:9' diff = diff + 2*pi;
        vxCmd_unitRange = (ylim - 3.1415926535897931) + 6.2831853071795862;
      }

      // '<S57>:1:12' cmd = diff;
      // '<S57>:1:13' meas = 0;
      yCmd = 0.0;
    }

    // Outputs for Atomic SubSystem: '<S52>/pidWithDebug'
    // ForEachSliceSelector generated from: '<S52>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   MATLAB Function: '<S52>/pickAttitudeCmdAndMeas'
    //   UnitDelay: '<S52>/Unit Delay'

    fcsModel_pidWithDebug
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       feedForwardCmd, vxCmd_unitRange, yCmd,
       fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       integratorReset, 0.0,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.ctrlParamsArray[ForEach_itr_l],
       fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE, &ylim,
       &rtb_BusCreator_f, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       pidWithDebug);

    // End of Outputs for SubSystem: '<S52>/pidWithDebug'

    // Update for UnitDelay: '<S52>/Unit Delay'
    fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE = ylim;

    // ForEachSliceAssignment generated from: '<S52>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[ForEach_itr_l] =
      rtb_BusCreator_f;

    // ForEachSliceAssignment generated from: '<S52>/angRateCmd '
    rtb_ImpAsg_InsertedFor_angRateC[ForEach_itr_l] = ylim;

    // ForEachSliceAssignment generated from: '<S52>/measFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[ForEach_itr_l] =
      plim;

    // ForEachSliceAssignment generated from: '<S52>/cmdFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[ForEach_itr_l] =
      rlim;
  }

  // End of Outputs for SubSystem: '<S10>/Attitude Control'

  // Switch: '<S9>/Switch'
  rtb_ImpAsg_InsertedFor_velCtrlO[0] = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_ImpAsg_InsertedFor_velCtrlO[1] = rtb_ImpAsg_InsertedFor_angRateC[1];

  // Switch: '<S10>/Switch' incorporates:
  //   Constant: '<S53>/Constant'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   RelationalOperator: '<S53>/Compare'
  //   Switch: '<S9>/Switch'

  if (flightMode == enumFlightMode::STABILIZE) {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] =
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2].cmd;
  } else {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] = rtb_ImpAsg_InsertedFor_angRateC[2];
  }

  // End of Switch: '<S10>/Switch'

  // MATLAB Function: '<S9>/EulerRates2BodyRates' incorporates:
  //   Inport: '<Root>/stateEstimate'
  //   Switch: '<S9>/Switch'

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
  rlim = fcsModel_U.stateEstimate.attitude_rad[1];

  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  plim = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(plim - 1.5707963267948966) <= 0.004245395477824045) || (std::abs
       (plim - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rlim = -1.0;
      } else {
        rlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rlim *= plim - 0.004245395477824045;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rlim = -1.0;
      } else {
        rlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rlim *= plim - 0.004245395477824045;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        rlim = -1.0;
      } else {
        rlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      rlim *= plim + 0.004245395477824045;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  plim = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  ylim = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  vxCmd_unitRange = std::cos(rlim);
  rtb_Transpose[0] = 1.0;
  rtb_Transpose[3] = 0.0;
  rtb_Transpose[6] = -std::sin(rlim);
  rtb_Transpose[1] = 0.0;
  rtb_Transpose[4] = ylim;
  rtb_Transpose[7] = plim * vxCmd_unitRange;
  rtb_Transpose[2] = 0.0;
  rtb_Transpose[5] = -plim;
  rtb_Transpose[8] = ylim * vxCmd_unitRange;

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
  for (tCmd = 0; tCmd < 3; tCmd++) {
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    rlim = rtb_Transpose[tCmd];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (rlim <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      rlim = 0.0;
    }

    rtb_Transpose[tCmd] = rlim;
    plim = rlim * rtb_ImpAsg_InsertedFor_velCtrlO[0];
    rlim = rtb_Transpose[tCmd + 3];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(rlim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      rlim = 0.0;
    }

    rtb_Transpose[tCmd + 3] = rlim;
    plim += rlim * rtb_ImpAsg_InsertedFor_velCtrlO[1];
    rlim = rtb_Transpose[tCmd + 6];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(rlim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      rlim = 0.0;
    }

    rtb_Transpose[tCmd + 6] = rlim;
    rtb_ImpAsg_InsertedFor_angRateC[tCmd] = rlim *
      rtb_ImpAsg_InsertedFor_velCtrlO[2] + plim;
  }

  // End of MATLAB Function: '<S9>/EulerRates2BodyRates'

  // BusCreator: '<S9>/Bus Creator' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_VectorConcatenate[0].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S9>/Bus Creator3' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_VectorConcatenate[1].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S9>/Bus Creator4' incorporates:
  //   Concatenate: '<S9>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[2].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_angRateC[2];
  rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_VectorConcatenate[2].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[2].trackingCtrlCmd = 0.0;

  // Outputs for Atomic SubSystem: '<S2>/Angular Rate Controller'
  // Outputs for Iterator SubSystem: '<S8>/For Each Subsystem' incorporates:
  //   ForEach: '<S11>/For Each'

  for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
    // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S11>/ctrlInputs' incorporates:
    //   BusCreator: '<S9>/Bus Creator1'
    //   Concatenate: '<S9>/Vector Concatenate'
    //   Inport: '<Root>/ctrlParams'
    //   UnitDelay: '<S11>/Unit Delay'

    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].cmd,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray
      [ForEach_itr_p], &ylim, 0.004, &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
      SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S11>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].meas,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.measSignalConditioningParamsArray
      [ForEach_itr_p], &plim, 0.004, &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
      SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S11>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S11>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, ylim, plim, rtb_VectorConcatenate[ForEach_itr_p].
                          integratorReset, 0.0,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.ctrlParamsArray
                          [ForEach_itr_p],
                          fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
                          UnitDelay_DSTATE, &rlim, &rtb_BusCreator_f, 0.004,
                          &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].pidWithDebug);

    // End of Outputs for SubSystem: '<S11>/pidWithDebug'

    // Update for UnitDelay: '<S11>/Unit Delay'
    fcsModel_DW.CoreSubsys_a[ForEach_itr_p].UnitDelay_DSTATE = rlim;

    // ForEachSliceAssignment generated from: '<S11>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[ForEach_itr_p]
      = rtb_BusCreator_f;

    // ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'
    rtb_ImpAsg_InsertedFor_angAccel[ForEach_itr_p] = rlim;

    // ForEachSliceAssignment generated from: '<S11>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[ForEach_itr_p] =
      plim;

    // ForEachSliceAssignment generated from: '<S11>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[ForEach_itr_p] =
      ylim;
  }

  // End of Outputs for SubSystem: '<S8>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S11>/angAccelCmd_radps2'

  for (tCmd = 0; tCmd < 3; tCmd++) {
    rtb_ImpAsg_InsertedFor_angRateC[tCmd] = 0.0;
    rtb_ImpAsg_InsertedFor_angRateC[tCmd] +=
      fcsModel_ConstP.Constant_Value_n[tCmd] * rtb_ImpAsg_InsertedFor_angAccel[0];
    rtb_ImpAsg_InsertedFor_angRateC[tCmd] +=
      fcsModel_ConstP.Constant_Value_n[tCmd + 3] *
      rtb_ImpAsg_InsertedFor_angAccel[1];
    rtb_ImpAsg_InsertedFor_angRateC[tCmd] +=
      fcsModel_ConstP.Constant_Value_n[tCmd + 6] *
      rtb_ImpAsg_InsertedFor_angAccel[2];
  }

  // End of Product: '<S2>/Matrix Multiply'

  // Gain: '<S1>/Gain' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'

  plim = rtb_ImpAsg_InsertedFor_angRateC[0];
  ylim = rtb_ImpAsg_InsertedFor_angRateC[1];
  vxCmd_unitRange = 0.7 * rtb_ImpAsg_InsertedFor_angRateC[2];

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);
  for (tCmd = 0; tCmd < 4; tCmd++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   BusCreator: '<S2>/Bus Creator1'
    //   Constant: '<S1>/Constant'
    //   Gain: '<S1>/Gain'

    yCmd = ((fcsModel_ConstP.Constant_Value_c[tCmd + 4] * plim +
             fcsModel_ConstP.Constant_Value_c[tCmd] *
             fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N) +
            fcsModel_ConstP.Constant_Value_c[tCmd + 8] * ylim) +
      fcsModel_ConstP.Constant_Value_c[tCmd + 12] * vxCmd_unitRange;

    // Saturate: '<S1>/Saturation'
    if (yCmd > 839601.76328711538) {
      yCmd = 839601.76328711538;
    } else if (yCmd < 0.0) {
      yCmd = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    rlim = 9.5492965855137211 * std::sqrt(yCmd) - 0.11372544828835565 *
      fcsModel_DW.DiscreteTransferFcn_states[tCmd];
    yCmd = 0.55686272414417781 * rlim + 0.55686272414417781 *
      fcsModel_DW.DiscreteTransferFcn_states[tCmd];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      yCmd = -1.0;
    }

    // End of Switch: '<S1>/Switch'

    // Outport: '<Root>/actuatorsCmds'
    fcsModel_Y.actuatorsCmds[tCmd] = yCmd;

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    DiscreteTransferFcn_tmp[tCmd] = rlim;

    // Product: '<S1>/Matrix Multiply' incorporates:
    //   Constant: '<S1>/Constant'

    rtb_DiscreteTransferFcn_e[tCmd] = yCmd;
  }

  // Outputs for Iterator SubSystem: '<S1>/For Each Subsystem' incorporates:
  //   ForEach: '<S7>/For Each'

  for (ForEach_itr_g = 0; ForEach_itr_g < 4; ForEach_itr_g++) {
    uint32_T rtb_Prelookup_o1;

    // ForEachSliceSelector generated from: '<S7>/propellerSpdCmds_rpm'
    rlim = rtb_DiscreteTransferFcn_e[ForEach_itr_g];

    // Saturate: '<S7>/Saturation'
    if (rlim > 9325.0) {
      // PreLookup: '<S7>/Prelookup'
      rlim = 9325.0;
    } else if (rlim < 2250.0) {
      // PreLookup: '<S7>/Prelookup'
      rlim = 2250.0;
    }

    // End of Saturate: '<S7>/Saturation'

    // PreLookup: '<S7>/Prelookup'
    rtb_Prelookup_o1 = plook_bincpag(rlim,
      &fcsModel_ConstP.Prelookup_BreakpointsData[0], 14U, &rlim,
      &fcsModel_DW.CoreSubsys[ForEach_itr_g].Prelookup_DWORK1);

    // ForEachSliceAssignment generated from: '<S7>/mtrPwmCmds' incorporates:
    //   Interpolation_n-D: '<S7>/Interpolation Using Prelookup'

    rtb_ImpAsg_InsertedFor_mtrPwmCm[ForEach_itr_g] = intrp1d_la(rtb_Prelookup_o1,
      rlim, &fcsModel_ConstP.InterpolationUsingPrelookup_Tab[0], 14U);
  }

  // End of Outputs for SubSystem: '<S1>/For Each Subsystem'

  // Outport: '<Root>/actuatorsPwmCmds' incorporates:
  //   ForEachSliceAssignment generated from: '<S7>/mtrPwmCmds'

  fcsModel_Y.actuatorsPwmCmds[0] = rtb_ImpAsg_InsertedFor_mtrPwmCm[0];
  fcsModel_Y.actuatorsPwmCmds[1] = rtb_ImpAsg_InsertedFor_mtrPwmCm[1];
  fcsModel_Y.actuatorsPwmCmds[2] = rtb_ImpAsg_InsertedFor_mtrPwmCm[2];
  fcsModel_Y.actuatorsPwmCmds[3] = rtb_ImpAsg_InsertedFor_mtrPwmCm[3];

  // RateTransition: '<Root>/Rate Transition'
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_Y.fcsDebug.outerLoopCtrlDebug = fcsModel_DW.RateTransition_Buffer0;

    // Switch: '<S3>/Switch3' incorporates:
    //   RateTransition: '<Root>/Rate Transition'
    //   Switch: '<S3>/Switch1'

    if (rtb_Compare_d) {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S156>/Bus Assignment'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    } else if (rtb_Compare_mi) {
      // Switch: '<S3>/Switch1' incorporates:
      //   BusCreator: '<S3>/Bus Creator'

      rtb_BusCreator_b_frcCmd_N = rtb_DiscreteTransferFcn;
    } else {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S156>/Bus Assignment'
      //   Switch: '<S3>/Switch1'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    }

    // End of Switch: '<S3>/Switch3'

    // BusCreator: '<S3>/Bus Creator' incorporates:
    //   BusCreator: '<S100>/Bus Creator'
    //   BusCreator: '<S171>/Bus Creator'
    //   BusCreator: '<S99>/Bus Creator'
    //   ForEachSliceAssignment generated from: '<S104>/cmd'
    //   ForEachSliceAssignment generated from: '<S157>/filtCmd'

    rtb_BusCreator_b_velCtrlDebug_c = rtb_ImpAsg_InsertedFor_filtCmd_[0];
    rtb_BusCreator_b_velCtrlDebug_m = rtb_ImpAsg_InsertedFor_filtMeas[0];
    rtb_BusCreator_b_velCtrlDebug_v = rtb_ImpAsg_InsertedFor_velCtrlF[0];
    rtb_BusCreator_b_velCtrlDebug_p = rtb_ImpAsg_InsertedFor_pidDebug[0];
    rtb_BusCreator_b_posCtrlDebug_c = rtb_ImpAsg_InsertedFor_cmd_at_i[0];
    rtb_BusCreator_b_posCtrlDebug_m = rtb_ImpAsg_InsertedFor_meas_at_[0];
    rtb_BusCreator_b_posCtrlDebug_p = rtb_ImpAsg_InsertedFor_pidDeb_m[0];
    rtb_BusCreator_b_velCtrlDebug_0 = rtb_ImpAsg_InsertedFor_filtCmd_[1];
    rtb_BusCreator_b_velCtrlDebug_1 = rtb_ImpAsg_InsertedFor_filtMeas[1];
    rtb_BusCreator_b_velCtrlDebug_2 = rtb_ImpAsg_InsertedFor_velCtrlF[1];
    rtb_BusCreator_b_velCtrlDebug_3 = rtb_ImpAsg_InsertedFor_pidDebug[1];
    rtb_BusCreator_b_posCtrlDebug_0 = rtb_ImpAsg_InsertedFor_cmd_at_i[1];
    rtb_BusCreator_b_posCtrlDebug_1 = rtb_ImpAsg_InsertedFor_meas_at_[1];
    rtb_BusCreator_b_posCtrlDebug_2 = rtb_ImpAsg_InsertedFor_pidDeb_m[1];
    rtb_BusCreator_b_velCtrlDebug_4 = rtb_ImpAsg_InsertedFor_filtCmd_[2];
    rtb_BusCreator_b_velCtrlDebug_5 = rtb_ImpAsg_InsertedFor_filtMeas[2];
    rtb_BusCreator_b_velCtrlDebug_6 = rtb_ImpAsg_InsertedFor_velCtrlF[2];
    rtb_BusCreator_b_velCtrlDebug_7 = rtb_ImpAsg_InsertedFor_pidDebug[2];
    rtb_BusCreator_b_posCtrlDebug_3 = rtb_ImpAsg_InsertedFor_cmd_at_i[2];
    rtb_BusCreator_b_posCtrlDebug_4 = rtb_ImpAsg_InsertedFor_meas_at_[2];
    rtb_BusCreator_b_posCtrlDebug_5 = rtb_ImpAsg_InsertedFor_pidDeb_m[2];
    rtb_BusCreator_b_zAccelCtrlDebu = rtb_DifferenceInputs2_a;
    rtb_BusCreator_b_zAccelCtrlDe_0 = rtb_DiscreteTransferFcn_k;
    rtb_BusCreator_b_zAccelCtrlDe_1 = rtb_BusCreator_em.output;
    rtb_BusCreator_b_zAccelCtrlDe_2 = rtb_BusCreator_em.proportionalOutput;
    rtb_BusCreator_b_zAccelCtrlDe_3 = rtb_BusCreator_em.integralOutput;
    rtb_BusCreator_b_zAccelCtrlDe_4 = rtb_BusCreator_em.derivativeOutput;
    rtb_BusCreator_b_xyBodyAccelCtr = fcsModel_rtZbusXyBodyAccelCtrIDebug;

    // Update for UnitDelay: '<S171>/Unit Delay'
    fcsModel_DW.UnitDelay_DSTATE = rtb_DifferenceInputs2_p;
  }

  // End of RateTransition: '<Root>/Rate Transition'

  // BusCreator: '<Root>/Bus Creator' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   Outport: '<Root>/fcsDebug'

  fcsModel_Y.fcsDebug.allocDebug.thrustCmd_N =
    fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N;
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
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.velCtrlFf[0] =
      rtb_BusCreator_b_velCtrlDebug_v;
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
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.velCtrlFf[1] =
      rtb_BusCreator_b_velCtrlDebug_2;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[1] =
      rtb_BusCreator_b_velCtrlDebug_3;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.cmd[1] =
      rtb_BusCreator_b_posCtrlDebug_0;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.meas[1] =
      rtb_BusCreator_b_posCtrlDebug_1;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.pidDebug[1] =
      rtb_BusCreator_b_posCtrlDebug_2;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.cmd[2] =
      rtb_BusCreator_b_velCtrlDebug_4;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.meas[2] =
      rtb_BusCreator_b_velCtrlDebug_5;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.velCtrlFf[2] =
      rtb_BusCreator_b_velCtrlDebug_6;
    fcsModel_DW.RateTransition_Buffer0.velCtrlDebug.pidDebug[2] =
      rtb_BusCreator_b_velCtrlDebug_7;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.cmd[2] =
      rtb_BusCreator_b_posCtrlDebug_3;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.meas[2] =
      rtb_BusCreator_b_posCtrlDebug_4;
    fcsModel_DW.RateTransition_Buffer0.posCtrlDebug.pidDebug[2] =
      rtb_BusCreator_b_posCtrlDebug_5;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.cmd =
      rtb_BusCreator_b_zAccelCtrlDebu;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.meas =
      rtb_BusCreator_b_zAccelCtrlDe_0;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.pidDebug.output =
      rtb_BusCreator_b_zAccelCtrlDe_1;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.pidDebug.proportionalOutput
      = rtb_BusCreator_b_zAccelCtrlDe_2;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.pidDebug.integralOutput =
      rtb_BusCreator_b_zAccelCtrlDe_3;
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug.pidDebug.derivativeOutput
      = rtb_BusCreator_b_zAccelCtrlDe_4;
    fcsModel_DW.RateTransition_Buffer0.xyBodyAccelCtrlDebug =
      rtb_BusCreator_b_xyBodyAccelCtr;
  }

  // End of Update for RateTransition: '<Root>/Rate Transition'
  rate_scheduler((&fcsModel_M));
}

// Model initialize function
void fcsModel::initialize()
{
  // Registration code

  // external outputs
  fcsModel_Y.fcsDebug = fcsModel_rtZbusFcsDebug;

  {
    // local scratch DWork variables
    int32_T ForEach_itr;
    int32_T ForEach_itr_i;
    int32_T ForEach_itr_l;
    int32_T ForEach_itr_p;
    std::array<real_T, 3> den;
    std::array<real_T, 3> num;

    // 'interpretRcInputs_function:24' throttle_is_up = false;
    // 'holdOutputAtCenter_function:6' last_input = 0;
    // SystemInitialize for Iterator SubSystem: '<S99>/NED Position Control'
    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // SystemInitialize for Iterator SubSystem: '<S99>/NED Position Control'
      // SystemInitialize for Atomic SubSystem: '<S104>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S104>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S104>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S104>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S104>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S104>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S99>/NED Position Control'
    }

    // End of SystemInitialize for SubSystem: '<S99>/NED Position Control'
    // SystemInitialize for Iterator SubSystem: '<S100>/For Each Subsystem'
    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // SystemInitialize for Iterator SubSystem: '<S100>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S157>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S157>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S157>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S157>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S157>/Signal Conditioning Block2' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock2);

      // End of SystemInitialize for SubSystem: '<S157>/Signal Conditioning Block2' 

      // SystemInitialize for Atomic SubSystem: '<S157>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S157>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S100>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S100>/For Each Subsystem'
    // SystemInitialize for Atomic SubSystem: '<S171>/Signal Conditioning Block1' 
    // SystemInitialize for MATLAB Function: '<S191>/Compute Filter Numerator And Denominator' 
    // 'holdOutputAtCenter_function:6' last_input = 0;
    ComputeFilterNumeratorAndD_Init(&num[0], &den[0]);

    // End of SystemInitialize for SubSystem: '<S171>/Signal Conditioning Block1' 

    // SystemInitialize for Atomic SubSystem: '<S171>/Signal Conditioning Block' 
    SignalConditioningBlock1_c_Init(&fcsModel_DW.SignalConditioningBlock);

    // End of SystemInitialize for SubSystem: '<S171>/Signal Conditioning Block' 

    // SystemInitialize for Atomic SubSystem: '<S171>/pidWithDebug'
    fcsModel_pidWithDebug_m_Init(&fcsModel_DW.pidWithDebug);

    // End of SystemInitialize for SubSystem: '<S171>/pidWithDebug'
    // SystemInitialize for Iterator SubSystem: '<S10>/Attitude Control'
    for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
      // SystemInitialize for Iterator SubSystem: '<S10>/Attitude Control'
      // SystemInitialize for Atomic SubSystem: '<S52>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S52>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S52>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S52>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S52>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S52>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S10>/Attitude Control'
    }

    // End of SystemInitialize for SubSystem: '<S10>/Attitude Control'
    // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
    // SystemInitialize for Iterator SubSystem: '<S8>/For Each Subsystem'
    for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
      // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
      // SystemInitialize for Iterator SubSystem: '<S8>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S11>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S11>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S11>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S11>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S11>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S11>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S8>/For Each Subsystem'
      // End of SystemInitialize for SubSystem: '<S2>/Angular Rate Controller'
    }

    // End of SystemInitialize for SubSystem: '<S8>/For Each Subsystem'
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
