//
// File: fcsModel.cpp
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

  { enumChirpTrigger::OFF,             // chirpTrigger
    enumChirpType::NONE,               // chirpType
    0.0                                // chirpSignal
  },                                   // sysIdDebug
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
  if ((fcsModel_M->Timing.TaskCounters.TID[1]) > 1) {// Sample time: [0.008s, 0.0s] 
    fcsModel_M->Timing.TaskCounters.TID[1] = 0;
  }
}

//
// Output and update for atomic system:
//    '<S24>/Discrete First Order Deriv Filter'
//    '<S68>/Discrete First Order Deriv Filter'
//    '<S130>/Discrete First Order Deriv Filter'
//    '<S185>/Discrete First Order Deriv Filter'
//    '<S224>/Discrete First Order Deriv Filter'
//
void fcsModel::f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
  rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
  rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW)
{
  real_T K;
  real_T normalizer;
  real_T num_tmp;

  // MATLAB Function: '<S55>/Compute Deriv Filter Numerator And Denominator'
  //  Call the main function
  // MATLAB Function 'Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator': '<S58>:1' 
  // '<S58>:1:4' [num, den] = computeFirstOrderDerivFilterNumAndDen_function(filterBandwidth_radps, sampleTime_s); 
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

  // DiscreteTransferFcn: '<S55>/Discrete Transfer Fcn'
  K = rtu_input - localDW->den[1] * localDW->DiscreteTransferFcn_states;
  *rty_filteredInputRate = localDW->num[0] * K + localDW->num[1] *
    localDW->DiscreteTransferFcn_states;

  // Update for DiscreteTransferFcn: '<S55>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states = K;
}

//
// System initialize for atomic system:
//    '<S21>/pidWithDebug'
//    '<S62>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_Init(DW_pidWithDebug_fcsModel_T *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S24>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S21>/pidWithDebug'
//    '<S62>/pidWithDebug'
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

  // Product: '<S56>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S56>/sample time'
  //
  //  About '<S56>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[1] * 0.004;

  // Sum: '<S24>/Sum'
  rtb_Sum_k = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S24>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_k,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_o, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S24>/Discrete First Order Deriv Filter'

  // Product: '<S24>/Product'
  rtb_Product5_o *= rtu_pidParamBus->Kd;

  // Product: '<S24>/Product1'
  rtb_UnitDelay_i = rtb_Sum_k * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S24>/Discrete-Time Integrator'
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S24>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S24>/Discrete-Time Integrator'

  rtb_Sum1_b = ((rtu_feedForward + rtb_Product5_o) + rtb_UnitDelay_i) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S57>/Switch2' incorporates:
  //   RelationalOperator: '<S57>/LowerRelop1'
  //   RelationalOperator: '<S57>/UpperRelop'
  //   Switch: '<S57>/Switch'

  if (rtb_Sum1_b > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_b < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S57>/Switch'
    rtb_Switch2_oz = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_oz = rtb_Sum1_b;
  }

  // End of Switch: '<S57>/Switch2'

  // Sum: '<S56>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S56>/Delay Input2'
  //
  //  Block description for '<S56>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S56>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1_k = rtb_Switch2_oz - localDW->DelayInput2_DSTATE;

  // Switch: '<S59>/Switch2' incorporates:
  //   RelationalOperator: '<S59>/LowerRelop1'

  if (rtb_UkYk1_k <= rtb_Switch2_p) {
    // Product: '<S56>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S56>/sample time'
    //
    //  About '<S56>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_p = rtu_pidParamBus->outputRateLimits[0] * 0.004;

    // Switch: '<S59>/Switch' incorporates:
    //   RelationalOperator: '<S59>/UpperRelop'

    if (rtb_UkYk1_k >= rtb_Switch2_p) {
      rtb_Switch2_p = rtb_UkYk1_k;
    }

    // End of Switch: '<S59>/Switch'
  }

  // End of Switch: '<S59>/Switch2'

  // Sum: '<S56>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S56>/Delay Input2'
  //
  //  Block description for '<S56>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S56>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Switch2_p + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S24>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S24>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_i;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_o;

  // Update for DiscreteIntegrator: '<S24>/Discrete-Time Integrator' incorporates:
  //   Product: '<S24>/Product2'
  //   Product: '<S24>/Product3'
  //   Product: '<S24>/Product5'
  //   Sum: '<S24>/Sum2'
  //   Sum: '<S24>/Sum3'
  //   Sum: '<S24>/Sum4'
  //   Sum: '<S24>/Sum5'
  //   UnitDelay: '<S24>/Unit Delay'
  //   UnitDelay: '<S24>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_k * rtu_pidParamBus->Ki) * 0.004;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S56>/Delay Input2'
  //
  //  Block description for '<S56>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S24>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_oz;

  // Update for UnitDelay: '<S24>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_b;
}

//
// Output and update for atomic system:
//    '<S40>/Compute Natural Frequency'
//    '<S41>/Compute Natural Frequency'
//    '<S25>/Compute Natural Frequency'
//    '<S26>/Compute Natural Frequency'
//    '<S84>/Compute Natural Frequency'
//    '<S85>/Compute Natural Frequency'
//    '<S69>/Compute Natural Frequency'
//    '<S70>/Compute Natural Frequency'
//    '<S146>/Compute Natural Frequency'
//    '<S147>/Compute Natural Frequency'
//    ...
//
void fcsModel::fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps,
  real_T rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps)
{
  real_T tmp;

  //  call the main function
  // MATLAB Function 'Discrete Second Order Filter/Compute Natural Frequency': '<S48>:1' 
  // '<S48>:1:4' naturalFrequency_radps = computeSecondOrderSystemNaturalFrequency_function(bandwidth_radps, dampingRatio_nd); 
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
//    '<S40>/Compute Numerator And Denominator'
//    '<S25>/Compute Numerator And Denominator'
//    '<S84>/Compute Numerator And Denominator'
//    '<S69>/Compute Numerator And Denominator'
//    '<S146>/Compute Numerator And Denominator'
//    '<S131>/Compute Numerator And Denominator'
//    '<S186>/Compute Numerator And Denominator'
//    '<S201>/Compute Numerator And Denominator'
//    '<S255>/Compute Numerator And Denominator'
//    '<S240>/Compute Numerator And Denominator'
//    '<S225>/Compute Numerator And Denominator'
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
  // MATLAB Function 'Discrete Second Order Deriv Filter/Compute Numerator And Denominator': '<S49>:1' 
  // '<S49>:1:4' [rateNum, accelNum, den] = computeSecondOrderDerivFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
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
//    '<S41>/Compute Filter Numerator And Denominator'
//    '<S26>/Compute Filter Numerator And Denominator'
//    '<S85>/Compute Filter Numerator And Denominator'
//    '<S70>/Compute Filter Numerator And Denominator'
//    '<S147>/Compute Filter Numerator And Denominator'
//    '<S132>/Compute Filter Numerator And Denominator'
//    '<S187>/Compute Filter Numerator And Denominator'
//    '<S202>/Compute Filter Numerator And Denominator'
//    '<S256>/Compute Filter Numerator And Denominator'
//    '<S241>/Compute Filter Numerator And Denominator'
//    '<S226>/Compute Filter Numerator And Denominator'
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
//    '<S41>/Compute Filter Numerator And Denominator'
//    '<S26>/Compute Filter Numerator And Denominator'
//    '<S85>/Compute Filter Numerator And Denominator'
//    '<S70>/Compute Filter Numerator And Denominator'
//    '<S147>/Compute Filter Numerator And Denominator'
//    '<S132>/Compute Filter Numerator And Denominator'
//    '<S187>/Compute Filter Numerator And Denominator'
//    '<S202>/Compute Filter Numerator And Denominator'
//    '<S256>/Compute Filter Numerator And Denominator'
//    '<S241>/Compute Filter Numerator And Denominator'
//    '<S226>/Compute Filter Numerator And Denominator'
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
  // MATLAB Function 'Discrete Second Order Filter/Compute Filter Numerator And Denominator': '<S50>:1' 
  // '<S50>:1:4' [num, den] = computeSecondOrderFilterNumAndDen_function(naturalFrequency_radps, dampingRatio_nd, sampleTime_s); 
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
//    '<S21>/Signal Conditioning Block1'
//    '<S21>/Signal Conditioning Block'
//    '<S62>/Signal Conditioning Block1'
//    '<S62>/Signal Conditioning Block'
//
void fcsModel::f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S41>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S21>/Signal Conditioning Block1'
//    '<S21>/Signal Conditioning Block'
//    '<S62>/Signal Conditioning Block1'
//    '<S62>/Signal Conditioning Block'
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

  // MATLAB Function: '<S40>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S40>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum_b[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S41>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_h);

  // MATLAB Function: '<S41>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2_h,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S41>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_j = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S45>/Switch2' incorporates:
  //   RelationalOperator: '<S45>/LowerRelop1'
  //   RelationalOperator: '<S45>/UpperRelop'
  //   Switch: '<S45>/Switch'

  if (rtb_DiscreteTransferFcn_j > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_j < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S45>/Switch'
    rtb_DiscreteTransferFcn_j = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S45>/Switch2'

  // Sum: '<S42>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S42>/Delay Input2'
  //
  //  Block description for '<S42>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S42>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_j -= localDW->DelayInput2_DSTATE;

  // Switch: '<S52>/Switch2' incorporates:
  //   Product: '<S42>/delta rise limit'
  //   SampleTimeMath: '<S42>/sample time'
  //
  //  About '<S42>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_h = rtu_params->filteredInputRateLimits[1] * 0.004;

  // Switch: '<S52>/Switch2' incorporates:
  //   RelationalOperator: '<S52>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_j <= rtb_Switch2_h) {
    // Product: '<S42>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S42>/sample time'
    //
    //  About '<S42>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_h = rtu_params->filteredInputRateLimits[0] * 0.004;

    // Switch: '<S52>/Switch' incorporates:
    //   RelationalOperator: '<S52>/UpperRelop'

    if (rtb_DiscreteTransferFcn_j >= rtb_Switch2_h) {
      // Switch: '<S52>/Switch2'
      rtb_Switch2_h = rtb_DiscreteTransferFcn_j;
    }

    // End of Switch: '<S52>/Switch'
  }

  // End of Switch: '<S52>/Switch2'

  // Sum: '<S42>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S42>/Delay Input2'
  //
  //  Block description for '<S42>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S42>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2_h + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S41>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S42>/Delay Input2'
  //
  //  Block description for '<S42>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_filteredInput;
}

//
// Output and update for atomic system:
//    '<S113>/holdOutputAtCenter1'
//    '<S113>/holdOutputAtCenter2'
//
void fcsModel::fcsModel_holdOutputAtCenter1(real_T rtu_input, real_T rtu_trigger,
  real_T *rty_output, boolean_T *rty_atCenter, DW_holdOutputAtCenter1_fcsMod_T
  *localDW)
{
  // MATLAB Function: '<S123>/holdOutputAtCenter'
  // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S126>:1'
  // '<S126>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
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

  // End of MATLAB Function: '<S123>/holdOutputAtCenter'
}

//
// System initialize for atomic system:
//    '<S114>/pidWithDebug'
//    '<S182>/pidWithDebug'
//    '<S167>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_m_Init(DW_pidWithDebug_fcsModel_i_T
  *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S130>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S114>/pidWithDebug'
//    '<S182>/pidWithDebug'
//    '<S167>/pidWithDebug'
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

  // Product: '<S162>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S162>/sample time'
  //
  //  About '<S162>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_d3 = rtu_pidParamBus->outputRateLimits[1] * 0.008;

  // Sum: '<S130>/Sum'
  rtb_Sum_c4 = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S130>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_c4,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_e, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S130>/Discrete First Order Deriv Filter'

  // Product: '<S130>/Product'
  rtb_Product5_e *= rtu_pidParamBus->Kd;

  // Product: '<S130>/Product1'
  rtb_UnitDelay_a = rtb_Sum_c4 * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S130>/Discrete-Time Integrator'
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S130>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S130>/Discrete-Time Integrator'

  rtb_Sum1_o = ((rtu_feedForward + rtb_Product5_e) + rtb_UnitDelay_a) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S163>/Switch2' incorporates:
  //   RelationalOperator: '<S163>/LowerRelop1'
  //   RelationalOperator: '<S163>/UpperRelop'
  //   Switch: '<S163>/Switch'

  if (rtb_Sum1_o > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_b = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1_o < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S163>/Switch'
    rtb_Switch2_b = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_b = rtb_Sum1_o;
  }

  // End of Switch: '<S163>/Switch2'

  // Sum: '<S162>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S162>/Delay Input2'
  //
  //  Block description for '<S162>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S162>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1_h = rtb_Switch2_b - localDW->DelayInput2_DSTATE;

  // Switch: '<S165>/Switch2' incorporates:
  //   RelationalOperator: '<S165>/LowerRelop1'

  if (rtb_UkYk1_h <= rtb_Switch2_d3) {
    // Product: '<S162>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S162>/sample time'
    //
    //  About '<S162>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_d3 = rtu_pidParamBus->outputRateLimits[0] * 0.008;

    // Switch: '<S165>/Switch' incorporates:
    //   RelationalOperator: '<S165>/UpperRelop'

    if (rtb_UkYk1_h >= rtb_Switch2_d3) {
      rtb_Switch2_d3 = rtb_UkYk1_h;
    }

    // End of Switch: '<S165>/Switch'
  }

  // End of Switch: '<S165>/Switch2'

  // Sum: '<S162>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S162>/Delay Input2'
  //
  //  Block description for '<S162>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S162>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Switch2_d3 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S130>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S130>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_a;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_e;

  // Update for DiscreteIntegrator: '<S130>/Discrete-Time Integrator' incorporates:
  //   Product: '<S130>/Product2'
  //   Product: '<S130>/Product3'
  //   Product: '<S130>/Product5'
  //   Sum: '<S130>/Sum2'
  //   Sum: '<S130>/Sum3'
  //   Sum: '<S130>/Sum4'
  //   Sum: '<S130>/Sum5'
  //   UnitDelay: '<S130>/Unit Delay'
  //   UnitDelay: '<S130>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_c4 * rtu_pidParamBus->Ki) * 0.008;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S162>/Delay Input2'
  //
  //  Block description for '<S162>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S130>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_b;

  // Update for UnitDelay: '<S130>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1_o;
}

//
// System initialize for atomic system:
//    '<S114>/Signal Conditioning Block1'
//    '<S114>/Signal Conditioning Block'
//    '<S182>/Signal Conditioning Block'
//    '<S167>/Signal Conditioning Block2'
//    '<S167>/Signal Conditioning Block1'
//    '<S167>/Signal Conditioning Block'
//
void fcsModel::SignalConditioningBlock1_c_Init(DW_SignalConditioningBlock1_g_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S147>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S114>/Signal Conditioning Block1'
//    '<S114>/Signal Conditioning Block'
//    '<S182>/Signal Conditioning Block'
//    '<S167>/Signal Conditioning Block2'
//    '<S167>/Signal Conditioning Block1'
//    '<S167>/Signal Conditioning Block'
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

  // MATLAB Function: '<S146>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_n);

  // MATLAB Function: '<S146>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2_n,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum_k[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S147>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2_n);

  // MATLAB Function: '<S147>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2_n,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S147>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_d = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S151>/Switch2' incorporates:
  //   RelationalOperator: '<S151>/LowerRelop1'
  //   RelationalOperator: '<S151>/UpperRelop'
  //   Switch: '<S151>/Switch'

  if (rtb_DiscreteTransferFcn_d > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_d < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S151>/Switch'
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S151>/Switch2'

  // Sum: '<S148>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S148>/Delay Input2'
  //
  //  Block description for '<S148>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S148>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_d -= localDW->DelayInput2_DSTATE;

  // Switch: '<S158>/Switch2' incorporates:
  //   Product: '<S148>/delta rise limit'
  //   SampleTimeMath: '<S148>/sample time'
  //
  //  About '<S148>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2_n = rtu_params->filteredInputRateLimits[1] * 0.008;

  // Switch: '<S158>/Switch2' incorporates:
  //   RelationalOperator: '<S158>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_d <= rtb_Switch2_n) {
    // Product: '<S148>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S148>/sample time'
    //
    //  About '<S148>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2_n = rtu_params->filteredInputRateLimits[0] * 0.008;

    // Switch: '<S158>/Switch' incorporates:
    //   RelationalOperator: '<S158>/UpperRelop'

    if (rtb_DiscreteTransferFcn_d >= rtb_Switch2_n) {
      // Switch: '<S158>/Switch2'
      rtb_Switch2_n = rtb_DiscreteTransferFcn_d;
    }

    // End of Switch: '<S158>/Switch'
  }

  // End of Switch: '<S158>/Switch2'

  // Sum: '<S148>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S148>/Delay Input2'
  //
  //  Block description for '<S148>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S148>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2_n + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S147>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S148>/Delay Input2'
  //
  //  Block description for '<S148>/Delay Input2':
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

  // MATLAB Function 'checkRcCmds': '<S275>:7'
  // '<S275>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S275>:7:3' if(rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle(1) && ... 
  // '<S275>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S275>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S275>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 1000) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 1000) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 1000) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 1000) {
          // '<S275>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S275>:7:8' else
          // '<S275>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S275>:7:8' else
        // '<S275>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S275>:7:8' else
      // '<S275>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S275>:7:8' else
    // '<S275>:7:9' isTrue = false;
    isTrue = false;
  }

  return isTrue;
}

real_T rt_urand_Upu32_Yd_f_pw(uint32_T *u)
{
  uint32_T hi;
  uint32_T lo;

  // Uniform random number generator (random number between 0 and 1)

  // #define IA      16807                      magic multiplier = 7^5
  // #define IM      2147483647                 modulus = 2^31-1
  // #define IQ      127773                     IM div IA
  // #define IR      2836                       IM modulo IA
  // #define S       4.656612875245797e-10      reciprocal of 2^31-1
  // test = IA * (seed % IQ) - IR * (seed/IQ)
  // seed = test < 0 ? (test + IM) : test
  // return (seed*S)

  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return static_cast<real_T>(*u) * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw(uint32_T *u)
{
  real_T si;
  real_T sr;
  real_T y;

  // Normal (Gaussian) random number generator
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = std::sqrt(-2.0 * std::log(si) / si) * sr;
  return y;
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
  std::array<real_T, 4> DiscreteTransferFcn_tmp_b;
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
  real_T rtb_DifferenceInputs2_i;
  real_T rtb_DifferenceInputs2_jb;
  real_T rtb_DiscreteTransferFcn;
  real_T rtb_DiscreteTransferFcn_k;
  real_T rtb_frcCmd_N;
  real_T tmp;
  real_T tmp_0;
  real_T tmp_1;
  real_T tmp_2;
  real_T vxCmd_unitRange;
  real_T vyCmd_unitRange;
  real_T yCmd;
  int32_T pCmd;
  int32_T rCmd;
  int32_T tCmd;
  boolean_T resetIntegrator;
  boolean_T rtb_Compare_d;
  boolean_T rtb_Compare_mi;
  enumChirpTrigger rtb_chirpTrigger;
  enumChirpType rtb_chirpType;
  enumFlightMode flightMode;
  enumStateMachine state;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // Math: '<S113>/Transpose' incorporates:
    //   Inport: '<Root>/stateEstimate'

    tCmd = 0;
    for (rCmd = 0; rCmd < 3; rCmd++) {
      rtb_Transpose[tCmd] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd];
      rtb_Transpose[tCmd + 1] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd + 3];
      rtb_Transpose[tCmd + 2] = fcsModel_U.stateEstimate.ned2FepDcm_nd[rCmd + 6];
      tCmd += 3;
    }

    // End of Math: '<S113>/Transpose'
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
    // Transition: '<S275>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S275>:1'
    // '<S275>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S275>:1:3' rcCheckFlag = checkRcCmds;
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
    }

    // '<S275>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S275>:3'
      // '<S275>:10:1' sf_internal_predicateOutput = after(60, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 15000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S275>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S275>:1'
        // '<S275>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S275>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S275>:1:4' resetIntegrator = true;

        // '<S275>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimitsThrottle(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 1000) {
        // Transition: '<S275>:12'
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S275>:11'
        // '<S275>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S275>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S275>:11:4' resetIntegrator = false;
      } else {
        // '<S275>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S275>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S275>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S275>:1'
      // '<S275>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 1 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S275>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S275>:3'
        // '<S275>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S275>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S275>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S275>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S275>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S275>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S275>:11'
      // '<S275>:20:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle (1); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1000) {
        // Transition: '<S275>:20'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S275>:3'
        // '<S275>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S275>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S275>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S275>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S275>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S275>:11:4' resetIntegrator = false;
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
  //   UnitDelay: '<S4>/Unit Delay'

  // Computes command and flight mode from the rc inputs
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S276>:1'
  // '<S276>:1:3' [flightMode, rcOutCmds, chirpTrigger, chirpType] = interpretRcInputs_function(rcCmds, prevChirpTrigger, rcParamsStruct); 
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
  // 'interpretRcInputs_function:24' if isempty(throttle_is_up)
  // 'interpretRcInputs_function:29' if(rcInCmds.rcSwitch3_nd >= 1500)
  if (fcsModel_U.rcCmdsIn.rcSwitch3_nd >= 1500) {
    // 'interpretRcInputs_function:30' if(prevChirpTrigger == enumChirpTrigger.OFF) 
    if (fcsModel_DW.UnitDelay_DSTATE_g == enumChirpTrigger::OFF) {
      // 'interpretRcInputs_function:31' chirpCount_ = chirpCount_ + uint8(1);
      tCmd = static_cast<int32_T>(fcsModel_DW.chirpCount_ + 1U);
      if (fcsModel_DW.chirpCount_ + 1U > 255U) {
        tCmd = 255;
      }

      fcsModel_DW.chirpCount_ = static_cast<uint8_T>(tCmd);
    }

    // 'interpretRcInputs_function:33' chirpTrigger = enumChirpTrigger.ON;
    rtb_chirpTrigger = enumChirpTrigger::ON;
  } else {
    // 'interpretRcInputs_function:34' else
    // 'interpretRcInputs_function:35' if(chirpCount_ == enumChirpType.FZ)
    if (fcsModel_DW.chirpCount_ == static_cast<int32_T>(enumChirpType::FZ)) {
      // 'interpretRcInputs_function:36' chirpCount_ = uint8(0);
      fcsModel_DW.chirpCount_ = 0U;
    }

    // 'interpretRcInputs_function:38' chirpTrigger = enumChirpTrigger.OFF;
    rtb_chirpTrigger = enumChirpTrigger::OFF;
  }

  // 'interpretRcInputs_function:41' switch chirpCount_
  switch (fcsModel_DW.chirpCount_) {
   case 0U:
    // 'interpretRcInputs_function:42' case 0
    // 'interpretRcInputs_function:43' chirpType = enumChirpType.NONE;
    rtb_chirpType = enumChirpType::NONE;
    break;

   case 1U:
    // 'interpretRcInputs_function:44' case 1
    // 'interpretRcInputs_function:45' chirpType = enumChirpType.MX;
    rtb_chirpType = enumChirpType::MX;
    break;

   case 2U:
    // 'interpretRcInputs_function:46' case 2
    // 'interpretRcInputs_function:47' chirpType = enumChirpType.MY;
    rtb_chirpType = enumChirpType::MY;
    break;

   case 3U:
    // 'interpretRcInputs_function:48' case 3
    // 'interpretRcInputs_function:49' chirpType = enumChirpType.MZ;
    rtb_chirpType = enumChirpType::MZ;
    break;

   case 4U:
    // 'interpretRcInputs_function:50' case 4
    // 'interpretRcInputs_function:51' chirpType = enumChirpType.FZ;
    rtb_chirpType = enumChirpType::FZ;
    break;

   default:
    // 'interpretRcInputs_function:52' otherwise
    // 'interpretRcInputs_function:53' chirpType = enumChirpType.NONE;
    rtb_chirpType = enumChirpType::NONE;
    break;
  }

  //  Used to directly set force commands in N that is fed into allocation in
  //  STABILIZE and ACRO flight modes
  // 'interpretRcInputs_function:59' rcOutCmds.throttleStick = 0;
  //  Used to command roll and pitch angles in STABILIZE flight mode
  // 'interpretRcInputs_function:62' rcOutCmds.rollStick = 0;
  // 'interpretRcInputs_function:63' rcOutCmds.pitchStick = 0;
  //  Used to control yaw rate in all flight modes
  // 'interpretRcInputs_function:66' rcOutCmds.yawStick = 0;
  //  Not used currently
  // 'interpretRcInputs_function:69' rcOutCmds.vxStick_mps = 0;
  // 'interpretRcInputs_function:70' rcOutCmds.vyStick_mps = 0;
  //  Used to command vertical velocity in ALT_CONTROL flight mode. In this
  //  mode rcOutCmds.throttleStick is ignored.
  // 'interpretRcInputs_function:74' rcOutCmds.vzStick_mps = 0;
  // Select mode and set the command limits
  // 'interpretRcInputs_function:77' if(rcInCmds.rcSwitch1_nd < 1100)
  if (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1100) {
    // 'interpretRcInputs_function:78' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:79' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:80' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:81' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:82' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:83' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
  } else if ((fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1100) &&
             (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1700)) {
    // 'interpretRcInputs_function:85' elseif (rcInCmds.rcSwitch1_nd >= 1100 && rcInCmds.rcSwitch1_nd < 1700) 
    // 'interpretRcInputs_function:86' flightMode = enumFlightMode.ALT_CONTROL;
    flightMode = enumFlightMode::ALT_CONTROL;

    // 'interpretRcInputs_function:87' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:88' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:89' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:90' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:91' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else if (fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1700) {
    // 'interpretRcInputs_function:92' elseif (rcInCmds.rcSwitch1_nd >= 1700)
    // 'interpretRcInputs_function:93' flightMode = enumFlightMode.POS_CONTROL;
    flightMode = enumFlightMode::POS_CONTROL;

    // 'interpretRcInputs_function:94' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:95' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:96' rlim = 0;
    rlim = 0.0;

    // 'interpretRcInputs_function:97' plim = 0;
    plim = 0.0;

    // 'interpretRcInputs_function:98' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else {
    // 'interpretRcInputs_function:99' else
    // 'interpretRcInputs_function:100' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:101' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:102' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:103' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:104' vxlim = rcParamsStruct.cmdLimits.vx_mps(2); 
    // 'interpretRcInputs_function:105' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
  }

  // 'interpretRcInputs_function:108' if double(rcInCmds.throttleCmd_nd) < rcParamsStruct.pwmLimits(1) 
  if (fcsModel_U.rcCmdsIn.throttleCmd_nd < 1000) {
    // This means either we haven't take off yet or we landed after a flight
    // and might take of again so set the throttle_is_up to false
    // 'interpretRcInputs_function:111' throttle_is_up = false;
    fcsModel_DW.throttle_is_up = false;
  } else if ((fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1450) &&
             (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 1300)) {
    // 'interpretRcInputs_function:112' elseif ((double(rcInCmds.throttleCmd_nd) <= rcParamsStruct.pwmThrottleMidHigh) && ... 
    // 'interpretRcInputs_function:113'         double(rcInCmds.throttleCmd_nd) >= rcParamsStruct.pwmThrottleMidLow) 
    // Since throttle stick on RC transmitter is half way up this probably
    // means the vehicle is flying so set throttle_is_up flag to true
    // 'interpretRcInputs_function:116' throttle_is_up = true;
    fcsModel_DW.throttle_is_up = true;
  } else {
    // 'interpretRcInputs_function:117' else
    // Otherwise preserve the previous flag
  }

  // 'interpretRcInputs_function:121' tCmd = min( rcParamsStruct.pwmLimitsThrottle(2), ... 
  // 'interpretRcInputs_function:122'     max( rcParamsStruct.pwmLimitsThrottle(1), double(rcInCmds.throttleCmd_nd) ) ); 
  tCmd = static_cast<int32_T>(std::fmin(1900.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.throttleCmd_nd))));

  //  In stabilize mode throttle stick starts at 0
  // 'interpretRcInputs_function:125' tCmd_unitRange = (rcParamsStruct.throttleUnitRangeMapCoeff.a*tCmd + ... 
  // 'interpretRcInputs_function:126'     rcParamsStruct.throttleUnitRangeMapCoeff.c)/rcParamsStruct.throttleUnitRangeMapCoeff.b; 
  //  Use throttle stick to set Vz commands. Always set but only used
  //  in ALT_CONTROL flight mode. Has different slopes about the center point
  //  as the center point is not always at 1500 which is PWM center
  // 'interpretRcInputs_function:131' if ((tCmd <= rcParamsStruct.pwmThrottleMidHigh) && ... 
  // 'interpretRcInputs_function:132'         tCmd >= rcParamsStruct.pwmThrottleMidLow) 
  if ((tCmd <= 1450) && (tCmd >= 1300)) {
    // 'interpretRcInputs_function:133' rcOutCmds.vzStick_mps = 0;
    fcsModel_DW.rcOutCmds.vzStick_mps = 0.0;
  } else if (tCmd < 1300) {
    // 'interpretRcInputs_function:134' elseif (tCmd < rcParamsStruct.pwmThrottleMidLow) 
    // 'interpretRcInputs_function:135' if (throttle_is_up)
    if (fcsModel_DW.throttle_is_up) {
      //  This means we are already flying, make lower half of throttle
      //  stick map to descent rates
      // 'interpretRcInputs_function:138' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:139'             rcParamsStruct.pwmToCmdThrottleIncptLow; 
      // ;
      // 'interpretRcInputs_function:140' rcOutCmds.vzStick_mps = (rcParamsStruct.vzLowRangeMapCoeff.a*vzCmd_unitRange + ... 
      // 'interpretRcInputs_function:141'             rcParamsStruct.vzLowRangeMapCoeff.c)/ ... 
      // 'interpretRcInputs_function:142'             rcParamsStruct.vzLowRangeMapCoeff.b; 
      fcsModel_DW.rcOutCmds.vzStick_mps = -(0.0033333333333333335 * static_cast<
        real_T>(tCmd) + -4.333333333333333);
    } else {
      // 'interpretRcInputs_function:143' else
      //  This means we haven't taken off yet or we landed and might take
      //  off again. Map the lower half of the throttle stick also to
      //  ascent rate as we want to climb up from the ground and only go
      //  upto 3/4 of the max ascent rate at mid point. Lowest position is
      //  throttle maps to -1 so add +1 make lowest throttle point 0 and
      //  increase from there on
      // 'interpretRcInputs_function:150' vzCmd_unitRange = ((rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:151'             rcParamsStruct.pwmToCmdThrottleIncptLow) + 1)*0.75; 
      // ;
      // 'interpretRcInputs_function:152' rcOutCmds.vzStick_mps = (rcParamsStruct.vzHighRangeMapCoeff.a*vzCmd_unitRange + ... 
      // 'interpretRcInputs_function:153'             rcParamsStruct.vzHighRangeMapCoeff.c)/ ... 
      // 'interpretRcInputs_function:154'             rcParamsStruct.vzHighRangeMapCoeff.b; 
      fcsModel_DW.rcOutCmds.vzStick_mps = ((0.0033333333333333335 * static_cast<
        real_T>(tCmd) + -4.333333333333333) + 1.0) * 0.75 * -1.5;
    }
  } else {
    // 'interpretRcInputs_function:156' else
    // 'interpretRcInputs_function:157' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeHigh*tCmd + ...; 
    // 'interpretRcInputs_function:158'         rcParamsStruct.pwmToCmdThrottleIncptHigh; 
    // ;
    // 'interpretRcInputs_function:159' rcOutCmds.vzStick_mps = (rcParamsStruct.vzHighRangeMapCoeff.a*vzCmd_unitRange + ... 
    // 'interpretRcInputs_function:160'         rcParamsStruct.vzHighRangeMapCoeff.c)/ ... 
    // 'interpretRcInputs_function:161'         rcParamsStruct.vzHighRangeMapCoeff.b; 
    fcsModel_DW.rcOutCmds.vzStick_mps = (0.0022222222222222222 *
      static_cast<real_T>(tCmd) + -3.2222222222222223) * -1.5;
  }

  //  Set roll, pitch and yaw stick
  // 'interpretRcInputs_function:165' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:166'     max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  rCmd = static_cast<int32_T>(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickXCmd_nd))));

  //  Use roll stick to set FEP Vy to be used for POS control mode
  // 'interpretRcInputs_function:169' if ((rCmd <= rcParamsStruct.pwmRollStickMidHigh) && ... 
  // 'interpretRcInputs_function:170'         rCmd >= rcParamsStruct.pwmRollStickMidLow) 
  if ((rCmd <= 1650) && (rCmd >= 1350)) {
    // 'interpretRcInputs_function:171' vyCmd_unitRange  = 0;
    vyCmd_unitRange = 0.0;
  } else {
    // 'interpretRcInputs_function:172' else
    // 'interpretRcInputs_function:173' vyCmd_unitRange = -3 + rCmd/500;
    vyCmd_unitRange = static_cast<real_T>(rCmd) / 500.0 + -3.0;
  }

  // 'interpretRcInputs_function:175' rCmd_unitRange = -3 + rCmd/500;
  // 'interpretRcInputs_function:178' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:179'     max( rcParamsStruct.pwmLimits(1),  double(rcInCmds.joystickYCmd_nd) ) ); 
  pCmd = static_cast<int32_T>(std::fmin(2000.0, std::fmax(1000.0,
    static_cast<real_T>(fcsModel_U.rcCmdsIn.joystickYCmd_nd))));

  //  Use pitch stick to set FEP Vx to be used for POS control mode
  // 'interpretRcInputs_function:182' if ((pCmd <= rcParamsStruct.pwmPitchStickMidHigh) && ... 
  // 'interpretRcInputs_function:183'         pCmd >= rcParamsStruct.pwmPitchStickMidLow) 
  if ((pCmd <= 1650) && (pCmd >= 1350)) {
    // 'interpretRcInputs_function:184' vxCmd_unitRange  = 0;
    vxCmd_unitRange = 0.0;
  } else {
    // 'interpretRcInputs_function:185' else
    // 'interpretRcInputs_function:186' vxCmd_unitRange = -3 + pCmd/500;
    vxCmd_unitRange = static_cast<real_T>(pCmd) / 500.0 + -3.0;
  }

  //  Reverse the pitch cmd
  // 'interpretRcInputs_function:190' pCmd_unitRange = -(-3 + pCmd/500);
  // 'interpretRcInputs_function:193' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:194'     max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  yCmd = std::fmin(2000.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickZCmd_nd)));

  //  Use yaw stick to also pick a Yaw angle in ALT or POS control mode
  // 'interpretRcInputs_function:197' if (flightMode == enumFlightMode.STABILIZE) 
  if (flightMode == enumFlightMode::STABILIZE) {
    // 'interpretRcInputs_function:198' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;

    // 'interpretRcInputs_function:199' else
    // 'interpretRcInputs_function:200' if ((yCmd <= rcParamsStruct.pwmYawStickMidHigh) && ... 
    // 'interpretRcInputs_function:201'             yCmd >= rcParamsStruct.pwmYawStickMidLow) 
  } else if ((yCmd <= 1600.0) && (yCmd >= 1400.0)) {
    // 'interpretRcInputs_function:202' yCmd_unitRange  = 0;
    yCmd = 0.0;
  } else {
    // 'interpretRcInputs_function:203' else
    // 'interpretRcInputs_function:204' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;
  }

  //  Usually expo is set in the Tx hence simply use a linear map here
  // 'interpretRcInputs_function:209' rcOutCmds.throttleStick = (rcParamsStruct.fCmdRangeMapCoeff.a*tCmd_unitRange + ... 
  // 'interpretRcInputs_function:210'     rcParamsStruct.fCmdRangeMapCoeff.c)/ ... 
  // 'interpretRcInputs_function:211'     rcParamsStruct.fCmdRangeMapCoeff.b;
  fcsModel_DW.rcOutCmds.throttleStick = (static_cast<real_T>(tCmd) + -1000.0) /
    900.0 * -20.0 + -15.0;

  // 'interpretRcInputs_function:212' rcOutCmds.rollStick = rCmd_unitRange*rlim; 
  fcsModel_DW.rcOutCmds.rollStick = (static_cast<real_T>(rCmd) / 500.0 + -3.0) *
    rlim;

  // 'interpretRcInputs_function:213' rcOutCmds.pitchStick = pCmd_unitRange*plim; 
  fcsModel_DW.rcOutCmds.pitchStick = -(static_cast<real_T>(pCmd) / 500.0 + -3.0)
    * plim;

  // 'interpretRcInputs_function:214' rcOutCmds.yawStick = yCmd_unitRange*ylim;
  fcsModel_DW.rcOutCmds.yawStick = yCmd * 1.0471975511965976;

  // 'interpretRcInputs_function:215' rcOutCmds.vxStick_mps = vxCmd_unitRange*vxlim; 
  fcsModel_DW.rcOutCmds.vxStick_mps = vxCmd_unitRange * 2.5;

  // 'interpretRcInputs_function:216' rcOutCmds.vyStick_mps = vyCmd_unitRange*vylim; 
  fcsModel_DW.rcOutCmds.vyStick_mps = vyCmd_unitRange * 2.5;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // SignalConversion generated from: '<S113>/Vector Concatenate1'
    rtb_VectorConcatenate1[0] = fcsModel_DW.rcOutCmds.vxStick_mps;

    // SignalConversion generated from: '<S113>/Vector Concatenate1'
    rtb_VectorConcatenate1[1] = fcsModel_DW.rcOutCmds.vyStick_mps;

    // Product: '<S113>/Matrix Multiply' incorporates:
    //   Math: '<S113>/Transpose'

    for (tCmd = 0; tCmd < 3; tCmd++) {
      rtb_MatrixMultiply[tCmd] = 0.0;
      rtb_MatrixMultiply[tCmd] += rtb_Transpose[tCmd] * rtb_VectorConcatenate1[0];
      rtb_MatrixMultiply[tCmd] += rtb_Transpose[tCmd + 3] *
        rtb_VectorConcatenate1[1];
    }

    // End of Product: '<S113>/Matrix Multiply'

    // Outputs for Atomic SubSystem: '<S113>/holdOutputAtCenter1'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[0],
      rtb_MatrixMultiply[0], &rtb_DifferenceInputs2_jb, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter1);

    // End of Outputs for SubSystem: '<S113>/holdOutputAtCenter1'

    // Concatenate: '<S113>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // Switch: '<S113>/Switch1' incorporates:
    //   RelationalOperator: '<S117>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S113>/Bus Assignment1' incorporates:
      //   Concatenate: '<S113>/Vector Concatenate'
      //   Constant: '<S113>/Constant1'

      rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S113>/Bus Assignment1' incorporates:
      //   Concatenate: '<S113>/Vector Concatenate'

      rtb_VectorConcatenate[0].feedForwardCmd = rtb_MatrixMultiply[0];
    }

    // End of Switch: '<S113>/Switch1'

    // BusAssignment: '<S113>/Bus Assignment1' incorporates:
    //   Concatenate: '<S113>/Vector Concatenate'
    //   Constant: '<S120>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S113>/Logical Operator2'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S120>/Compare'

    rtb_VectorConcatenate[0].cmd = rtb_DifferenceInputs2_jb;
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedPos_m[0];
    rtb_VectorConcatenate[0].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S113>/holdOutputAtCenter2'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[1],
      rtb_MatrixMultiply[1], &rtb_DifferenceInputs2_jb, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter2);

    // End of Outputs for SubSystem: '<S113>/holdOutputAtCenter2'

    // Concatenate: '<S113>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // Switch: '<S113>/Switch2' incorporates:
    //   RelationalOperator: '<S118>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S113>/Bus Assignment2' incorporates:
      //   Concatenate: '<S113>/Vector Concatenate'
      //   Constant: '<S113>/Constant3'

      rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S113>/Bus Assignment2' incorporates:
      //   Concatenate: '<S113>/Vector Concatenate'

      rtb_VectorConcatenate[1].feedForwardCmd = rtb_MatrixMultiply[1];
    }

    // End of Switch: '<S113>/Switch2'

    // BusAssignment: '<S113>/Bus Assignment2' incorporates:
    //   Concatenate: '<S113>/Vector Concatenate'
    //   Constant: '<S121>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S113>/Logical Operator3'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S121>/Compare'

    rtb_VectorConcatenate[1].cmd = rtb_DifferenceInputs2_jb;
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedPos_m[1];
    rtb_VectorConcatenate[1].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S113>/holdOutputAtCenter'
    // MATLAB Function: '<S122>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S125>:1'
    // '<S125>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
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

    // End of Outputs for SubSystem: '<S113>/holdOutputAtCenter'

    // Switch: '<S113>/Switch' incorporates:
    //   RelationalOperator: '<S115>/Compare'

    // 'holdOutputAtCenter_function:17' output = last_input;
    if (rtb_Compare_d) {
      // Sum: '<S242>/Difference Inputs2' incorporates:
      //   Constant: '<S113>/Constant5'
      //
      //  Block description for '<S242>/Difference Inputs2':
      //
      //   Add in CPU

      rlim = 0.0;
    } else {
      // Sum: '<S242>/Difference Inputs2'
      //
      //  Block description for '<S242>/Difference Inputs2':
      //
      //   Add in CPU

      rlim = fcsModel_DW.rcOutCmds.vzStick_mps;
    }

    // End of Switch: '<S113>/Switch'

    // Outputs for Atomic SubSystem: '<S113>/holdOutputAtCenter'
    // Gain: '<S113>/Gain' incorporates:
    //   MATLAB Function: '<S122>/holdOutputAtCenter'

    rtb_DifferenceInputs2_i = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S113>/holdOutputAtCenter'

    // Gain: '<S113>/Gain1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = -fcsModel_U.stateEstimate.aglEst_m;

    // Concatenate: '<S113>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S113>/Bus Assignment3' incorporates:
    //   Concatenate: '<S113>/Vector Concatenate'
    //   Constant: '<S116>/Constant'
    //   Constant: '<S119>/Constant'
    //   Gain: '<S113>/Gain'
    //   Gain: '<S113>/Gain1'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S113>/Logical Operator'
    //   Logic: '<S113>/Logical Operator1'
    //   MATLAB Function: '<S122>/holdOutputAtCenter'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S116>/Compare'
    //   RelationalOperator: '<S119>/Compare'

    rtb_VectorConcatenate[2].feedForwardCmd = rlim;

    // Outputs for Atomic SubSystem: '<S113>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].cmd = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S113>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.aglEst_m;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && (flightMode != enumFlightMode::
      POS_CONTROL)));

    // Outputs for Iterator SubSystem: '<S109>/NED Position Control' incorporates:
    //   ForEach: '<S114>/For Each'

    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // Outputs for Atomic SubSystem: '<S114>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S114>/ctrlInputs' incorporates:
      //   BusAssignment: '<S113>/Bus Assignment'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S114>/Unit Delay'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_DifferenceInputs2_jb, 0.008,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S114>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S114>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_frcCmd_N, 0.008,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S114>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S114>/pidWithDebug'
      fcsModel_pidWithDebug_j(rtb_VectorConcatenate[ForEach_itr_i].
        feedForwardCmd, rtb_DifferenceInputs2_jb, rtb_frcCmd_N,
        rtb_VectorConcatenate[ForEach_itr_i].integratorReset, 0.0,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.ctrlParamsArray[
        ForEach_itr_i], fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE,
        &rtb_DifferenceInputs2_i, &rtb_BusCreator_em, 0.008,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].pidWithDebug);

      // End of Outputs for SubSystem: '<S114>/pidWithDebug'

      // Update for UnitDelay: '<S114>/Unit Delay'
      fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE =
        rtb_DifferenceInputs2_i;

      // ForEachSliceAssignment generated from: '<S114>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDeb_m[ForEach_itr_i] = rtb_BusCreator_em;

      // ForEachSliceAssignment generated from: '<S114>/neVelCmd_mps'
      rtb_ImpAsg_InsertedFor_neVelCmd[ForEach_itr_i] = rtb_DifferenceInputs2_i;

      // ForEachSliceAssignment generated from: '<S114>/meas'
      rtb_ImpAsg_InsertedFor_meas_at_[ForEach_itr_i] = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S114>/cmd'
      rtb_ImpAsg_InsertedFor_cmd_at_i[ForEach_itr_i] = rtb_DifferenceInputs2_jb;
    }

    // End of Outputs for SubSystem: '<S109>/NED Position Control'

    // RelationalOperator: '<S112>/Compare' incorporates:
    //   Constant: '<S112>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_mi = (flightMode != enumFlightMode::POS_CONTROL);

    // Logic: '<S108>/Logical Operator2'
    rtb_Compare_d = (resetIntegrator || rtb_Compare_mi);

    // Concatenate: '<S108>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S108>/Bus Assignment' incorporates:
    //   Concatenate: '<S108>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[0];
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedVel_mps[0];
    rtb_VectorConcatenate[0].integratorReset = rtb_Compare_d;

    // Concatenate: '<S108>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S108>/Bus Assignment1' incorporates:
    //   Concatenate: '<S108>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[1];
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedVel_mps[1];
    rtb_VectorConcatenate[1].integratorReset = rtb_Compare_d;

    // Concatenate: '<S108>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S108>/Bus Assignment2' incorporates:
    //   Concatenate: '<S108>/Vector Concatenate'
    //   Constant: '<S111>/Constant'
    //   Gain: '<S108>/Gain'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S108>/Logical Operator'
    //   Logic: '<S108>/Logical Operator1'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S111>/Compare'

    rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.climbRateEst_mps;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && rtb_Compare_mi));

    // SignalConversion generated from: '<S110>/For Each Subsystem'
    rtb_VectorConcatenate1[0] = 0.0;
    rtb_VectorConcatenate1[1] = 0.0;
    rtb_VectorConcatenate1[2] = 0.0;

    // Outputs for Iterator SubSystem: '<S110>/For Each Subsystem' incorporates:
    //   ForEach: '<S167>/For Each'

    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // Outputs for Atomic SubSystem: '<S167>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S167>/ctrlInputs' incorporates:
      //   BusAssignment: '<S108>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr], &rtb_DifferenceInputs2_i, 0.008,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S167>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S167>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr], &rlim, 0.008, &fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S167>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S167>/Signal Conditioning Block2'
      // ForEachSliceSelector generated from: '<S167>/nedAccel_mps2' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   Inport: '<Root>/stateEstimate'

      fcsM_SignalConditioningBlock1_f
        (fcsModel_U.stateEstimate.nedAccel_mps2[ForEach_itr],
         &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelSignalConditioningParamsArray
         [ForEach_itr], &rtb_DifferenceInputs2_jb, 0.008,
         &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock2);

      // End of Outputs for SubSystem: '<S167>/Signal Conditioning Block2'

      // Product: '<S167>/Product' incorporates:
      //   ForEachSliceSelector generated from: '<S167>/accelFbGain'
      //   Inport: '<Root>/ctrlParams'

      rtb_DifferenceInputs2_jb *=
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelFbGainsArray
        [ForEach_itr];

      // Product: '<S167>/Product1' incorporates:
      //   BusAssignment: '<S108>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S167>/ctrlInputs'
      //   ForEachSliceSelector generated from: '<S167>/ffGain'
      //   Inport: '<Root>/ctrlParams'

      rtb_frcCmd_N =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ffGainsArray[ForEach_itr]
        * rtb_VectorConcatenate[ForEach_itr].cmd;

      // Outputs for Atomic SubSystem: '<S167>/pidWithDebug'
      // Sum: '<S167>/Sum' incorporates:
      //   BusAssignment: '<S108>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S167>/ctrlInputs'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S167>/Unit Delay'

      fcsModel_pidWithDebug_j(rtb_frcCmd_N - rtb_DifferenceInputs2_jb,
        rtb_DifferenceInputs2_i, rlim, rtb_VectorConcatenate[ForEach_itr].
        integratorReset, rtb_VectorConcatenate1[ForEach_itr],
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ctrlParamsArray[
        ForEach_itr], fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE,
        &rtb_frcCmd_N, &rtb_BusCreator_em, 0.008,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].pidWithDebug);

      // End of Outputs for SubSystem: '<S167>/pidWithDebug'

      // Update for UnitDelay: '<S167>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S167>/velCtrlOut '
      rtb_ImpAsg_InsertedFor_velCtrlO[ForEach_itr] = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S167>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr] = rtb_BusCreator_em;

      // ForEachSliceAssignment generated from: '<S167>/velCtrlFf'
      rtb_ImpAsg_InsertedFor_velCtrlF[ForEach_itr] = rtb_DifferenceInputs2_jb;

      // ForEachSliceAssignment generated from: '<S167>/filtMeas'
      rtb_ImpAsg_InsertedFor_filtMeas[ForEach_itr] = rlim;

      // ForEachSliceAssignment generated from: '<S167>/filtCmd'
      rtb_ImpAsg_InsertedFor_filtCmd_[ForEach_itr] = rtb_DifferenceInputs2_i;
    }

    // End of Outputs for SubSystem: '<S110>/For Each Subsystem'

    // DiscreteTransferFcn: '<S179>/Discrete Transfer Fcn' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   Inport: '<Root>/stateEstimate'
    //   Math: '<S179>/Transpose'
    //   Math: '<S179>/Transpose1'

    fcsModel_DW.DiscreteTransferFcn_tmp = fcsModel_U.stateEstimate.attitude_rad
      [2] -
      fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.firstOrderHeadingFilterDen
      [1] * fcsModel_DW.DiscreteTransferFcn_states;
    vyCmd_unitRange =
      fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.firstOrderHeadingFilterNum
      [0] * fcsModel_DW.DiscreteTransferFcn_tmp +
      fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.firstOrderHeadingFilterNum
      [1] * fcsModel_DW.DiscreteTransferFcn_states;

    // SignalConversion generated from: '<S180>/ SFunction ' incorporates:
    //   Constant: '<S179>/g'
    //   MATLAB Function: '<S179>/NE Accel Cmds To Roll Pitch Cmds'
    //   Sum: '<S179>/Sum'

    rtb_VectorConcatenate1[0] = rtb_ImpAsg_InsertedFor_velCtrlO[0];
    rtb_VectorConcatenate1[1] = rtb_ImpAsg_InsertedFor_velCtrlO[1];
    rtb_VectorConcatenate1[2] = rtb_ImpAsg_InsertedFor_velCtrlO[2] + -9.806;

    // MATLAB Function: '<S179>/NE Accel Cmds To Roll Pitch Cmds' incorporates:
    //   Constant: '<S179>/g'
    //   SignalConversion generated from: '<S180>/ SFunction '
    //   Sum: '<S179>/Sum'

    // MATLAB Function 'Velocity Controller/Assemble Inner Loop Inputs/nedAccelToRollPitchCmd/kinematicInversion/NE Accel Cmds To Roll Pitch Cmds': '<S180>:1' 
    // '<S180>:1:3' [rollCmd_rad,pitchCmd_rad] = ...
    // '<S180>:1:4'     accelsToDesiredRollPitchAngle_function(accelCmds_mps2, heading_rad); 
    // ACCELSTODESIREDROLLPITCHANGLE_FUNCTION converts NED accel commands to
    // desired roll and pitch angle
    // 'accelsToDesiredRollPitchAngle_function:5' att_eps = 1e-5;
    // 'accelsToDesiredRollPitchAngle_function:6' axOverAz = accelCmds_mps2(1)/accelCmds_mps2(3); 
    rtb_DifferenceInputs2_jb = rtb_ImpAsg_InsertedFor_velCtrlO[0] /
      (rtb_ImpAsg_InsertedFor_velCtrlO[2] + -9.806);

    // 'accelsToDesiredRollPitchAngle_function:7' ayOverAz = accelCmds_mps2(2)/accelCmds_mps2(3); 
    rtb_frcCmd_N = rtb_ImpAsg_InsertedFor_velCtrlO[1] /
      (rtb_ImpAsg_InsertedFor_velCtrlO[2] + -9.806);

    // 'accelsToDesiredRollPitchAngle_function:9' cos_heading = cos(heading_rad); 
    rtb_DifferenceInputs2_i = std::cos(vyCmd_unitRange);

    // 'accelsToDesiredRollPitchAngle_function:10' sin_heading = sin(heading_rad); 
    rlim = std::sin(vyCmd_unitRange);

    // 'accelsToDesiredRollPitchAngle_function:11' pitchCmd_rad = atan(cos_heading*axOverAz + sin_heading*ayOverAz); 
    vyCmd_unitRange = std::atan(rtb_DifferenceInputs2_i *
      rtb_DifferenceInputs2_jb + rlim * rtb_frcCmd_N);

    // 'accelsToDesiredRollPitchAngle_function:13' rollCmd_rad = atan(cos(pitchCmd_rad)*( sin_heading*axOverAz - ... 
    // 'accelsToDesiredRollPitchAngle_function:14'                    cos_heading*ayOverAz ) ); 
    rtb_DifferenceInputs2_jb = std::atan((rlim * rtb_DifferenceInputs2_jb -
      rtb_DifferenceInputs2_i * rtb_frcCmd_N) * std::cos(vyCmd_unitRange));

    // 'accelsToDesiredRollPitchAngle_function:15' if(abs(rollCmd_rad) <= att_eps) 
    if (std::abs(rtb_DifferenceInputs2_jb) <= 1.0E-5) {
      // 'accelsToDesiredRollPitchAngle_function:16' rollCmd_rad = 0;
      rtb_DifferenceInputs2_jb = 0.0;
    }

    // 'accelsToDesiredRollPitchAngle_function:19' if(abs(pitchCmd_rad) <= att_eps) 
    if (std::abs(vyCmd_unitRange) <= 1.0E-5) {
      // 'accelsToDesiredRollPitchAngle_function:20' pitchCmd_rad = 0;
      vyCmd_unitRange = 0.0;
    }

    // BusAssignment: '<S166>/Bus Assignment1' incorporates:
    //   Inport: '<Root>/stateEstimate'
    //   MATLAB Function: '<S179>/NE Accel Cmds To Roll Pitch Cmds'

    std::memset(&rtb_BusAssignment1_a, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment1_a.cmd = rtb_DifferenceInputs2_jb;
    rtb_BusAssignment1_a.meas = fcsModel_U.stateEstimate.attitude_rad[0];
    rtb_BusAssignment1_a.integratorReset = resetIntegrator;

    // BusAssignment: '<S166>/Bus Assignment2' incorporates:
    //   Inport: '<Root>/stateEstimate'
    //   MATLAB Function: '<S179>/NE Accel Cmds To Roll Pitch Cmds'

    std::memset(&rtb_BusAssignment2, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment2.cmd = vyCmd_unitRange;
    rtb_BusAssignment2.meas = fcsModel_U.stateEstimate.attitude_rad[1];
    rtb_BusAssignment2.integratorReset = resetIntegrator;

    // Logic: '<S166>/Logical Operator' incorporates:
    //   Constant: '<S169>/Constant'
    //   Constant: '<S170>/Constant'
    //   Logic: '<S166>/Logical Operator1'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S169>/Compare'
    //   RelationalOperator: '<S170>/Compare'

    rtb_Compare_d = (resetIntegrator || ((flightMode != enumFlightMode::
      ALT_CONTROL) && (flightMode != enumFlightMode::POS_CONTROL)));

    // Outputs for Atomic SubSystem: '<S166>/holdOutputAtCenter'
    // MATLAB Function: '<S172>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S177>:1'
    // '<S177>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
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

    // End of Outputs for SubSystem: '<S166>/holdOutputAtCenter'

    // BusAssignment: '<S166>/Bus Assignment4'
    // 'holdOutputAtCenter_function:17' output = last_input;
    std::memset(&rtb_BusAssignment4, 0, sizeof(busCtrlInputs));

    // Switch: '<S166>/Switch' incorporates:
    //   RelationalOperator: '<S168>/Compare'

    if (rtb_Compare_mi) {
      // BusAssignment: '<S166>/Bus Assignment4' incorporates:
      //   Constant: '<S166>/Constant5'

      rtb_BusAssignment4.feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S166>/Bus Assignment4'
      rtb_BusAssignment4.feedForwardCmd = fcsModel_DW.rcOutCmds.yawStick;
    }

    // End of Switch: '<S166>/Switch'

    // Outputs for Atomic SubSystem: '<S166>/holdOutputAtCenter'
    // BusAssignment: '<S166>/Bus Assignment4' incorporates:
    //   Inport: '<Root>/stateEstimate'
    //   MATLAB Function: '<S172>/holdOutputAtCenter'

    rtb_BusAssignment4.cmd = fcsModel_DW.last_input;

    // End of Outputs for SubSystem: '<S166>/holdOutputAtCenter'
    rtb_BusAssignment4.meas = fcsModel_U.stateEstimate.attitude_rad[2];
    rtb_BusAssignment4.integratorReset = rtb_Compare_d;

    // BusAssignment: '<S166>/Bus Assignment3' incorporates:
    //   Inport: '<Root>/stateEstimate'

    vyCmd_unitRange = fcsModel_U.stateEstimate.nedAccel_mps2[2];

    // Outputs for Atomic SubSystem: '<S182>/Signal Conditioning Block1'
    // MATLAB Function: '<S201>/Compute Natural Frequency' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMode_ComputeNaturalFrequency
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.filterBandwidth_radps,
       fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
       &rtb_DifferenceInputs2_jb);

    // MATLAB Function: '<S201>/Compute Numerator And Denominator' incorporates:
    //   Inport: '<Root>/ctrlParams'

    ComputeNumeratorAndDenominator(rtb_DifferenceInputs2_jb,
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
      &rtb_VectorConcatenate1[0], &rtb_MatrixMultiply[0],
      &rtb_ImpAsg_InsertedFor_neVelCmd[0], 0.008);

    // DiscreteTransferFcn: '<S201>/Discrete Transfer Fcn' incorporates:
    //   BusAssignment: '<S166>/Bus Assignment3'
    //   Inport: '<Root>/stateEstimate'

    rtb_DifferenceInputs2_i = (fcsModel_U.stateEstimate.nedAccel_mps2[2] -
      fcsModel_DW.DiscreteTransferFcn_states_n[0] *
      rtb_ImpAsg_InsertedFor_neVelCmd[1]) -
      fcsModel_DW.DiscreteTransferFcn_states_n[1] *
      rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_DiscreteTransferFcn = (rtb_VectorConcatenate1[0] *
      rtb_DifferenceInputs2_i + fcsModel_DW.DiscreteTransferFcn_states_n[0] *
      rtb_VectorConcatenate1[1]) + fcsModel_DW.DiscreteTransferFcn_states_n[1] *
      rtb_VectorConcatenate1[2];

    // Switch: '<S207>/Switch2' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   RelationalOperator: '<S207>/LowerRelop1'
    //   RelationalOperator: '<S207>/UpperRelop'
    //   Switch: '<S207>/Switch'

    if (rtb_DiscreteTransferFcn >
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [1]) {
      rtb_DiscreteTransferFcn =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [1];
    } else if (rtb_DiscreteTransferFcn <
               fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
               [0]) {
      // Switch: '<S207>/Switch'
      rtb_DiscreteTransferFcn =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [0];
    }

    // End of Switch: '<S207>/Switch2'

    // Sum: '<S204>/Difference Inputs1' incorporates:
    //   UnitDelay: '<S204>/Delay Input2'
    //
    //  Block description for '<S204>/Difference Inputs1':
    //
    //   Add in CPU
    //
    //  Block description for '<S204>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DifferenceInputs2_jb = rtb_DiscreteTransferFcn -
      fcsModel_DW.DelayInput2_DSTATE;

    // Product: '<S204>/delta rise limit' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   SampleTimeMath: '<S204>/sample time'
    //
    //  About '<S204>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_frcCmd_N =
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputAccelLimits
      [1] * 0.008;

    // Switch: '<S214>/Switch2' incorporates:
    //   RelationalOperator: '<S214>/LowerRelop1'

    if (rtb_DifferenceInputs2_jb <= rtb_frcCmd_N) {
      // Product: '<S204>/delta fall limit' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   SampleTimeMath: '<S204>/sample time'
      //
      //  About '<S204>/sample time':
      //   y = K where K = ( w * Ts )

      rtb_frcCmd_N =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputAccelLimits
        [0] * 0.008;

      // Switch: '<S214>/Switch' incorporates:
      //   RelationalOperator: '<S214>/UpperRelop'

      if (rtb_DifferenceInputs2_jb >= rtb_frcCmd_N) {
        rtb_frcCmd_N = rtb_DifferenceInputs2_jb;
      }

      // End of Switch: '<S214>/Switch'
    }

    // End of Switch: '<S214>/Switch2'

    // Sum: '<S204>/Difference Inputs2' incorporates:
    //   UnitDelay: '<S204>/Delay Input2'
    //
    //  Block description for '<S204>/Difference Inputs2':
    //
    //   Add in CPU
    //
    //  Block description for '<S204>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DiscreteTransferFcn = rtb_frcCmd_N + fcsModel_DW.DelayInput2_DSTATE;

    // MATLAB Function: '<S202>/Compute Natural Frequency' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMode_ComputeNaturalFrequency
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.filterBandwidth_radps,
       fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
       &rtb_DifferenceInputs2_jb);

    // MATLAB Function: '<S202>/Compute Filter Numerator And Denominator' incorporates:
    //   Inport: '<Root>/ctrlParams'

    ComputeFilterNumeratorAndDenomi(rtb_DifferenceInputs2_jb,
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filterParams.dampingRatio_nd,
      &rtb_VectorConcatenate1[0], &rtb_MatrixMultiply[0], 0.008);

    // DiscreteTransferFcn: '<S202>/Discrete Transfer Fcn' incorporates:
    //   BusAssignment: '<S166>/Bus Assignment3'
    //   Inport: '<Root>/stateEstimate'

    rlim = (fcsModel_U.stateEstimate.nedAccel_mps2[2] -
            fcsModel_DW.DiscreteTransferFcn_states_nh[0] * rtb_MatrixMultiply[1])
      - fcsModel_DW.DiscreteTransferFcn_states_nh[1] * rtb_MatrixMultiply[2];
    rtb_DiscreteTransferFcn_k = (rtb_VectorConcatenate1[0] * rlim +
      fcsModel_DW.DiscreteTransferFcn_states_nh[0] * rtb_VectorConcatenate1[1])
      + fcsModel_DW.DiscreteTransferFcn_states_nh[1] * rtb_VectorConcatenate1[2];

    // Switch: '<S206>/Switch2' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   RelationalOperator: '<S206>/LowerRelop1'
    //   RelationalOperator: '<S206>/UpperRelop'
    //   Switch: '<S206>/Switch'

    if (rtb_DiscreteTransferFcn_k >
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [1]) {
      rtb_DiscreteTransferFcn_k =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [1];
    } else if (rtb_DiscreteTransferFcn_k <
               fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
               [0]) {
      // Switch: '<S206>/Switch'
      rtb_DiscreteTransferFcn_k =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputLimits
        [0];
    }

    // End of Switch: '<S206>/Switch2'

    // Sum: '<S203>/Difference Inputs1' incorporates:
    //   UnitDelay: '<S203>/Delay Input2'
    //
    //  Block description for '<S203>/Difference Inputs1':
    //
    //   Add in CPU
    //
    //  Block description for '<S203>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DifferenceInputs2_jb = rtb_DiscreteTransferFcn_k -
      fcsModel_DW.DelayInput2_DSTATE_e;

    // Product: '<S203>/delta rise limit' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   SampleTimeMath: '<S203>/sample time'
    //
    //  About '<S203>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_frcCmd_N =
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
      [1] * 0.008;

    // Switch: '<S213>/Switch2' incorporates:
    //   RelationalOperator: '<S213>/LowerRelop1'

    if (rtb_DifferenceInputs2_jb <= rtb_frcCmd_N) {
      // Product: '<S203>/delta fall limit' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   SampleTimeMath: '<S203>/sample time'
      //
      //  About '<S203>/sample time':
      //   y = K where K = ( w * Ts )

      rtb_frcCmd_N =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.measSignalConditioningParams.filteredInputRateLimits
        [0] * 0.008;

      // Switch: '<S213>/Switch' incorporates:
      //   RelationalOperator: '<S213>/UpperRelop'

      if (rtb_DifferenceInputs2_jb >= rtb_frcCmd_N) {
        rtb_frcCmd_N = rtb_DifferenceInputs2_jb;
      }

      // End of Switch: '<S213>/Switch'
    }

    // End of Switch: '<S213>/Switch2'

    // Sum: '<S203>/Difference Inputs2' incorporates:
    //   UnitDelay: '<S203>/Delay Input2'
    //
    //  Block description for '<S203>/Difference Inputs2':
    //
    //   Add in CPU
    //
    //  Block description for '<S203>/Delay Input2':
    //
    //   Store in Global RAM

    rtb_DiscreteTransferFcn_k = rtb_frcCmd_N + fcsModel_DW.DelayInput2_DSTATE_e;

    // Update for DiscreteTransferFcn: '<S201>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states_n[1] =
      fcsModel_DW.DiscreteTransferFcn_states_n[0];
    fcsModel_DW.DiscreteTransferFcn_states_n[0] = rtb_DifferenceInputs2_i;

    // Update for UnitDelay: '<S204>/Delay Input2'
    //
    //  Block description for '<S204>/Delay Input2':
    //
    //   Store in Global RAM

    fcsModel_DW.DelayInput2_DSTATE = rtb_DiscreteTransferFcn;

    // Update for DiscreteTransferFcn: '<S202>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states_nh[1] =
      fcsModel_DW.DiscreteTransferFcn_states_nh[0];
    fcsModel_DW.DiscreteTransferFcn_states_nh[0] = rlim;

    // Update for UnitDelay: '<S203>/Delay Input2'
    //
    //  Block description for '<S203>/Delay Input2':
    //
    //   Store in Global RAM

    fcsModel_DW.DelayInput2_DSTATE_e = rtb_DiscreteTransferFcn_k;

    // End of Outputs for SubSystem: '<S182>/Signal Conditioning Block1'

    // Product: '<S182>/Product' incorporates:
    //   Gain: '<S182>/Gain'
    //   Inport: '<Root>/ctrlParams'

    rtb_DifferenceInputs2_i = -rtb_DiscreteTransferFcn *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ctrlParams.Kd;

    // Sum: '<S182>/Sum' incorporates:
    //   BusAssignment: '<S166>/Bus Assignment3'
    //   Inport: '<Root>/ctrlParams'
    //   Product: '<S182>/Product1'

    rtb_frcCmd_N = rtb_ImpAsg_InsertedFor_velCtrlO[2] *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ffGain +
      rtb_DifferenceInputs2_i;

    // Outputs for Atomic SubSystem: '<S182>/Signal Conditioning Block'
    // BusAssignment: '<S166>/Bus Assignment3' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   UnitDelay: '<S182>/Unit Delay'

    fcsM_SignalConditioningBlock1_f(rtb_ImpAsg_InsertedFor_velCtrlO[2],
      &fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.cmdSignalConditioningParams,
      &rtb_DifferenceInputs2_i, 0.008, &fcsModel_DW.SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S182>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S182>/pidWithDebug'
    fcsModel_pidWithDebug_j(rtb_frcCmd_N, rtb_DifferenceInputs2_i,
      rtb_DiscreteTransferFcn_k, rtb_Compare_d, 0.0,
      &fcsModel_U.ctrlParams.outerLoopCtrlParams.zAccelCtrlParams.ctrlParams,
      fcsModel_DW.UnitDelay_DSTATE, &rtb_DifferenceInputs2_jb,
      &rtb_BusCreator_em, 0.008, &fcsModel_DW.pidWithDebug);

    // End of Outputs for SubSystem: '<S182>/pidWithDebug'

    // Product: '<S166>/Divide1' incorporates:
    //   Constant: '<S178>/g'
    //   Inport: '<Root>/ctrlParams'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S166>/Product5'
    //   Product: '<S178>/Divide3'
    //   Sum: '<S166>/Sum2'
    //   Trigonometry: '<S166>/Sin2'
    //   Trigonometry: '<S166>/Sin3'

    rtb_frcCmd_N = 1.0 / (std::cos(fcsModel_U.stateEstimate.attitude_rad[0]) *
                          std::cos(fcsModel_U.stateEstimate.attitude_rad[1])) *
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.baseMass_kg *
       -9.806 + rtb_DifferenceInputs2_jb);

    // RelationalOperator: '<S105>/Compare' incorporates:
    //   Constant: '<S105>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_d = (flightMode == enumFlightMode::POS_CONTROL);

    // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
    //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
    //   Constant: '<S3>/Constant'
    //   Inport: '<Root>/stateEstimate'

    std::memcpy(&rtb_VectorConcatenate[0],
                &fcsModel_ConstP.pooled3.attCtrlInputs.ctrlInputsArray[0], 3U *
                sizeof(busCtrlInputs));

    // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S107>:1' 
    // '<S107>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
    // '<S107>:1:3' outDebug = throttleCmd_N;
    rtb_DiscreteTransferFcn = fcsModel_DW.rcOutCmds.throttleStick;

    //  This is a stop gap setup where we are only assuming that rate control
    //  is active and therefore not setting up attCtrlInputs for Euler angle
    //  control
    // '<S107>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
    rtb_VectorConcatenate[0].cmd = fcsModel_DW.rcOutCmds.rollStick;

    // '<S107>:1:8' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.attitude_rad[0];

    // '<S107>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
    rtb_VectorConcatenate[1].cmd = fcsModel_DW.rcOutCmds.pitchStick;

    // '<S107>:1:10' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.attitude_rad[1];

    // '<S107>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
    rtb_VectorConcatenate[2].cmd = fcsModel_DW.rcOutCmds.yawStick;

    // '<S107>:1:12' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
    rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.attitude_rad[2];

    // '<S107>:1:14' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[0].integratorReset = resetIntegrator;

    // '<S107>:1:15' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[1].integratorReset = resetIntegrator;

    // '<S107>:1:16' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = true; 
    rtb_VectorConcatenate[2].integratorReset = true;

    // RelationalOperator: '<S104>/Compare' incorporates:
    //   Constant: '<S104>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_mi = (flightMode != enumFlightMode::ALT_CONTROL);

    // Switch: '<S3>/Switch2' incorporates:
    //   Switch: '<S3>/Switch'

    if (rtb_Compare_d) {
      // Switch: '<S3>/Switch2' incorporates:
      //   BusAssignment: '<S166>/Bus Assignment'
      //   Concatenate: '<S166>/Vector Concatenate'

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
      //   BusAssignment: '<S166>/Bus Assignment'
      //   Concatenate: '<S166>/Vector Concatenate'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Switch: '<S3>/Switch'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      std::memcpy(&fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0],
                  &rtb_VectorConcatenate[0], sizeof(busCtrlInputs) << 1U);
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment4;
    }

    // End of Switch: '<S3>/Switch2'
  }

  // Outputs for Iterator SubSystem: '<S20>/Attitude Control' incorporates:
  //   ForEach: '<S62>/For Each'

  for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
    // Outputs for Atomic SubSystem: '<S62>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S62>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].cmd,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.cmdSignalConditioningParamsArray
       [ForEach_itr_l], &rlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S62>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S62>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].meas,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.measSignalConditioningParamsArray
       [ForEach_itr_l], &plim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S62>/Signal Conditioning Block1'

    // MATLAB Function: '<S62>/pickAttitudeCmdAndMeas' incorporates:
    //   Constant: '<S20>/Constant'
    //   ForEachSliceSelector generated from: '<S62>/index'

    vxCmd_unitRange = rlim;
    yCmd = plim;

    //  Passes cmd and meas as it is for roll and pitch channel
    //  but for yaw channel computes shortest angular distance between cmd Yaw
    //  and meas Yaw and overwrites Yaw cmd with that error and sets the meas Yaw to 
    //  zero for PID block
    // MATLAB Function 'Attitude Controller/Attitude Control/pickAttitudeCmdAndMeas': '<S67>:1' 
    // '<S67>:1:6' if index == cast(3, 'uint8')
    if (fcsModel_ConstP.Constant_Value_e[ForEach_itr_l] == 3) {
      // '<S67>:1:7' diff = mod(( cmd - meas + pi ), 2*pi) - pi;
      vxCmd_unitRange = (rlim - plim) + 3.1415926535897931;
      if (vxCmd_unitRange == 0.0) {
        vyCmd_unitRange = 0.0;
      } else {
        vyCmd_unitRange = std::fmod(vxCmd_unitRange, 6.2831853071795862);
        resetIntegrator = (vyCmd_unitRange == 0.0);
        if (!resetIntegrator) {
          yCmd = std::abs(vxCmd_unitRange / 6.2831853071795862);
          resetIntegrator = (std::abs(yCmd - std::floor(yCmd + 0.5)) <=
                             2.2204460492503131E-16 * yCmd);
        }

        if (resetIntegrator) {
          vyCmd_unitRange = 0.0;
        } else if (vxCmd_unitRange < 0.0) {
          vyCmd_unitRange += 6.2831853071795862;
        }
      }

      vxCmd_unitRange = vyCmd_unitRange - 3.1415926535897931;

      // '<S67>:1:8' if diff < -pi
      if (vyCmd_unitRange - 3.1415926535897931 < -3.1415926535897931) {
        // '<S67>:1:9' diff = diff + 2*pi;
        vxCmd_unitRange = (vyCmd_unitRange - 3.1415926535897931) +
          6.2831853071795862;
      }

      // '<S67>:1:12' cmd = diff;
      // '<S67>:1:13' meas = 0;
      yCmd = 0.0;
    }

    // Outputs for Atomic SubSystem: '<S62>/pidWithDebug'
    // ForEachSliceSelector generated from: '<S62>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   MATLAB Function: '<S62>/pickAttitudeCmdAndMeas'
    //   UnitDelay: '<S62>/Unit Delay'

    fcsModel_pidWithDebug
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       feedForwardCmd, vxCmd_unitRange, yCmd,
       fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       integratorReset, 0.0,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.ctrlParamsArray[ForEach_itr_l],
       fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE,
       &vyCmd_unitRange, &rtb_BusCreator_f, 0.004,
       &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].pidWithDebug);

    // End of Outputs for SubSystem: '<S62>/pidWithDebug'

    // Update for UnitDelay: '<S62>/Unit Delay'
    fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE = vyCmd_unitRange;

    // ForEachSliceAssignment generated from: '<S62>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[ForEach_itr_l] =
      rtb_BusCreator_f;

    // ForEachSliceAssignment generated from: '<S62>/angRateCmd '
    rtb_ImpAsg_InsertedFor_angRateC[ForEach_itr_l] = vyCmd_unitRange;

    // ForEachSliceAssignment generated from: '<S62>/measFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[ForEach_itr_l] =
      plim;

    // ForEachSliceAssignment generated from: '<S62>/cmdFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[ForEach_itr_l] =
      rlim;
  }

  // End of Outputs for SubSystem: '<S20>/Attitude Control'

  // Switch: '<S19>/Switch'
  rtb_ImpAsg_InsertedFor_velCtrlO[0] = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_ImpAsg_InsertedFor_velCtrlO[1] = rtb_ImpAsg_InsertedFor_angRateC[1];

  // Switch: '<S20>/Switch' incorporates:
  //   Constant: '<S63>/Constant'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   RelationalOperator: '<S63>/Compare'
  //   Switch: '<S19>/Switch'

  if (flightMode == enumFlightMode::STABILIZE) {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] =
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2].cmd;
  } else {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] = rtb_ImpAsg_InsertedFor_angRateC[2];
  }

  // End of Switch: '<S20>/Switch'

  // MATLAB Function: '<S19>/EulerRates2BodyRates' incorporates:
  //   Inport: '<Root>/stateEstimate'
  //   Switch: '<S19>/Switch'

  // MATLAB Function 'EulerRates2BodyRates': '<S61>:1'
  // '<S61>:1:3' bodyRates_radps = eulerRates2bodyRates_function(taitBryanRates_radps,shipOrientation_rad); 
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
  vyCmd_unitRange = fcsModel_U.stateEstimate.attitude_rad[1];

  // 'eulerRates2bodyRates_function:20' eps = 10^(-12);
  // 'eulerRates2bodyRates_function:21' limit = pi/740;
  // Check for pm pi/2 rotation to avoid NaNs
  // 'eulerRates2bodyRates_function:24' if( abs( abs(pitch)- pi/2 ) <= limit || abs( abs(pitch) - 3*pi/2 ) <= limit) 
  rlim = std::abs(fcsModel_U.stateEstimate.attitude_rad[1]);
  if ((std::abs(rlim - 1.5707963267948966) <= 0.004245395477824045) || (std::abs
       (rlim - 4.71238898038469) <= 0.004245395477824045)) {
    // 'eulerRates2bodyRates_function:25' if((abs(pitch)- pi/2) <= 0 || (abs(pitch) - 3*pi/2) <= 0) 
    if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) - 1.5707963267948966 <=
        0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      vyCmd_unitRange = (rlim - 0.004245395477824045) * plim;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      vyCmd_unitRange = (rlim - 0.004245395477824045) * plim;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        plim = -1.0;
      } else {
        plim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      vyCmd_unitRange = (rlim + 0.004245395477824045) * plim;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  rlim = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  plim = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  vxCmd_unitRange = std::cos(vyCmd_unitRange);
  rtb_Transpose[0] = 1.0;
  rtb_Transpose[3] = 0.0;
  rtb_Transpose[6] = -std::sin(vyCmd_unitRange);
  rtb_Transpose[1] = 0.0;
  rtb_Transpose[4] = plim;
  rtb_Transpose[7] = rlim * vxCmd_unitRange;
  rtb_Transpose[2] = 0.0;
  rtb_Transpose[5] = -rlim;
  rtb_Transpose[8] = plim * vxCmd_unitRange;

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

  // End of MATLAB Function: '<S19>/EulerRates2BodyRates'

  // BusCreator: '<S19>/Bus Creator' incorporates:
  //   Concatenate: '<S19>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_VectorConcatenate[0].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S19>/Bus Creator3' incorporates:
  //   Concatenate: '<S19>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_VectorConcatenate[1].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S19>/Bus Creator4' incorporates:
  //   Concatenate: '<S19>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[2].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_angRateC[2];
  rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_VectorConcatenate[2].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[2].trackingCtrlCmd = 0.0;

  // Outputs for Atomic SubSystem: '<S2>/Angular Rate Controller'
  // Outputs for Iterator SubSystem: '<S18>/For Each Subsystem' incorporates:
  //   ForEach: '<S21>/For Each'

  for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
    // Outputs for Atomic SubSystem: '<S21>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S21>/ctrlInputs' incorporates:
    //   BusCreator: '<S19>/Bus Creator1'
    //   Concatenate: '<S19>/Vector Concatenate'
    //   Inport: '<Root>/ctrlParams'
    //   UnitDelay: '<S21>/Unit Delay'

    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].cmd,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.cmdSignalConditioningParamsArray
      [ForEach_itr_p], &vyCmd_unitRange, 0.004,
      &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S21>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S21>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1(rtb_VectorConcatenate[ForEach_itr_p].meas,
      &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.measSignalConditioningParamsArray
      [ForEach_itr_p], &rlim, 0.004, &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
      SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S21>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S21>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, vyCmd_unitRange, rlim,
                          rtb_VectorConcatenate[ForEach_itr_p].integratorReset,
                          0.0,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.ctrlParamsArray
                          [ForEach_itr_p],
                          fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
                          UnitDelay_DSTATE, &plim, &rtb_BusCreator_f, 0.004,
                          &fcsModel_DW.CoreSubsys_a[ForEach_itr_p].pidWithDebug);

    // End of Outputs for SubSystem: '<S21>/pidWithDebug'

    // Update for UnitDelay: '<S21>/Unit Delay'
    fcsModel_DW.CoreSubsys_a[ForEach_itr_p].UnitDelay_DSTATE = plim;

    // ForEachSliceAssignment generated from: '<S21>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[ForEach_itr_p]
      = rtb_BusCreator_f;

    // ForEachSliceAssignment generated from: '<S21>/angAccelCmd_radps2'
    rtb_ImpAsg_InsertedFor_angAccel[ForEach_itr_p] = plim;

    // ForEachSliceAssignment generated from: '<S21>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[ForEach_itr_p] =
      rlim;

    // ForEachSliceAssignment generated from: '<S21>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[ForEach_itr_p] =
      vyCmd_unitRange;
  }

  // End of Outputs for SubSystem: '<S18>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S21>/angAccelCmd_radps2'

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

  // RelationalOperator: '<S9>/Compare' incorporates:
  //   Constant: '<S9>/Constant'

  resetIntegrator = (rtb_chirpTrigger == enumChirpTrigger::ON);

  // MATLAB Function: '<S15>/Generate Chirp' incorporates:
  //   Inport: '<Root>/ctrlParams'
  //   UnitDelay: '<S15>/Unit Delay'
  //   UnitDelay: '<S15>/Unit Delay1'

  // MATLAB Function 'sysIdInputGeneration/chirpInjection/Generate Chirp': '<S17>:1' 
  // '<S17>:1:4' [chirpSignal, chirpTime_s, theta_rad] = generateChirp_function(trigger, amp, fStart_Hz, fEnd_Hz, tRec_s, ... 
  // '<S17>:1:5'                                        fadeInTime_s, fadeOutTime_s, c1, c2, inElpsdChirpTime_s, ... 
  // '<S17>:1:6'                                        inTheta_rad, sampleTime_s); 
  // GENERATECHIRP_FUNCTION generates chirp for sysId purposes
  //
  // Inputs:
  // trigger: Binary signal, if it is 0 this function outputs 0
  // amp: Chirp amplitude
  // fStart_Hz: Chirp starting frequency
  // fEnd_Hz: Chirp end frequency
  // tRec_s: Total chirp duratin
  // fadeInTime_s: Duration at the beginning during which amplitude will be increase to 1 
  // fadeOutTime_s: Duration at the end during which amplitude will decrease to 0 
  // c1: Chirp parameter 1
  // c2: Chirp parameter 2
  // inElpsdDhirpTime_s: Previous value of elapsed time since the chirp was triggered else it is set to 0 
  // inTheta_rad: Previous chirp phase angle
  // sampleTime_s: Sample time
  //
  // Outputs:
  // chirpSignal: Chirp signal at the current time step
  // chirpTime_s: Elapsed time since chirp was triggerd else it is set to 0
  // theta_rad: Chirp phase angle
  //  trim time
  // 'generateChirp_function:26' if(trigger == 0)
  if (!resetIntegrator) {
    // 'generateChirp_function:27' chirpSignal = 0;
    vyCmd_unitRange = 0.0;

    // 'generateChirp_function:28' chirpTime_s = 0;
    plim = 0.0;

    // 'generateChirp_function:29' theta_rad = 0;
    rlim = 0.0;
  } else {
    //  Chirp is triggered
    //  Convert frequencies to rad/s
    // 'generateChirp_function:36' fStart_radps = 2*pi*fStart_Hz;
    rlim = 6.2831853071795862 *
      fcsModel_U.ctrlParams.sysIdInjectionParams.fStart_Hz;

    // 'generateChirp_function:37' fEnd_radps = 2*pi*fEnd_Hz;
    //  Find the slope of amplitude during the fade in and fade out time
    // 'generateChirp_function:40' fadeInSlope = 1/fadeInTime_s;
    // 'generateChirp_function:41' fadeOutSlope = 1/fadeOutTime_s;
    //  Find the time when the fade out starts
    // 'generateChirp_function:44' tFadeOutStart_s = tRec_s - fadeOutTime_s;
    vyCmd_unitRange = fcsModel_U.ctrlParams.sysIdInjectionParams.tRec_s -
      fcsModel_U.ctrlParams.sysIdInjectionParams.fadeOutTime_s;

    //  Time required for 2 full cycles at starting frequency
    // 'generateChirp_function:47' timeLowCycles_s = 2*(2*pi/fStart_radps);
    plim = 6.2831853071795862 / rlim * 2.0;

    //  Do fade in and fade out at appropriate times
    // 'generateChirp_function:50' if(inElpsdChirpTime_s <= fadeInTime_s)
    if (fcsModel_DW.UnitDelay1_DSTATE <=
        fcsModel_U.ctrlParams.sysIdInjectionParams.fadeInTime_s) {
      // 'generateChirp_function:51' a = fadeInSlope*inElpsdChirpTime_s;
      vyCmd_unitRange = 1.0 /
        fcsModel_U.ctrlParams.sysIdInjectionParams.fadeInTime_s *
        fcsModel_DW.UnitDelay1_DSTATE;
    } else if (fcsModel_DW.UnitDelay1_DSTATE >= vyCmd_unitRange) {
      // 'generateChirp_function:52' elseif(inElpsdChirpTime_s >= tFadeOutStart_s) 
      // 'generateChirp_function:53' a = 1 - fadeOutSlope*(inElpsdChirpTime_s - tFadeOutStart_s); 
      vyCmd_unitRange = 1.0 - 1.0 /
        fcsModel_U.ctrlParams.sysIdInjectionParams.fadeOutTime_s *
        (fcsModel_DW.UnitDelay1_DSTATE - vyCmd_unitRange);
    } else {
      // 'generateChirp_function:54' else
      // 'generateChirp_function:55' a = 1;
      vyCmd_unitRange = 1.0;
    }

    //  Create chirps
    // 'generateChirp_function:59' chirpSignal = amp*a*sin(inTheta_rad);
    vyCmd_unitRange *= std::sin(fcsModel_DW.UnitDelay_DSTATE_j);

    //  Increment theta
    // 'generateChirp_function:62' if(inElpsdChirpTime_s < timeLowCycles_s)
    if (fcsModel_DW.UnitDelay1_DSTATE < plim) {
      // 'generateChirp_function:63' theta_rad = wrapTo2Pi_(inTheta_rad + fStart_radps*sampleTime_s); 
      rlim = rlim * 0.004 + fcsModel_DW.UnitDelay_DSTATE_j;

      // WRAPTO2PI_ custom implementation of matlab's wrapTp2Pi
      // 'wrapTo2Pi_:3' twoPi = 2*pi;
      // 'wrapTo2Pi_:4' outAng = inAng -  twoPi* floor( inAng / twoPi );
      rlim -= std::floor(rlim / 6.2831853071795862) * 6.2831853071795862;
    } else {
      // 'generateChirp_function:64' else
      // 'generateChirp_function:65' k = c2*(exp(c1*(inElpsdChirpTime_s - timeLowCycles_s)/(tRec_s - timeLowCycles_s)) - 1); 
      // 'generateChirp_function:66' theta_rad = wrapTo2Pi_(inTheta_rad + (fStart_radps + k*(fEnd_radps - fStart_radps))*sampleTime_s); 
      rlim = ((std::exp((fcsModel_DW.UnitDelay1_DSTATE - plim) *
                        fcsModel_U.ctrlParams.sysIdInjectionParams.c1 /
                        (fcsModel_U.ctrlParams.sysIdInjectionParams.tRec_s -
                         plim)) - 1.0) *
              fcsModel_U.ctrlParams.sysIdInjectionParams.c2 *
              (6.2831853071795862 *
               fcsModel_U.ctrlParams.sysIdInjectionParams.fEnd_hz - rlim) + rlim)
        * 0.004 + fcsModel_DW.UnitDelay_DSTATE_j;

      // WRAPTO2PI_ custom implementation of matlab's wrapTp2Pi
      // 'wrapTo2Pi_:3' twoPi = 2*pi;
      // 'wrapTo2Pi_:4' outAng = inAng -  twoPi* floor( inAng / twoPi );
      rlim -= std::floor(rlim / 6.2831853071795862) * 6.2831853071795862;
    }

    //  Increment chirp time
    // 'generateChirp_function:70' chirpTime_s = inElpsdChirpTime_s + sampleTime_s; 
    plim = fcsModel_DW.UnitDelay1_DSTATE + 0.004;

    // 'generateChirp_function:72' if(inElpsdChirpTime_s > tRec_s)
    if (fcsModel_DW.UnitDelay1_DSTATE >
        fcsModel_U.ctrlParams.sysIdInjectionParams.tRec_s) {
      // 'generateChirp_function:73' chirpSignal = 0;
      vyCmd_unitRange = 0.0;

      // 'generateChirp_function:74' theta_rad = 0;
      rlim = 0.0;
    }
  }

  // End of MATLAB Function: '<S15>/Generate Chirp'

  // DiscreteTransferFcn: '<S15>/Discrete Transfer Fcn' incorporates:
  //   Gain: '<S16>/Output'
  //   Inport: '<Root>/ctrlParams'
  //   Math: '<S8>/Transpose1'
  //   RandomNumber: '<S16>/White Noise'

  vxCmd_unitRange = 0.5 * fcsModel_DW.NextOutput -
    fcsModel_U.ctrlParams.sysIdInjectionParams.filterDen[1] *
    fcsModel_DW.DiscreteTransferFcn_states_c;

  // Product: '<S15>/Product' incorporates:
  //   DiscreteTransferFcn: '<S15>/Discrete Transfer Fcn'
  //   Inport: '<Root>/ctrlParams'
  //   Math: '<S8>/Transpose'
  //   Sum: '<S15>/Sum'

  yCmd = resetIntegrator ?
    (fcsModel_U.ctrlParams.sysIdInjectionParams.filterNum[0] * vxCmd_unitRange +
     fcsModel_U.ctrlParams.sysIdInjectionParams.filterNum[1] *
     fcsModel_DW.DiscreteTransferFcn_states_c) + vyCmd_unitRange : 0.0;

  // SignalConversion generated from: '<S1>/Matrix Multiply' incorporates:
  //   BusCreator: '<S2>/Bus Creator1'
  //   Constant: '<S10>/Constant'
  //   Constant: '<S11>/Constant'
  //   Constant: '<S12>/Constant'
  //   Constant: '<S13>/Constant'
  //   Gain: '<S1>/Gain'
  //   Inport: '<Root>/ctrlParams'
  //   Product: '<S8>/Product'
  //   Product: '<S8>/Product1'
  //   Product: '<S8>/Product2'
  //   Product: '<S8>/Product3'
  //   RelationalOperator: '<S10>/Compare'
  //   RelationalOperator: '<S11>/Compare'
  //   RelationalOperator: '<S12>/Compare'
  //   RelationalOperator: '<S13>/Compare'
  //   Sum: '<S8>/Sum'
  //   Sum: '<S8>/Sum1'
  //   Sum: '<S8>/Sum2'
  //   Sum: '<S8>/Sum3'

  tmp = (rtb_chirpType == enumChirpType::FZ ?
         fcsModel_U.ctrlParams.sysIdInjectionParams.fzAmp * yCmd : 0.0) +
    fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N;
  tmp_0 = (rtb_chirpType == enumChirpType::MX ?
           fcsModel_U.ctrlParams.sysIdInjectionParams.mxAmp * yCmd : 0.0) +
    rtb_ImpAsg_InsertedFor_angRateC[0];
  tmp_1 = (rtb_chirpType == enumChirpType::MY ?
           fcsModel_U.ctrlParams.sysIdInjectionParams.myAmp * yCmd : 0.0) +
    rtb_ImpAsg_InsertedFor_angRateC[1];
  tmp_2 = (rtb_chirpType == enumChirpType::MY ?
           fcsModel_U.ctrlParams.sysIdInjectionParams.myAmp * yCmd : 0.0) +
    rtb_ImpAsg_InsertedFor_angRateC[2] * 0.7;

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);
  for (tCmd = 0; tCmd < 4; tCmd++) {
    real_T rtb_DiscreteTransferFcn_m;

    // Product: '<S1>/Matrix Multiply' incorporates:
    //   Constant: '<S1>/Constant'
    //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

    rtb_DiscreteTransferFcn_m = ((fcsModel_ConstP.Constant_Value_c[tCmd + 4] *
      tmp_0 + fcsModel_ConstP.Constant_Value_c[tCmd] * tmp) +
      fcsModel_ConstP.Constant_Value_c[tCmd + 8] * tmp_1) +
      fcsModel_ConstP.Constant_Value_c[tCmd + 12] * tmp_2;

    // Saturate: '<S1>/Saturation'
    if (rtb_DiscreteTransferFcn_m > 839601.76328711538) {
      rtb_DiscreteTransferFcn_m = 839601.76328711538;
    } else if (rtb_DiscreteTransferFcn_m < 0.0) {
      rtb_DiscreteTransferFcn_m = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    vyCmd_unitRange = 9.5492965855137211 * std::sqrt(rtb_DiscreteTransferFcn_m)
      - 0.11372544828835565 * fcsModel_DW.DiscreteTransferFcn_states_d[tCmd];
    rtb_DiscreteTransferFcn_m = 0.55686272414417781 * vyCmd_unitRange +
      0.55686272414417781 * fcsModel_DW.DiscreteTransferFcn_states_d[tCmd];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      rtb_DiscreteTransferFcn_m = -1.0;
    }

    // End of Switch: '<S1>/Switch'

    // Outport: '<Root>/actuatorsCmds'
    fcsModel_Y.actuatorsCmds[tCmd] = rtb_DiscreteTransferFcn_m;

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    DiscreteTransferFcn_tmp_b[tCmd] = vyCmd_unitRange;

    // Product: '<S1>/Matrix Multiply' incorporates:
    //   Constant: '<S1>/Constant'
    //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

    rtb_DiscreteTransferFcn_e[tCmd] = rtb_DiscreteTransferFcn_m;
  }

  // Outputs for Iterator SubSystem: '<S1>/For Each Subsystem' incorporates:
  //   ForEach: '<S7>/For Each'

  for (ForEach_itr_g = 0; ForEach_itr_g < 4; ForEach_itr_g++) {
    uint32_T rtb_Prelookup_o1;

    // ForEachSliceSelector generated from: '<S7>/propellerSpdCmds_rpm' incorporates:
    //   Switch: '<S1>/Switch'

    vyCmd_unitRange = rtb_DiscreteTransferFcn_e[ForEach_itr_g];

    // Saturate: '<S7>/Saturation'
    if (vyCmd_unitRange > 9325.0) {
      // PreLookup: '<S7>/Prelookup'
      vyCmd_unitRange = 9325.0;
    } else if (vyCmd_unitRange < 2250.0) {
      // PreLookup: '<S7>/Prelookup'
      vyCmd_unitRange = 2250.0;
    }

    // End of Saturate: '<S7>/Saturation'

    // PreLookup: '<S7>/Prelookup'
    rtb_Prelookup_o1 = plook_bincpag(vyCmd_unitRange,
      &fcsModel_ConstP.Prelookup_BreakpointsData[0], 14U, &vyCmd_unitRange,
      &fcsModel_DW.CoreSubsys[ForEach_itr_g].Prelookup_DWORK1);

    // ForEachSliceAssignment generated from: '<S7>/mtrPwmCmds' incorporates:
    //   Interpolation_n-D: '<S7>/Interpolation Using Prelookup'

    rtb_ImpAsg_InsertedFor_mtrPwmCm[ForEach_itr_g] = intrp1d_la(rtb_Prelookup_o1,
      vyCmd_unitRange, &fcsModel_ConstP.InterpolationUsingPrelookup_Tab[0], 14U);
  }

  // End of Outputs for SubSystem: '<S1>/For Each Subsystem'

  // Outport: '<Root>/actuatorsPwmCmds' incorporates:
  //   ForEachSliceAssignment generated from: '<S7>/mtrPwmCmds'

  fcsModel_Y.actuatorsPwmCmds[0] = rtb_ImpAsg_InsertedFor_mtrPwmCm[0];
  fcsModel_Y.actuatorsPwmCmds[1] = rtb_ImpAsg_InsertedFor_mtrPwmCm[1];
  fcsModel_Y.actuatorsPwmCmds[2] = rtb_ImpAsg_InsertedFor_mtrPwmCm[2];
  fcsModel_Y.actuatorsPwmCmds[3] = rtb_ImpAsg_InsertedFor_mtrPwmCm[3];

  // BusCreator: '<Root>/Bus Creator1'
  fcsModel_Y.fcsDebug.sysIdDebug.chirpTrigger = rtb_chirpTrigger;
  fcsModel_Y.fcsDebug.sysIdDebug.chirpType = rtb_chirpType;
  fcsModel_Y.fcsDebug.sysIdDebug.chirpSignal = yCmd;

  // RateTransition: '<Root>/Rate Transition'
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_Y.fcsDebug.outerLoopCtrlDebug = fcsModel_DW.RateTransition_Buffer0;

    // Switch: '<S3>/Switch3' incorporates:
    //   RateTransition: '<Root>/Rate Transition'
    //   Switch: '<S3>/Switch1'

    if (rtb_Compare_d) {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S166>/Bus Assignment'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    } else if (rtb_Compare_mi) {
      // Switch: '<S3>/Switch1' incorporates:
      //   BusCreator: '<S3>/Bus Creator'

      rtb_BusCreator_b_frcCmd_N = rtb_DiscreteTransferFcn;
    } else {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S166>/Bus Assignment'
      //   Switch: '<S3>/Switch1'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    }

    // End of Switch: '<S3>/Switch3'

    // BusCreator: '<S3>/Bus Creator' incorporates:
    //   BusCreator: '<S109>/Bus Creator'
    //   BusCreator: '<S110>/Bus Creator'
    //   BusCreator: '<S182>/Bus Creator'
    //   ForEachSliceAssignment generated from: '<S114>/cmd'
    //   ForEachSliceAssignment generated from: '<S167>/filtCmd'

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
    rtb_BusCreator_b_zAccelCtrlDebu = rtb_DifferenceInputs2_i;
    rtb_BusCreator_b_zAccelCtrlDe_0 = rtb_DiscreteTransferFcn_k;
    rtb_BusCreator_b_zAccelCtrlDe_1 = rtb_BusCreator_em.output;
    rtb_BusCreator_b_zAccelCtrlDe_2 = rtb_BusCreator_em.proportionalOutput;
    rtb_BusCreator_b_zAccelCtrlDe_3 = rtb_BusCreator_em.integralOutput;
    rtb_BusCreator_b_zAccelCtrlDe_4 = rtb_BusCreator_em.derivativeOutput;
    rtb_BusCreator_b_xyBodyAccelCtr = fcsModel_rtZbusXyBodyAccelCtrIDebug;
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

  // Update for UnitDelay: '<S4>/Unit Delay'
  fcsModel_DW.UnitDelay_DSTATE_g = rtb_chirpTrigger;

  // Update for RateTransition: '<Root>/Rate Transition' incorporates:
  //   BusCreator: '<S3>/Bus Creator'
  //
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // Update for DiscreteTransferFcn: '<S179>/Discrete Transfer Fcn'
    fcsModel_DW.DiscreteTransferFcn_states = fcsModel_DW.DiscreteTransferFcn_tmp;

    // Update for UnitDelay: '<S182>/Unit Delay'
    fcsModel_DW.UnitDelay_DSTATE = rtb_DifferenceInputs2_jb;
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

  // Update for UnitDelay: '<S15>/Unit Delay1'
  fcsModel_DW.UnitDelay1_DSTATE = plim;

  // Update for UnitDelay: '<S15>/Unit Delay'
  fcsModel_DW.UnitDelay_DSTATE_j = rlim;

  // Update for RandomNumber: '<S16>/White Noise'
  fcsModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&fcsModel_DW.RandSeed);

  // Update for DiscreteTransferFcn: '<S15>/Discrete Transfer Fcn'
  fcsModel_DW.DiscreteTransferFcn_states_c = vxCmd_unitRange;

  // Update for DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
  fcsModel_DW.DiscreteTransferFcn_states_d[0] = DiscreteTransferFcn_tmp_b[0];
  fcsModel_DW.DiscreteTransferFcn_states_d[1] = DiscreteTransferFcn_tmp_b[1];
  fcsModel_DW.DiscreteTransferFcn_states_d[2] = DiscreteTransferFcn_tmp_b[2];
  fcsModel_DW.DiscreteTransferFcn_states_d[3] = DiscreteTransferFcn_tmp_b[3];
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

    // InitializeConditions for RandomNumber: '<S16>/White Noise'
    fcsModel_DW.RandSeed = 1529675776U;
    fcsModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&fcsModel_DW.RandSeed);

    // 'interpretRcInputs_function:25' throttle_is_up = false;
    // 'interpretRcInputs_function:26' chirpCount_ = uint8(0);
    // 'holdOutputAtCenter_function:6' last_input = 0;
    // SystemInitialize for Iterator SubSystem: '<S109>/NED Position Control'
    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // SystemInitialize for Iterator SubSystem: '<S109>/NED Position Control'
      // SystemInitialize for Atomic SubSystem: '<S114>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S114>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S114>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S114>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S114>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S114>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S109>/NED Position Control'
    }

    // End of SystemInitialize for SubSystem: '<S109>/NED Position Control'
    // SystemInitialize for Iterator SubSystem: '<S110>/For Each Subsystem'
    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // SystemInitialize for Iterator SubSystem: '<S110>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S167>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S167>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S167>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S167>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S167>/Signal Conditioning Block2' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock2);

      // End of SystemInitialize for SubSystem: '<S167>/Signal Conditioning Block2' 

      // SystemInitialize for Atomic SubSystem: '<S167>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S167>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S110>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S110>/For Each Subsystem'
    // SystemInitialize for Atomic SubSystem: '<S182>/Signal Conditioning Block1' 
    // SystemInitialize for MATLAB Function: '<S202>/Compute Filter Numerator And Denominator' 
    // 'holdOutputAtCenter_function:6' last_input = 0;
    ComputeFilterNumeratorAndD_Init(&num[0], &den[0]);

    // End of SystemInitialize for SubSystem: '<S182>/Signal Conditioning Block1' 

    // SystemInitialize for Atomic SubSystem: '<S182>/Signal Conditioning Block' 
    SignalConditioningBlock1_c_Init(&fcsModel_DW.SignalConditioningBlock);

    // End of SystemInitialize for SubSystem: '<S182>/Signal Conditioning Block' 

    // SystemInitialize for Atomic SubSystem: '<S182>/pidWithDebug'
    fcsModel_pidWithDebug_m_Init(&fcsModel_DW.pidWithDebug);

    // End of SystemInitialize for SubSystem: '<S182>/pidWithDebug'
    // SystemInitialize for Iterator SubSystem: '<S20>/Attitude Control'
    for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
      // SystemInitialize for Iterator SubSystem: '<S20>/Attitude Control'
      // SystemInitialize for Atomic SubSystem: '<S62>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S62>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S62>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S62>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S62>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S62>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S20>/Attitude Control'
    }

    // End of SystemInitialize for SubSystem: '<S20>/Attitude Control'
    // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
    // SystemInitialize for Iterator SubSystem: '<S18>/For Each Subsystem'
    for (ForEach_itr_p = 0; ForEach_itr_p < 3; ForEach_itr_p++) {
      // SystemInitialize for Atomic SubSystem: '<S2>/Angular Rate Controller'
      // SystemInitialize for Iterator SubSystem: '<S18>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S21>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S21>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S21>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S21>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S21>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys_a[ForEach_itr_p].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S21>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S18>/For Each Subsystem'
      // End of SystemInitialize for SubSystem: '<S2>/Angular Rate Controller'
    }

    // End of SystemInitialize for SubSystem: '<S18>/For Each Subsystem'
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
