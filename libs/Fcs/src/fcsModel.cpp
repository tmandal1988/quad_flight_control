//
// File: fcsModel.cpp
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.100
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Mon Apr  8 11:00:11 2024
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

const busZaccelCtrlDebug fcsModel_rtZbusZaccelCtrlDebug{
  0.0,                                 // cmd
  0.0,                                 // meas

  {
    0.0,                               // output
    0.0,                               // proportionalOutput
    0.0,                               // integralOutput
    0.0                                // derivativeOutput
  }                                    // pidDebug
} ;                                    // busZaccelCtrlDebug ground

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
    }                                  // zAccelCtrlDebug
  },                                   // outerLoopCtrlDebug
  enumStateMachine::INACTIVE,          // state
  enumFlightMode::STABILIZE            // flightMode
};

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
//    '<S119>/Discrete First Order Deriv Filter'
//    '<S165>/Discrete First Order Deriv Filter'
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
// System initialize for atomic system:
//    '<S10>/pidWithDebug'
//    '<S51>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_Init(DW_pidWithDebug_fcsModel_T *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S13>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S10>/pidWithDebug'
//    '<S51>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc, const
  busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
  busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s, DW_pidWithDebug_fcsModel_T
  *localDW)
{
  real_T rtb_Product5_o;
  real_T rtb_Sum1;
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
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S13>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S13>/Discrete-Time Integrator'

  rtb_Sum1 = ((rtu_feedForward + rtb_Product5_o) + rtb_UnitDelay_i) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S46>/Switch2' incorporates:
  //   RelationalOperator: '<S46>/LowerRelop1'
  //   RelationalOperator: '<S46>/UpperRelop'
  //   Switch: '<S46>/Switch'

  if (rtb_Sum1 > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2 = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1 < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S46>/Switch'
    rtb_Switch2 = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2 = rtb_Sum1;
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

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
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
  localDW->UnitDelay1_DSTATE = rtb_Sum1;
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
//    '<S135>/Compute Natural Frequency'
//    '<S136>/Compute Natural Frequency'
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
//    '<S135>/Compute Numerator And Denominator'
//    '<S120>/Compute Numerator And Denominator'
//    '<S196>/Compute Numerator And Denominator'
//    '<S181>/Compute Numerator And Denominator'
//    '<S166>/Compute Numerator And Denominator'
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
//    '<S136>/Compute Filter Numerator And Denominator'
//    '<S121>/Compute Filter Numerator And Denominator'
//    '<S197>/Compute Filter Numerator And Denominator'
//    '<S182>/Compute Filter Numerator And Denominator'
//    '<S167>/Compute Filter Numerator And Denominator'
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
//    '<S136>/Compute Filter Numerator And Denominator'
//    '<S121>/Compute Filter Numerator And Denominator'
//    '<S197>/Compute Filter Numerator And Denominator'
//    '<S182>/Compute Filter Numerator And Denominator'
//    '<S167>/Compute Filter Numerator And Denominator'
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
//    '<S51>/Signal Conditioning Block1'
//    '<S51>/Signal Conditioning Block'
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
//    '<S51>/Signal Conditioning Block1'
//    '<S51>/Signal Conditioning Block'
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
//    '<S102>/holdOutputAtCenter1'
//    '<S102>/holdOutputAtCenter2'
//
void fcsModel::fcsModel_holdOutputAtCenter1(real_T rtu_input, real_T rtu_trigger,
  real_T *rty_output, boolean_T *rty_atCenter, DW_holdOutputAtCenter1_fcsMod_T
  *localDW)
{
  // MATLAB Function: '<S112>/holdOutputAtCenter'
  // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S115>:1'
  // '<S115>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
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

  // End of MATLAB Function: '<S112>/holdOutputAtCenter'
}

//
// System initialize for atomic system:
//    '<S103>/pidWithDebug'
//    '<S156>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_m_Init(DW_pidWithDebug_fcsModel_i_T
  *localDW)
{
  // InitializeConditions for DiscreteIntegrator: '<S119>/Discrete-Time Integrator' 
  localDW->DiscreteTimeIntegrator_IC_LOADI = 1U;
}

//
// Output and update for atomic system:
//    '<S103>/pidWithDebug'
//    '<S156>/pidWithDebug'
//
void fcsModel::fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
  real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc, const
  busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
  busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
  DW_pidWithDebug_fcsModel_i_T *localDW)
{
  real_T rtb_Product5_e;
  real_T rtb_Sum1;
  real_T rtb_Sum_b;
  real_T rtb_Switch2;
  real_T rtb_Switch2_n;
  real_T rtb_UkYk1;
  real_T rtb_UnitDelay_a;

  // Product: '<S151>/delta rise limit' incorporates:
  //   SampleTimeMath: '<S151>/sample time'
  //
  //  About '<S151>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_pidParamBus->outputRateLimits[1] * 0.02;

  // Sum: '<S119>/Sum'
  rtb_Sum_b = rtu_cmd - rtu_meas;

  // Outputs for Atomic SubSystem: '<S119>/Discrete First Order Deriv Filter'
  f_DiscreteFirstOrderDerivFilter(rtb_Sum_b,
    rtu_pidParamBus->filterBandwidth_radps, &rtb_Product5_e, rtp_sampleTime_s,
    &localDW->DiscreteFirstOrderDerivFilter);

  // End of Outputs for SubSystem: '<S119>/Discrete First Order Deriv Filter'

  // Product: '<S119>/Product'
  rtb_Product5_e *= rtu_pidParamBus->Kd;

  // Product: '<S119>/Product1'
  rtb_UnitDelay_a = rtb_Sum_b * rtu_pidParamBus->Kp;

  // DiscreteIntegrator: '<S119>/Discrete-Time Integrator'
  if (localDW->DiscreteTimeIntegrator_IC_LOADI != 0) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  if (rtu_integratorReset || (localDW->DiscreteTimeIntegrator_PrevRese != 0)) {
    localDW->DiscreteTimeIntegrator_DSTATE = rtu_integratorIc;
  }

  // Sum: '<S119>/Sum1' incorporates:
  //   DiscreteIntegrator: '<S119>/Discrete-Time Integrator'

  rtb_Sum1 = ((rtu_feedForward + rtb_Product5_e) + rtb_UnitDelay_a) +
    localDW->DiscreteTimeIntegrator_DSTATE;

  // Switch: '<S152>/Switch2' incorporates:
  //   RelationalOperator: '<S152>/LowerRelop1'
  //   RelationalOperator: '<S152>/UpperRelop'
  //   Switch: '<S152>/Switch'

  if (rtb_Sum1 > rtu_pidParamBus->outputLimits[1]) {
    rtb_Switch2_n = rtu_pidParamBus->outputLimits[1];
  } else if (rtb_Sum1 < rtu_pidParamBus->outputLimits[0]) {
    // Switch: '<S152>/Switch'
    rtb_Switch2_n = rtu_pidParamBus->outputLimits[0];
  } else {
    rtb_Switch2_n = rtb_Sum1;
  }

  // End of Switch: '<S152>/Switch2'

  // Sum: '<S151>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S151>/Delay Input2'
  //
  //  Block description for '<S151>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S151>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_UkYk1 = rtb_Switch2_n - localDW->DelayInput2_DSTATE;

  // Switch: '<S154>/Switch2' incorporates:
  //   RelationalOperator: '<S154>/LowerRelop1'

  if (rtb_UkYk1 <= rtb_Switch2) {
    // Product: '<S151>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S151>/sample time'
    //
    //  About '<S151>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_pidParamBus->outputRateLimits[0] * 0.02;

    // Switch: '<S154>/Switch' incorporates:
    //   RelationalOperator: '<S154>/UpperRelop'

    if (rtb_UkYk1 >= rtb_Switch2) {
      rtb_Switch2 = rtb_UkYk1;
    }

    // End of Switch: '<S154>/Switch'
  }

  // End of Switch: '<S154>/Switch2'

  // Sum: '<S151>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S151>/Delay Input2'
  //
  //  Block description for '<S151>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S151>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_ctrlCmd = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // BusCreator: '<S119>/Bus Creator' incorporates:
  //   DiscreteIntegrator: '<S119>/Discrete-Time Integrator'

  rty_pidDebug->output = *rty_ctrlCmd;
  rty_pidDebug->proportionalOutput = rtb_UnitDelay_a;
  rty_pidDebug->integralOutput = localDW->DiscreteTimeIntegrator_DSTATE;
  rty_pidDebug->derivativeOutput = rtb_Product5_e;

  // Update for DiscreteIntegrator: '<S119>/Discrete-Time Integrator' incorporates:
  //   Product: '<S119>/Product2'
  //   Product: '<S119>/Product3'
  //   Product: '<S119>/Product5'
  //   Sum: '<S119>/Sum2'
  //   Sum: '<S119>/Sum3'
  //   Sum: '<S119>/Sum4'
  //   Sum: '<S119>/Sum5'
  //   UnitDelay: '<S119>/Unit Delay'
  //   UnitDelay: '<S119>/Unit Delay1'

  localDW->DiscreteTimeIntegrator_IC_LOADI = 0U;
  localDW->DiscreteTimeIntegrator_DSTATE += (((rtu_trackingCtrlCmd -
    localDW->UnitDelay_DSTATE) * rtu_pidParamBus->Kt +
    (localDW->UnitDelay_DSTATE - localDW->UnitDelay1_DSTATE) *
    rtu_pidParamBus->Kb) + rtb_Sum_b * rtu_pidParamBus->Ki) * 0.02;
  localDW->DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>
    (rtu_integratorReset);

  // Update for UnitDelay: '<S151>/Delay Input2'
  //
  //  Block description for '<S151>/Delay Input2':
  //
  //   Store in Global RAM

  localDW->DelayInput2_DSTATE = *rty_ctrlCmd;

  // Update for UnitDelay: '<S119>/Unit Delay'
  localDW->UnitDelay_DSTATE = rtb_Switch2_n;

  // Update for UnitDelay: '<S119>/Unit Delay1'
  localDW->UnitDelay1_DSTATE = rtb_Sum1;
}

//
// System initialize for atomic system:
//    '<S103>/Signal Conditioning Block1'
//    '<S103>/Signal Conditioning Block'
//    '<S156>/Signal Conditioning Block2'
//    '<S156>/Signal Conditioning Block1'
//    '<S156>/Signal Conditioning Block'
//
void fcsModel::SignalConditioningBlock1_c_Init(DW_SignalConditioningBlock1_g_T
  *localDW)
{
  // SystemInitialize for MATLAB Function: '<S136>/Compute Filter Numerator And Denominator' 
  ComputeFilterNumeratorAndD_Init(&localDW->num[0], &localDW->den[0]);
}

//
// Output and update for atomic system:
//    '<S103>/Signal Conditioning Block1'
//    '<S103>/Signal Conditioning Block'
//    '<S156>/Signal Conditioning Block2'
//    '<S156>/Signal Conditioning Block1'
//    '<S156>/Signal Conditioning Block'
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

  // MATLAB Function: '<S135>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S135>/Compute Numerator And Denominator'
  ComputeNumeratorAndDenominator(rtb_Switch2,
    rtu_params->filterParams.dampingRatio_nd, &rtb_rateNum[0], &rtb_accelNum[0],
    &rtb_den[0], rtp_sampleTime_s);

  // MATLAB Function: '<S136>/Compute Natural Frequency'
  fcsMode_ComputeNaturalFrequency(rtu_params->filterParams.filterBandwidth_radps,
    rtu_params->filterParams.dampingRatio_nd, &rtb_Switch2);

  // MATLAB Function: '<S136>/Compute Filter Numerator And Denominator'
  ComputeFilterNumeratorAndDenomi(rtb_Switch2,
    rtu_params->filterParams.dampingRatio_nd, &localDW->num[0], &localDW->den[0],
    rtp_sampleTime_s);

  // DiscreteTransferFcn: '<S136>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_tmp = (rtu_input -
    localDW->DiscreteTransferFcn_states[0] * localDW->den[1]) -
    localDW->DiscreteTransferFcn_states[1] * localDW->den[2];
  rtb_DiscreteTransferFcn_d = (localDW->num[0] *
    localDW->DiscreteTransferFcn_tmp + localDW->DiscreteTransferFcn_states[0] *
    localDW->num[1]) + localDW->DiscreteTransferFcn_states[1] * localDW->num[2];

  // Switch: '<S140>/Switch2' incorporates:
  //   RelationalOperator: '<S140>/LowerRelop1'
  //   RelationalOperator: '<S140>/UpperRelop'
  //   Switch: '<S140>/Switch'

  if (rtb_DiscreteTransferFcn_d > rtu_params->filteredInputLimits[1]) {
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[1];
  } else if (rtb_DiscreteTransferFcn_d < rtu_params->filteredInputLimits[0]) {
    // Switch: '<S140>/Switch'
    rtb_DiscreteTransferFcn_d = rtu_params->filteredInputLimits[0];
  }

  // End of Switch: '<S140>/Switch2'

  // Sum: '<S137>/Difference Inputs1' incorporates:
  //   UnitDelay: '<S137>/Delay Input2'
  //
  //  Block description for '<S137>/Difference Inputs1':
  //
  //   Add in CPU
  //
  //  Block description for '<S137>/Delay Input2':
  //
  //   Store in Global RAM

  rtb_DiscreteTransferFcn_d -= localDW->DelayInput2_DSTATE;

  // Switch: '<S147>/Switch2' incorporates:
  //   Product: '<S137>/delta rise limit'
  //   SampleTimeMath: '<S137>/sample time'
  //
  //  About '<S137>/sample time':
  //   y = K where K = ( w * Ts )

  rtb_Switch2 = rtu_params->filteredInputRateLimits[1] * 0.02;

  // Switch: '<S147>/Switch2' incorporates:
  //   RelationalOperator: '<S147>/LowerRelop1'

  if (rtb_DiscreteTransferFcn_d <= rtb_Switch2) {
    // Product: '<S137>/delta fall limit' incorporates:
    //   SampleTimeMath: '<S137>/sample time'
    //
    //  About '<S137>/sample time':
    //   y = K where K = ( w * Ts )

    rtb_Switch2 = rtu_params->filteredInputRateLimits[0] * 0.02;

    // Switch: '<S147>/Switch' incorporates:
    //   RelationalOperator: '<S147>/UpperRelop'

    if (rtb_DiscreteTransferFcn_d >= rtb_Switch2) {
      // Switch: '<S147>/Switch2'
      rtb_Switch2 = rtb_DiscreteTransferFcn_d;
    }

    // End of Switch: '<S147>/Switch'
  }

  // End of Switch: '<S147>/Switch2'

  // Sum: '<S137>/Difference Inputs2' incorporates:
  //   UnitDelay: '<S137>/Delay Input2'
  //
  //  Block description for '<S137>/Difference Inputs2':
  //
  //   Add in CPU
  //
  //  Block description for '<S137>/Delay Input2':
  //
  //   Store in Global RAM

  *rty_filteredInput = rtb_Switch2 + localDW->DelayInput2_DSTATE;

  // Update for DiscreteTransferFcn: '<S136>/Discrete Transfer Fcn'
  localDW->DiscreteTransferFcn_states[1] = localDW->DiscreteTransferFcn_states[0];
  localDW->DiscreteTransferFcn_states[0] = localDW->DiscreteTransferFcn_tmp;

  // Update for UnitDelay: '<S137>/Delay Input2'
  //
  //  Block description for '<S137>/Delay Input2':
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

  // MATLAB Function 'checkRcCmds': '<S216>:7'
  // '<S216>:7:2' pwmLowVal = paramsStruct.pwmLimits(1);
  // '<S216>:7:3' if(rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle(1) && ... 
  // '<S216>:7:4'        rcCmds.joystickYCmd_nd <= pwmLowVal && ...
  // '<S216>:7:5'        rcCmds.joystickXCmd_nd <= pwmLowVal && ...
  // '<S216>:7:6'        rcCmds.joystickZCmd_nd <= pwmLowVal)
  if (BusConversion_InsertedFor_Chart->throttleCmd_nd <= 1000) {
    if (BusConversion_InsertedFor_Chart->joystickYCmd_nd <= 1000) {
      if (BusConversion_InsertedFor_Chart->joystickXCmd_nd <= 1000) {
        if (BusConversion_InsertedFor_Chart->joystickZCmd_nd <= 1000) {
          // '<S216>:7:7' isTrue = true;
          isTrue = true;
        } else {
          // '<S216>:7:8' else
          // '<S216>:7:9' isTrue = false;
          isTrue = false;
        }
      } else {
        // '<S216>:7:8' else
        // '<S216>:7:9' isTrue = false;
        isTrue = false;
      }
    } else {
      // '<S216>:7:8' else
      // '<S216>:7:9' isTrue = false;
      isTrue = false;
    }
  } else {
    // '<S216>:7:8' else
    // '<S216>:7:9' isTrue = false;
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
  int32_T ForEach_itr;
  int32_T ForEach_itr_i;
  std::array<real_T, 4> DiscreteTransferFcn_tmp;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angAccel;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_angRateC;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_cmd_at_i;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtCmd_;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_filtMeas;
  std::array<real_T, 3> rtb_ImpAsg_InsertedFor_meas_at_;
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
  busPidDebug rtb_BusCreator_og;
  busZaccelCtrlDebug rtb_BusCreator_b_zAccelCtrlDebu;
  real_T pCmd;
  real_T plim;
  real_T rCmd;
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
  real_T rtb_Product6;
  real_T rtb_frcCmd_N;
  real_T tCmd;
  real_T tCmd_unitRange;
  real_T tlim;
  real_T vxCmd_unitRange;
  real_T vyCmd_unitRange;
  real_T yCmd;
  int32_T ii;
  boolean_T resetIntegrator;
  boolean_T rtb_Compare_d;
  boolean_T rtb_Compare_f;
  enumFlightMode flightMode;
  enumStateMachine state;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // Math: '<S102>/Transpose' incorporates:
    //   Inport: '<Root>/stateEstimate'

    ii = 0;
    for (int32_T i{0}; i < 3; i++) {
      rtb_Transpose[ii] = fcsModel_U.stateEstimate.ned2FepDcm_nd[i];
      rtb_Transpose[ii + 1] = fcsModel_U.stateEstimate.ned2FepDcm_nd[i + 3];
      rtb_Transpose[ii + 2] = fcsModel_U.stateEstimate.ned2FepDcm_nd[i + 6];
      ii += 3;
    }

    // End of Math: '<S102>/Transpose'
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
    // Transition: '<S216>:2'
    fcsModel_DW.durationCounter_1 = 0;
    fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

    // Entry 'INACTIVE': '<S216>:1'
    // '<S216>:1:2' state = enumStateMachine.INACTIVE;
    state = enumStateMachine::INACTIVE;

    // '<S216>:1:3' rcCheckFlag = checkRcCmds;
    fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
    if (!fcsModel_DW.rcCheckFlag) {
      fcsModel_DW.durationCounter_1_j = 0;
    }

    // '<S216>:1:4' resetIntegrator = true;
    resetIntegrator = true;
  } else {
    switch (fcsModel_DW.is_c1_rcInterpreter) {
     case fcsModel_IN_ARM_MTRS:
      // During 'ARM_MTRS': '<S216>:3'
      // '<S216>:10:1' sf_internal_predicateOutput = after(60, sec) || duration(rcCheckFlag == true, sec) >= 5; 
      if (fcsModel_DW.temporalCounter_i1 >= 15000U) {
        resetIntegrator = true;
      } else {
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        resetIntegrator = (fcsModel_DW.durationCounter_1_j >= 1250);
      }

      if (resetIntegrator) {
        // Transition: '<S216>:10'
        fcsModel_DW.durationCounter_1 = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INACTIVE;

        // Entry 'INACTIVE': '<S216>:1'
        // '<S216>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S216>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S216>:1:4' resetIntegrator = true;

        // '<S216>:12:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd > paramsStruct.pwmLimitsThrottle(1); 
      } else if (fcsModel_U.rcCmdsIn.throttleCmd_nd > 1000) {
        // Transition: '<S216>:12'
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_INFLIGHT;

        // Entry 'INFLIGHT': '<S216>:11'
        // '<S216>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S216>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S216>:11:4' resetIntegrator = false;
      } else {
        // '<S216>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S216>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S216>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     case fcsModel_IN_INACTIVE:
      // During 'INACTIVE': '<S216>:1'
      // '<S216>:5:1' sf_internal_predicateOutput = duration(rcCheckFlag, sec) >= 1 && rcCmds.throttleCmd_nd >= 900; 
      if (!fcsModel_DW.rcCheckFlag) {
        fcsModel_DW.durationCounter_1 = 0;
      }

      if ((fcsModel_DW.durationCounter_1 >= 250) &&
          (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 900)) {
        // Transition: '<S216>:5'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S216>:3'
        // '<S216>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S216>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S216>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S216>:1:2' state = enumStateMachine.INACTIVE;
        state = enumStateMachine::INACTIVE;

        // '<S216>:1:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S216>:1:4' resetIntegrator = true;
        resetIntegrator = true;
      }
      break;

     default:
      // During 'INFLIGHT': '<S216>:11'
      // '<S216>:20:1' sf_internal_predicateOutput = rcCmds.throttleCmd_nd <= paramsStruct.pwmLimitsThrottle (1); 
      if (fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1000) {
        // Transition: '<S216>:20'
        fcsModel_DW.durationCounter_1_j = 0;
        fcsModel_DW.is_c1_rcInterpreter = fcsModel_IN_ARM_MTRS;
        fcsModel_DW.temporalCounter_i1 = 0U;

        // Entry 'ARM_MTRS': '<S216>:3'
        // '<S216>:3:2' state = enumStateMachine.MTR_ARMED;
        state = enumStateMachine::MTR_ARMED;

        // '<S216>:3:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
        }

        // '<S216>:3:4' resetIntegrator = true;
        resetIntegrator = true;
      } else {
        // '<S216>:11:2' state = enumStateMachine.INFLIGHT;
        state = enumStateMachine::INFLIGHT;

        // '<S216>:11:3' rcCheckFlag = checkRcCmds;
        fcsModel_DW.rcCheckFlag = fcsModel_checkRcCmds(&fcsModel_U.rcCmdsIn);
        if (!fcsModel_DW.rcCheckFlag) {
          fcsModel_DW.durationCounter_1 = 0;
          fcsModel_DW.durationCounter_1_j = 0;
        }

        // '<S216>:11:4' resetIntegrator = false;
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
  // MATLAB Function 'rcInterpreter/Interpret RC In Cmds': '<S217>:1'
  // '<S217>:1:3' [flightMode, rcOutCmds] = interpretRcInputs_function(rcCmds, expo, rcParamsStruct); 
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

    // 'interpretRcInputs_function:49' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -45.0;

    // 'interpretRcInputs_function:50' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:51' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:52' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:53' vxlim = rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:54' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:55' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
  } else if ((fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1100) &&
             (fcsModel_U.rcCmdsIn.rcSwitch1_nd < 1700)) {
    // 'interpretRcInputs_function:57' elseif (rcInCmds.rcSwitch1_nd >= 1100 && rcInCmds.rcSwitch1_nd < 1700) 
    // 'interpretRcInputs_function:58' flightMode = enumFlightMode.ALT_CONTROL;
    flightMode = enumFlightMode::ALT_CONTROL;

    // 'interpretRcInputs_function:59' vxlim = rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:60' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:61' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:62' tlim = 0;
    tlim = 0.0;

    // 'interpretRcInputs_function:63' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:64' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:65' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else if (fcsModel_U.rcCmdsIn.rcSwitch1_nd >= 1700) {
    // 'interpretRcInputs_function:66' elseif (rcInCmds.rcSwitch1_nd >= 1700)
    // 'interpretRcInputs_function:67' flightMode = enumFlightMode.POS_CONTROL;
    flightMode = enumFlightMode::POS_CONTROL;

    // 'interpretRcInputs_function:68' vxlim = rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:69' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:70' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:71' tlim = 0;
    tlim = 0.0;

    // 'interpretRcInputs_function:72' rlim = 0;
    rlim = 0.0;

    // 'interpretRcInputs_function:73' plim = 0;
    plim = 0.0;

    // 'interpretRcInputs_function:74' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
  } else {
    // 'interpretRcInputs_function:75' else
    // 'interpretRcInputs_function:76' flightMode = enumFlightMode.STABILIZE;
    flightMode = enumFlightMode::STABILIZE;

    // 'interpretRcInputs_function:77' tlim = -rcParamsStruct.cmdLimits.zForce_N(2); 
    tlim = -45.0;

    // 'interpretRcInputs_function:78' rlim = rcParamsStruct.cmdLimits.roll_rad(2); 
    rlim = 0.78539816339744828;

    // 'interpretRcInputs_function:79' plim = rcParamsStruct.cmdLimits.pitch_rad(2); 
    plim = 0.78539816339744828;

    // 'interpretRcInputs_function:80' ylim = rcParamsStruct.cmdLimits.yawRate_radps(2); 
    // 'interpretRcInputs_function:81' vxlim = rcParamsStruct.cmdLimits.vz_mps(2); 
    // 'interpretRcInputs_function:82' vylim = rcParamsStruct.cmdLimits.vy_mps(2); 
    // 'interpretRcInputs_function:83' vzlim = -rcParamsStruct.cmdLimits.vz_mps(2); 
  }

  // 'interpretRcInputs_function:87' if double(rcInCmds.throttleCmd_nd) < rcParamsStruct.pwmLimitsThrottle(1) 
  if (fcsModel_U.rcCmdsIn.throttleCmd_nd < 1000) {
    // This means either we haven't take off yet or we landed after a flight
    // and might take of again so set the throttle_is_up to false
    // 'interpretRcInputs_function:90' throttle_is_up = false;
    fcsModel_DW.throttle_is_up = false;
  } else if ((fcsModel_U.rcCmdsIn.throttleCmd_nd <= 1470) &&
             (fcsModel_U.rcCmdsIn.throttleCmd_nd >= 1260)) {
    // 'interpretRcInputs_function:91' elseif ((double(rcInCmds.throttleCmd_nd) <= rcParamsStruct.pwmThrottleMidHigh) && ... 
    // 'interpretRcInputs_function:92'         double(rcInCmds.throttleCmd_nd) >= rcParamsStruct.pwmThrottleMidLow) 
    // Since throttle stick on RC transmitter is half way up this probably
    // means the vehicle is flying so set throttle_is_up flag to true
    // 'interpretRcInputs_function:95' throttle_is_up = true;
    fcsModel_DW.throttle_is_up = true;
  } else {
    // 'interpretRcInputs_function:96' else
    // Otherwise preserve the previous flag
  }

  // 'interpretRcInputs_function:100' tCmd = min( rcParamsStruct.pwmLimitsThrottle(2), ... 
  // 'interpretRcInputs_function:101'         max( rcParamsStruct.pwmLimitsThrottle(1), double(rcInCmds.throttleCmd_nd) ) ); 
  tCmd = std::fmin(1900.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.throttleCmd_nd)));

  // 'interpretRcInputs_function:103' if (flightMode == enumFlightMode.STABILIZE || flightMode == enumFlightMode.ACRO) 
  if (flightMode == enumFlightMode::STABILIZE) {
    //  In stabilize mode throttle stick starts at 0
    // 'interpretRcInputs_function:105' tCmd_unitRange = -1 + tCmd/1000;
    tCmd_unitRange = tCmd / 1000.0 + -1.0;
  } else {
    // 'interpretRcInputs_function:106' else
    //  In other modes throttle is used with symmetry around middle stick
    //  position. NOT USED CURRENTLY. Defined so that tCmd_unitRange is
    //  defined in all execution path
    // 'interpretRcInputs_function:110' tCmd_unitRange = -1 + tCmd/500;
    tCmd_unitRange = tCmd / 500.0 + -1.0;
  }

  //  Use throttle stick to set Vz commands. Always set but only used
  //  in ALT_CONTROL flight mode. Has different slopes about the center point
  //  as the center point is not always at 1500 which is PWM center
  // 'interpretRcInputs_function:117' if ((tCmd <= rcParamsStruct.pwmThrottleMidHigh) && ... 
  // 'interpretRcInputs_function:118'         tCmd >= rcParamsStruct.pwmThrottleMidLow) 
  if ((tCmd <= 1470.0) && (tCmd >= 1260.0)) {
    // 'interpretRcInputs_function:119' vzCmd_unitRange  = 0;
    tCmd = 0.0;
  } else if (tCmd < 1260.0) {
    // 'interpretRcInputs_function:120' elseif (tCmd < rcParamsStruct.pwmThrottleMidLow) 
    // 'interpretRcInputs_function:121' if (throttle_is_up)
    if (fcsModel_DW.throttle_is_up) {
      //  This means we are already flying, make lower half of throttle
      //  stick map to descent rates
      // 'interpretRcInputs_function:124' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:125'         rcParamsStruct.pwmToCmdThrottleIncptLow; 
      tCmd = 0.0038461538461538464 * tCmd + -4.8461538461538458;

      // ;
    } else {
      // 'interpretRcInputs_function:126' else
      //  This means we haven't taken off yet or we landed and might take
      //  off again. Map the lower half of the throttle stick also to
      //  ascent rate as we want to climb up from the ground and only go
      //  upto 3/4 of the max ascent rate at mid point. Lowest position is
      //  throttle maps to -1 so add +1 make lowest throttle point 0 and
      //  increase from there on
      // 'interpretRcInputs_function:133' vzCmd_unitRange = ((rcParamsStruct.pwmToCmdThrottleSlopeLow*tCmd + ...; 
      // 'interpretRcInputs_function:134'             rcParamsStruct.pwmToCmdThrottleIncptLow) + 1)*0.75; 
      tCmd = ((0.0038461538461538464 * tCmd + -4.8461538461538458) + 1.0) * 0.75;

      // ;
    }
  } else {
    // 'interpretRcInputs_function:136' else
    // 'interpretRcInputs_function:137' vzCmd_unitRange = rcParamsStruct.pwmToCmdThrottleSlopeHigh*tCmd + ...; 
    // 'interpretRcInputs_function:138'         rcParamsStruct.pwmToCmdThrottleIncptHigh; 
    tCmd = 0.0024271844660194173 * tCmd + -3.5679611650485437;

    // ;
  }

  //  Set roll, pitch and yaw stick
  // 'interpretRcInputs_function:142' rCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:143'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickXCmd_nd) ) ); 
  rCmd = std::fmin(2000.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickXCmd_nd)));

  //  Use roll stick to set FEP Vy to be used for POS control mode
  // 'interpretRcInputs_function:146' if ((rCmd <= rcParamsStruct.pwmRollStickMidHigh) && ... 
  // 'interpretRcInputs_function:147'         rCmd >= rcParamsStruct.pwmRollStickMidLow) 
  if ((rCmd <= 1650.0) && (rCmd >= 1350.0)) {
    // 'interpretRcInputs_function:148' vyCmd_unitRange  = 0;
    vyCmd_unitRange = 0.0;
  } else {
    // 'interpretRcInputs_function:149' else
    // 'interpretRcInputs_function:150' vyCmd_unitRange = -3 + rCmd/500;
    vyCmd_unitRange = rCmd / 500.0 + -3.0;
  }

  // 'interpretRcInputs_function:152' rCmd_unitRange = -3 + rCmd/500;
  rCmd = rCmd / 500.0 + -3.0;

  // 'interpretRcInputs_function:155' pCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:156'         max( rcParamsStruct.pwmLimits(1),  double(rcInCmds.joystickYCmd_nd) ) ); 
  pCmd = std::fmin(2000.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickYCmd_nd)));

  //  Use pitch stick to set FEP Vx to be used for POS control mode
  // 'interpretRcInputs_function:159' if ((pCmd <= rcParamsStruct.pwmPitchStickMidHigh) && ... 
  // 'interpretRcInputs_function:160'         pCmd >= rcParamsStruct.pwmPitchStickMidLow) 
  if ((pCmd <= 1650.0) && (pCmd >= 1350.0)) {
    // 'interpretRcInputs_function:161' vxCmd_unitRange  = 0;
    vxCmd_unitRange = 0.0;
  } else {
    // 'interpretRcInputs_function:162' else
    // 'interpretRcInputs_function:163' vxCmd_unitRange = -3 + pCmd/500;
    vxCmd_unitRange = pCmd / 500.0 + -3.0;
  }

  //  Reverse the pitch cmd
  // 'interpretRcInputs_function:167' pCmd_unitRange = -(-3 + pCmd/500);
  pCmd = -(pCmd / 500.0 + -3.0);

  // 'interpretRcInputs_function:170' yCmd = min( rcParamsStruct.pwmLimits(2), ... 
  // 'interpretRcInputs_function:171'         max( rcParamsStruct.pwmLimits(1), double(rcInCmds.joystickZCmd_nd) ) ); 
  yCmd = std::fmin(2000.0, std::fmax(1000.0, static_cast<real_T>
    (fcsModel_U.rcCmdsIn.joystickZCmd_nd)));

  //  Use yaw stick to also pick a Yaw angle in ALT or POS control mode
  // 'interpretRcInputs_function:174' if (flightMode == enumFlightMode.STABILIZE) 
  if (flightMode == enumFlightMode::STABILIZE) {
    // 'interpretRcInputs_function:175' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;

    // 'interpretRcInputs_function:176' else
    // 'interpretRcInputs_function:177' if ((yCmd <= rcParamsStruct.pwmYawStickMidHigh) && ... 
    // 'interpretRcInputs_function:178'         yCmd >= rcParamsStruct.pwmYawStickMidLow) 
  } else if ((yCmd <= 1600.0) && (yCmd >= 1400.0)) {
    // 'interpretRcInputs_function:179' yCmd_unitRange  = 0;
    yCmd = 0.0;
  } else {
    // 'interpretRcInputs_function:180' else
    // 'interpretRcInputs_function:181' yCmd_unitRange = -3 + yCmd/500;
    yCmd = yCmd / 500.0 + -3.0;
  }

  // 'interpretRcInputs_function:186' if expo
  // 'interpretRcInputs_function:209' else
  //  Usually expo is set in the Tx hence simply use a linear map here
  // 'interpretRcInputs_function:211' rcOutCmds.throttleStick = tCmd_unitRange*tlim; 
  fcsModel_DW.rcOutCmds.throttleStick = tCmd_unitRange * tlim;

  // 'interpretRcInputs_function:212' rcOutCmds.rollStick = rCmd_unitRange*rlim; 
  fcsModel_DW.rcOutCmds.rollStick = rCmd * rlim;

  // 'interpretRcInputs_function:213' rcOutCmds.pitchStick = pCmd_unitRange*plim; 
  fcsModel_DW.rcOutCmds.pitchStick = pCmd * plim;

  // 'interpretRcInputs_function:214' rcOutCmds.yawStick = yCmd_unitRange*ylim;
  fcsModel_DW.rcOutCmds.yawStick = yCmd * 1.0471975511965976;

  // 'interpretRcInputs_function:215' rcOutCmds.vzStick_mps = vzCmd_unitRange*vzlim; 
  fcsModel_DW.rcOutCmds.vzStick_mps = tCmd * -2.0;

  // 'interpretRcInputs_function:216' rcOutCmds.vxStick_mps = vxCmd_unitRange*vxlim; 
  fcsModel_DW.rcOutCmds.vxStick_mps = vxCmd_unitRange * 2.0;

  // 'interpretRcInputs_function:217' rcOutCmds.vyStick_mps = vyCmd_unitRange*vylim; 
  fcsModel_DW.rcOutCmds.vyStick_mps = vyCmd_unitRange * 2.5;
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    // Product: '<S102>/Matrix Multiply' incorporates:
    //   Math: '<S102>/Transpose'
    //   SignalConversion generated from: '<S102>/Vector Concatenate1'

    for (ii = 0; ii < 3; ii++) {
      rtb_MatrixMultiply[ii] = 0.0;
      rtb_MatrixMultiply[ii] += rtb_Transpose[ii] *
        fcsModel_DW.rcOutCmds.vxStick_mps;
      rtb_MatrixMultiply[ii] += rtb_Transpose[ii + 3] *
        fcsModel_DW.rcOutCmds.vyStick_mps;
    }

    // End of Product: '<S102>/Matrix Multiply'

    // Outputs for Atomic SubSystem: '<S102>/holdOutputAtCenter1'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[0],
      rtb_MatrixMultiply[0], &rtb_frcCmd_N, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter1);

    // End of Outputs for SubSystem: '<S102>/holdOutputAtCenter1'

    // Concatenate: '<S102>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // Switch: '<S102>/Switch1' incorporates:
    //   RelationalOperator: '<S106>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S102>/Bus Assignment1' incorporates:
      //   Concatenate: '<S102>/Vector Concatenate'
      //   Constant: '<S102>/Constant1'

      rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S102>/Bus Assignment1' incorporates:
      //   Concatenate: '<S102>/Vector Concatenate'

      rtb_VectorConcatenate[0].feedForwardCmd = rtb_MatrixMultiply[0];
    }

    // End of Switch: '<S102>/Switch1'

    // BusAssignment: '<S102>/Bus Assignment1' incorporates:
    //   Concatenate: '<S102>/Vector Concatenate'
    //   Constant: '<S109>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S102>/Logical Operator2'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S109>/Compare'

    rtb_VectorConcatenate[0].cmd = rtb_frcCmd_N;
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedPos_m[0];
    rtb_VectorConcatenate[0].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S102>/holdOutputAtCenter2'
    // Inport: '<Root>/stateEstimate'
    fcsModel_holdOutputAtCenter1(fcsModel_U.stateEstimate.nedPos_m[1],
      rtb_MatrixMultiply[1], &rtb_frcCmd_N, &rtb_Compare_d,
      &fcsModel_DW.holdOutputAtCenter2);

    // End of Outputs for SubSystem: '<S102>/holdOutputAtCenter2'

    // Concatenate: '<S102>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // Switch: '<S102>/Switch2' incorporates:
    //   RelationalOperator: '<S107>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S102>/Bus Assignment2' incorporates:
      //   Concatenate: '<S102>/Vector Concatenate'
      //   Constant: '<S102>/Constant3'

      rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S102>/Bus Assignment2' incorporates:
      //   Concatenate: '<S102>/Vector Concatenate'

      rtb_VectorConcatenate[1].feedForwardCmd = rtb_MatrixMultiply[1];
    }

    // End of Switch: '<S102>/Switch2'

    // BusAssignment: '<S102>/Bus Assignment2' incorporates:
    //   Concatenate: '<S102>/Vector Concatenate'
    //   Constant: '<S110>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S102>/Logical Operator3'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S110>/Compare'

    rtb_VectorConcatenate[1].cmd = rtb_frcCmd_N;
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedPos_m[1];
    rtb_VectorConcatenate[1].integratorReset = (resetIntegrator || (flightMode
      != enumFlightMode::POS_CONTROL));

    // Outputs for Atomic SubSystem: '<S102>/holdOutputAtCenter'
    // MATLAB Function: '<S111>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S114>:1'
    // '<S114>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
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

    // End of Outputs for SubSystem: '<S102>/holdOutputAtCenter'

    // Switch: '<S102>/Switch' incorporates:
    //   Constant: '<S102>/Constant5'
    //   RelationalOperator: '<S104>/Compare'

    // 'holdOutputAtCenter_function:17' output = last_input;
    if (rtb_Compare_d) {
      rtb_frcCmd_N = 0.0;
    } else {
      rtb_frcCmd_N = fcsModel_DW.rcOutCmds.vzStick_mps;
    }

    // End of Switch: '<S102>/Switch'

    // Outputs for Atomic SubSystem: '<S102>/holdOutputAtCenter'
    // Gain: '<S102>/Gain' incorporates:
    //   MATLAB Function: '<S111>/holdOutputAtCenter'

    rtb_Product6 = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S102>/holdOutputAtCenter'

    // Gain: '<S102>/Gain1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    tlim = -fcsModel_U.stateEstimate.aglEst_m;

    // Concatenate: '<S102>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S102>/Bus Assignment3' incorporates:
    //   Concatenate: '<S102>/Vector Concatenate'
    //   Constant: '<S105>/Constant'
    //   Constant: '<S108>/Constant'
    //   Gain: '<S102>/Gain'
    //   Gain: '<S102>/Gain1'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S102>/Logical Operator'
    //   Logic: '<S102>/Logical Operator1'
    //   MATLAB Function: '<S111>/holdOutputAtCenter'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S105>/Compare'
    //   RelationalOperator: '<S108>/Compare'

    rtb_VectorConcatenate[2].feedForwardCmd = rtb_frcCmd_N;

    // Outputs for Atomic SubSystem: '<S102>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].cmd = -fcsModel_DW.last_input_c;

    // End of Outputs for SubSystem: '<S102>/holdOutputAtCenter'
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.aglEst_m;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && (flightMode != enumFlightMode::
      POS_CONTROL)));

    // Outputs for Iterator SubSystem: '<S98>/NED Position Control' incorporates:
    //   ForEach: '<S103>/For Each'

    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // Outputs for Atomic SubSystem: '<S103>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S103>/ctrlInputs' incorporates:
      //   BusAssignment: '<S102>/Bus Assignment'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S103>/Unit Delay'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr_i], &rtb_frcCmd_N, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S103>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S103>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr_i].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr_i], &tlim, 0.02, &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S103>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S103>/pidWithDebug'
      fcsModel_pidWithDebug_j(rtb_VectorConcatenate[ForEach_itr_i].
        feedForwardCmd, rtb_frcCmd_N, tlim, rtb_VectorConcatenate[ForEach_itr_i]
        .integratorReset, 0.0,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.posCtrlParams.ctrlParamsArray[
        ForEach_itr_i], fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE,
        &rtb_Product6, &rtb_BusCreator_og, 0.02,
        &fcsModel_DW.CoreSubsys_g[ForEach_itr_i].pidWithDebug);

      // End of Outputs for SubSystem: '<S103>/pidWithDebug'

      // Update for UnitDelay: '<S103>/Unit Delay'
      fcsModel_DW.CoreSubsys_g[ForEach_itr_i].UnitDelay_DSTATE = rtb_Product6;

      // ForEachSliceAssignment generated from: '<S103>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDeb_m[ForEach_itr_i] = rtb_BusCreator_og;

      // ForEachSliceAssignment generated from: '<S103>/neVelCmd_mps'
      rtb_ImpAsg_InsertedFor_neVelCmd[ForEach_itr_i] = rtb_Product6;

      // ForEachSliceAssignment generated from: '<S103>/meas'
      rtb_ImpAsg_InsertedFor_meas_at_[ForEach_itr_i] = tlim;

      // ForEachSliceAssignment generated from: '<S103>/cmd'
      rtb_ImpAsg_InsertedFor_cmd_at_i[ForEach_itr_i] = rtb_frcCmd_N;
    }

    // End of Outputs for SubSystem: '<S98>/NED Position Control'

    // RelationalOperator: '<S101>/Compare' incorporates:
    //   Constant: '<S101>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_f = (flightMode != enumFlightMode::POS_CONTROL);

    // Logic: '<S97>/Logical Operator2'
    rtb_Compare_d = (resetIntegrator || rtb_Compare_f);

    // Concatenate: '<S97>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[0], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S97>/Bus Assignment' incorporates:
    //   Concatenate: '<S97>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[0];
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.nedVel_mps[0];
    rtb_VectorConcatenate[0].integratorReset = rtb_Compare_d;

    // Concatenate: '<S97>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[1], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S97>/Bus Assignment1' incorporates:
    //   Concatenate: '<S97>/Vector Concatenate'
    //   Inport: '<Root>/stateEstimate'

    rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[1];
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.nedVel_mps[1];
    rtb_VectorConcatenate[1].integratorReset = rtb_Compare_d;

    // Gain: '<S97>/Gain' incorporates:
    //   Inport: '<Root>/stateEstimate'

    tlim = -fcsModel_U.stateEstimate.climbRateEst_mps;

    // Concatenate: '<S97>/Vector Concatenate'
    std::memset(&rtb_VectorConcatenate[2], 0, sizeof(busCtrlInputs));

    // BusAssignment: '<S97>/Bus Assignment2' incorporates:
    //   Concatenate: '<S97>/Vector Concatenate'
    //   Constant: '<S100>/Constant'
    //   Gain: '<S97>/Gain'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S97>/Logical Operator'
    //   Logic: '<S97>/Logical Operator1'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S100>/Compare'

    rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_neVelCmd[2];
    rtb_VectorConcatenate[2].meas = -fcsModel_U.stateEstimate.climbRateEst_mps;
    rtb_VectorConcatenate[2].integratorReset = (resetIntegrator || ((flightMode
      != enumFlightMode::ALT_CONTROL) && rtb_Compare_f));

    // SignalConversion generated from: '<S99>/For Each Subsystem'
    rtb_VectorConcatenate1[0] = 0.0;
    rtb_VectorConcatenate1[1] = 0.0;
    rtb_VectorConcatenate1[2] = 0.0;

    // Outputs for Iterator SubSystem: '<S99>/For Each Subsystem' incorporates:
    //   ForEach: '<S156>/For Each'

    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // Outputs for Atomic SubSystem: '<S156>/Signal Conditioning Block'
      // ForEachSliceSelector generated from: '<S156>/ctrlInputs' incorporates:
      //   BusAssignment: '<S97>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Inport: '<Root>/ctrlParams'

      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].cmd,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.cmdSignalConditioningParamsArray
        [ForEach_itr], &tlim, 0.02, &fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of Outputs for SubSystem: '<S156>/Signal Conditioning Block'

      // Outputs for Atomic SubSystem: '<S156>/Signal Conditioning Block1'
      fcsM_SignalConditioningBlock1_f(rtb_VectorConcatenate[ForEach_itr].meas,
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.measSignalConditioningParamsArray
        [ForEach_itr], &rtb_Product6, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock1);

      // End of Outputs for SubSystem: '<S156>/Signal Conditioning Block1'

      // Outputs for Atomic SubSystem: '<S156>/Signal Conditioning Block2'
      // ForEachSliceSelector generated from: '<S156>/nedAccel_mps2' incorporates:
      //   Inport: '<Root>/ctrlParams'
      //   Inport: '<Root>/stateEstimate'

      fcsM_SignalConditioningBlock1_f
        (fcsModel_U.stateEstimate.nedAccel_mps2[ForEach_itr],
         &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelSignalConditioningParamsArray
         [ForEach_itr], &rtb_frcCmd_N, 0.02,
         &fcsModel_DW.CoreSubsys_i[ForEach_itr].SignalConditioningBlock2);

      // End of Outputs for SubSystem: '<S156>/Signal Conditioning Block2'

      // Product: '<S156>/Product' incorporates:
      //   ForEachSliceSelector generated from: '<S156>/accelFbGain'
      //   Inport: '<Root>/ctrlParams'

      rlim =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelFbGainsArray
        [ForEach_itr] * rtb_frcCmd_N;

      // Product: '<S156>/Product1' incorporates:
      //   ForEachSliceSelector generated from: '<S156>/ffGain'
      //   Inport: '<Root>/ctrlParams'

      rtb_frcCmd_N =
        fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ffGainsArray[ForEach_itr]
        * tlim;

      // Outputs for Atomic SubSystem: '<S156>/pidWithDebug'
      // Sum: '<S156>/Sum' incorporates:
      //   BusAssignment: '<S97>/Bus Assignment3'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   ForEachSliceSelector generated from: '<S156>/ctrlInputs'
      //   Inport: '<Root>/ctrlParams'
      //   UnitDelay: '<S156>/Unit Delay'

      fcsModel_pidWithDebug_j(rtb_frcCmd_N - rlim, tlim, rtb_Product6,
        rtb_VectorConcatenate[ForEach_itr].integratorReset,
        rtb_VectorConcatenate1[ForEach_itr],
        &fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.ctrlParamsArray[
        ForEach_itr], fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE,
        &rtb_frcCmd_N, &rtb_BusCreator_og, 0.02,
        &fcsModel_DW.CoreSubsys_i[ForEach_itr].pidWithDebug);

      // End of Outputs for SubSystem: '<S156>/pidWithDebug'

      // Update for UnitDelay: '<S156>/Unit Delay'
      fcsModel_DW.CoreSubsys_i[ForEach_itr].UnitDelay_DSTATE = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S156>/velCtrlOut '
      rtb_ImpAsg_InsertedFor_velCtrlO[ForEach_itr] = rtb_frcCmd_N;

      // ForEachSliceAssignment generated from: '<S156>/pidDebug'
      rtb_ImpAsg_InsertedFor_pidDebug[ForEach_itr] = rtb_BusCreator_og;

      // ForEachSliceAssignment generated from: '<S156>/velCtrlFf'
      rtb_ImpAsg_InsertedFor_velCtrlF[ForEach_itr] = rlim;

      // ForEachSliceAssignment generated from: '<S156>/filtMeas'
      rtb_ImpAsg_InsertedFor_filtMeas[ForEach_itr] = rtb_Product6;

      // ForEachSliceAssignment generated from: '<S156>/filtCmd'
      rtb_ImpAsg_InsertedFor_filtCmd_[ForEach_itr] = tlim;
    }

    // End of Outputs for SubSystem: '<S99>/For Each Subsystem'

    // Trigonometry: '<S155>/Sin' incorporates:
    //   Inport: '<Root>/stateEstimate'

    tlim = std::sin(fcsModel_U.stateEstimate.attitude_rad[2]);

    // Trigonometry: '<S155>/Sin1' incorporates:
    //   Inport: '<Root>/stateEstimate'

    rtb_frcCmd_N = std::cos(fcsModel_U.stateEstimate.attitude_rad[2]);

    // BusAssignment: '<S155>/Bus Assignment1' incorporates:
    //   Constant: '<S155>/Constant4'
    //   Inport: '<Root>/ctrlParams'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S155>/Product1'
    //   Product: '<S155>/Product2'
    //   Product: '<S155>/Product6'
    //   Sum: '<S155>/Sum'

    std::memset(&rtb_BusAssignment1_a, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment1_a.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] * tlim -
      rtb_frcCmd_N * rtb_ImpAsg_InsertedFor_velCtrlO[1]) *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelCmdToAttitudeCmdScale_nd
      [0] * -0.10197838058331635;
    rtb_BusAssignment1_a.meas = fcsModel_U.stateEstimate.attitude_rad[0];
    rtb_BusAssignment1_a.integratorReset = resetIntegrator;

    // BusAssignment: '<S155>/Bus Assignment2' incorporates:
    //   Constant: '<S155>/Constant4'
    //   Inport: '<Root>/ctrlParams'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S155>/Product3'
    //   Product: '<S155>/Product4'
    //   Product: '<S155>/Product7'
    //   Sum: '<S155>/Sum1'

    std::memset(&rtb_BusAssignment2, 0, sizeof(busCtrlInputs));
    rtb_BusAssignment2.cmd = (rtb_ImpAsg_InsertedFor_velCtrlO[0] * rtb_frcCmd_N
      + tlim * rtb_ImpAsg_InsertedFor_velCtrlO[1]) *
      fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.accelCmdToAttitudeCmdScale_nd
      [1] * -0.10197838058331635;
    rtb_BusAssignment2.meas = fcsModel_U.stateEstimate.attitude_rad[1];
    rtb_BusAssignment2.integratorReset = resetIntegrator;

    // Outputs for Atomic SubSystem: '<S155>/holdOutputAtCenter'
    // MATLAB Function: '<S160>/holdOutputAtCenter' incorporates:
    //   Inport: '<Root>/stateEstimate'

    // MATLAB Function 'holdOutputAtCenter/holdOutputAtCenter': '<S161>:1'
    // '<S161>:1:2' [output, atCenter] = holdOutputAtCenter_function(input, trigger, params); 
    // HOLDOUTPUTATCENTER_FUNCTION holds the output constant at last input if the 
    // trigger value is within user defined delta from the center
    // 'holdOutputAtCenter_function:5' if isempty(last_input)
    // 'holdOutputAtCenter_function:9' if(trigger <= (params.center + params.posDeltaFromCenter) && ... 
    // 'holdOutputAtCenter_function:10'         trigger >=(params.center - params.negDeltaFromCenter)) 
    if ((fcsModel_DW.rcOutCmds.yawStick <= 0.017453292519943295) &&
        (fcsModel_DW.rcOutCmds.yawStick >= -0.017453292519943295)) {
      // 'holdOutputAtCenter_function:11' atCenter = true;
      rtb_Compare_d = true;
    } else {
      // 'holdOutputAtCenter_function:12' else
      // 'holdOutputAtCenter_function:13' atCenter = false;
      rtb_Compare_d = false;

      // 'holdOutputAtCenter_function:14' last_input = input;
      fcsModel_DW.last_input = fcsModel_U.stateEstimate.attitude_rad[2];
    }

    // End of Outputs for SubSystem: '<S155>/holdOutputAtCenter'

    // BusAssignment: '<S155>/Bus Assignment4'
    // 'holdOutputAtCenter_function:17' output = last_input;
    std::memset(&rtb_BusAssignment4, 0, sizeof(busCtrlInputs));

    // Switch: '<S155>/Switch' incorporates:
    //   RelationalOperator: '<S157>/Compare'

    if (rtb_Compare_d) {
      // BusAssignment: '<S155>/Bus Assignment4' incorporates:
      //   Constant: '<S155>/Constant5'

      rtb_BusAssignment4.feedForwardCmd = 0.0;
    } else {
      // BusAssignment: '<S155>/Bus Assignment4'
      rtb_BusAssignment4.feedForwardCmd = fcsModel_DW.rcOutCmds.yawStick;
    }

    // End of Switch: '<S155>/Switch'

    // Outputs for Atomic SubSystem: '<S155>/holdOutputAtCenter'
    // BusAssignment: '<S155>/Bus Assignment4' incorporates:
    //   Constant: '<S158>/Constant'
    //   Constant: '<S159>/Constant'
    //   Inport: '<Root>/stateEstimate'
    //   Logic: '<S155>/Logical Operator'
    //   Logic: '<S155>/Logical Operator1'
    //   MATLAB Function: '<S160>/holdOutputAtCenter'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'
    //   RelationalOperator: '<S158>/Compare'
    //   RelationalOperator: '<S159>/Compare'

    rtb_BusAssignment4.cmd = fcsModel_DW.last_input;

    // End of Outputs for SubSystem: '<S155>/holdOutputAtCenter'
    rtb_BusAssignment4.meas = fcsModel_U.stateEstimate.attitude_rad[2];
    rtb_BusAssignment4.integratorReset = (resetIntegrator || ((flightMode !=
      enumFlightMode::ALT_CONTROL) && (flightMode != enumFlightMode::POS_CONTROL)));

    // Product: '<S155>/Divide1' incorporates:
    //   Constant: '<S155>/g'
    //   Inport: '<Root>/ctrlParams'
    //   Inport: '<Root>/stateEstimate'
    //   Product: '<S155>/Divide3'
    //   Product: '<S155>/Product5'
    //   Sum: '<S155>/Sum2'
    //   Trigonometry: '<S155>/Sin2'
    //   Trigonometry: '<S155>/Sin3'

    rtb_frcCmd_N = 1.0 / (std::cos(fcsModel_U.stateEstimate.attitude_rad[0]) *
                          std::cos(fcsModel_U.stateEstimate.attitude_rad[1])) *
      (fcsModel_U.ctrlParams.outerLoopCtrlParams.velCtrlParams.baseMass_kg *
       -9.806 + rtb_ImpAsg_InsertedFor_velCtrlO[2]);

    // RelationalOperator: '<S94>/Compare' incorporates:
    //   Constant: '<S94>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_d = (flightMode == enumFlightMode::POS_CONTROL);

    // MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus' incorporates:
    //   BusCreator generated from: '<S3>/assembleOuterLoopToInnerLoopBus'
    //   Constant: '<S3>/Constant'
    //   Inport: '<Root>/stateEstimate'

    std::memcpy(&rtb_VectorConcatenate[0],
                &fcsModel_ConstP.pooled3.attCtrlInputs.ctrlInputsArray[0], 3U *
                sizeof(busCtrlInputs));

    // MATLAB Function 'Outer Loop Controller/assembleOuterLoopToInnerLoopBus': '<S96>:1' 
    // '<S96>:1:2' outBus.outerLoopCmds.thrustCmd_N = throttleCmd_N;
    // '<S96>:1:3' outDebug = throttleCmd_N;
    rtb_Product6 = fcsModel_DW.rcOutCmds.throttleStick;

    //  This is a stop gap setup where we are only assuming that rate control
    //  is active and therefore not setting up attCtrlInputs for Euler angle
    //  control
    // '<S96>:1:7' outBus.attCtrlInputs.ctrlInputsArray(1).cmd = rcOutCmds.rollStick; 
    rtb_VectorConcatenate[0].cmd = fcsModel_DW.rcOutCmds.rollStick;

    // '<S96>:1:8' outBus.attCtrlInputs.ctrlInputsArray(1).meas = stateEstimate.attitude_rad(1); 
    rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.attitude_rad[0];

    // '<S96>:1:9' outBus.attCtrlInputs.ctrlInputsArray(2).cmd = rcOutCmds.pitchStick; 
    rtb_VectorConcatenate[1].cmd = fcsModel_DW.rcOutCmds.pitchStick;

    // '<S96>:1:10' outBus.attCtrlInputs.ctrlInputsArray(2).meas = stateEstimate.attitude_rad(2); 
    rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.attitude_rad[1];

    // '<S96>:1:11' outBus.attCtrlInputs.ctrlInputsArray(3).cmd = rcOutCmds.yawStick; 
    rtb_VectorConcatenate[2].cmd = fcsModel_DW.rcOutCmds.yawStick;

    // '<S96>:1:12' outBus.attCtrlInputs.ctrlInputsArray(3).meas = stateEstimate.attitude_rad(3); 
    rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.attitude_rad[2];

    // '<S96>:1:14' outBus.attCtrlInputs.ctrlInputsArray(1).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[0].integratorReset = resetIntegrator;

    // '<S96>:1:15' outBus.attCtrlInputs.ctrlInputsArray(2).integratorReset = resetIntegrator; 
    rtb_VectorConcatenate[1].integratorReset = resetIntegrator;

    // '<S96>:1:16' outBus.attCtrlInputs.ctrlInputsArray(3).integratorReset = true; 
    rtb_VectorConcatenate[2].integratorReset = true;

    // RelationalOperator: '<S93>/Compare' incorporates:
    //   Constant: '<S93>/Constant'
    //   MATLAB Function: '<S4>/Interpret RC In Cmds'

    rtb_Compare_f = (flightMode != enumFlightMode::ALT_CONTROL);

    // Switch: '<S3>/Switch2' incorporates:
    //   Switch: '<S3>/Switch'

    if (rtb_Compare_d) {
      // Switch: '<S3>/Switch2' incorporates:
      //   BusAssignment: '<S155>/Bus Assignment'
      //   Concatenate: '<S155>/Vector Concatenate'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0] =
        rtb_BusAssignment1_a;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[1] = rtb_BusAssignment2;
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment4;
    } else if (rtb_Compare_f) {
      // Switch: '<S3>/Switch2' incorporates:
      //   MATLAB Function: '<S3>/assembleOuterLoopToInnerLoopBus'
      //   Switch: '<S3>/Switch'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N =
        fcsModel_DW.rcOutCmds.throttleStick;
      std::memcpy(&fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0],
                  &rtb_VectorConcatenate[0], 3U * sizeof(busCtrlInputs));
    } else {
      // Switch: '<S3>/Switch2' incorporates:
      //   BusAssignment: '<S155>/Bus Assignment'
      //   Concatenate: '<S155>/Vector Concatenate'
      //   Concatenate: '<S3>/Vector Concatenate'
      //   Switch: '<S3>/Switch'

      fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N = rtb_frcCmd_N;
      std::memcpy(&fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0],
                  &rtb_VectorConcatenate[0], sizeof(busCtrlInputs) << 1U);
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2] = rtb_BusAssignment4;
    }

    // End of Switch: '<S3>/Switch2'
  }

  // Outputs for Iterator SubSystem: '<S9>/Attitude Control' incorporates:
  //   ForEach: '<S51>/For Each'

  for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
    // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block'
    // ForEachSliceSelector generated from: '<S51>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'

    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].cmd,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.cmdSignalConditioningParamsArray
       [ForEach_itr_l], &tlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock);

    // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block'

    // Outputs for Atomic SubSystem: '<S51>/Signal Conditioning Block1'
    fcsMod_SignalConditioningBlock1
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].meas,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.measSignalConditioningParamsArray
       [ForEach_itr_l], &rlim, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S51>/Signal Conditioning Block1'

    // MATLAB Function: '<S51>/pickAttitudeCmdAndMeas' incorporates:
    //   Constant: '<S9>/Constant'
    //   ForEachSliceSelector generated from: '<S51>/index'

    tCmd = tlim;
    tCmd_unitRange = rlim;

    //  Passes cmd and meas as it is for roll and pitch channel
    //  but for yaw channel computes shortest angular distance between cmd Yaw
    //  and meas Yaw and overwrites Yaw cmd with that error and sets the meas Yaw to 
    //  zero for PID block
    // MATLAB Function 'Attitude Controller/Attitude Control/pickAttitudeCmdAndMeas': '<S56>:1' 
    // '<S56>:1:6' if index == cast(3, 'uint8')
    if (fcsModel_ConstP.Constant_Value_e[ForEach_itr_l] == 3) {
      // '<S56>:1:7' diff = mod(( cmd - meas + pi ), 2*pi) - pi;
      tCmd = (tlim - rlim) + 3.1415926535897931;
      if (tCmd == 0.0) {
        plim = 0.0;
      } else {
        plim = std::fmod(tCmd, 6.2831853071795862);
        resetIntegrator = (plim == 0.0);
        if (!resetIntegrator) {
          tCmd_unitRange = std::abs(tCmd / 6.2831853071795862);
          resetIntegrator = (std::abs(tCmd_unitRange - std::floor(tCmd_unitRange
            + 0.5)) <= 2.2204460492503131E-16 * tCmd_unitRange);
        }

        if (resetIntegrator) {
          plim = 0.0;
        } else if (tCmd < 0.0) {
          plim += 6.2831853071795862;
        }
      }

      tCmd = plim - 3.1415926535897931;

      // '<S56>:1:8' if diff < -pi
      if (plim - 3.1415926535897931 < -3.1415926535897931) {
        // '<S56>:1:9' diff = diff + 2*pi;
        tCmd = (plim - 3.1415926535897931) + 6.2831853071795862;
      }

      // '<S56>:1:12' cmd = diff;
      // '<S56>:1:13' meas = 0;
      tCmd_unitRange = 0.0;
    }

    // Outputs for Atomic SubSystem: '<S51>/pidWithDebug'
    // ForEachSliceSelector generated from: '<S51>/ctrlInputs' incorporates:
    //   Inport: '<Root>/ctrlParams'
    //   MATLAB Function: '<S51>/pickAttitudeCmdAndMeas'
    //   UnitDelay: '<S51>/Unit Delay'

    fcsModel_pidWithDebug
      (fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       feedForwardCmd, tCmd, tCmd_unitRange,
       fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[ForEach_itr_l].
       integratorReset, 0.0,
       &fcsModel_U.ctrlParams.innerLoopCtrlParams.attCtrlParams.ctrlParamsArray[ForEach_itr_l],
       fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE, &plim,
       &rtb_BusCreator_og, 0.004, &fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
       pidWithDebug);

    // End of Outputs for SubSystem: '<S51>/pidWithDebug'

    // Update for UnitDelay: '<S51>/Unit Delay'
    fcsModel_DW.CoreSubsys_p[ForEach_itr_l].UnitDelay_DSTATE = plim;

    // ForEachSliceAssignment generated from: '<S51>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[ForEach_itr_l] =
      rtb_BusCreator_og;

    // ForEachSliceAssignment generated from: '<S51>/angRateCmd '
    rtb_ImpAsg_InsertedFor_angRateC[ForEach_itr_l] = plim;

    // ForEachSliceAssignment generated from: '<S51>/measFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[ForEach_itr_l] =
      rlim;

    // ForEachSliceAssignment generated from: '<S51>/cmdFlt'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[ForEach_itr_l] =
      tlim;
  }

  // End of Outputs for SubSystem: '<S9>/Attitude Control'

  // Switch: '<S8>/Switch'
  rtb_ImpAsg_InsertedFor_velCtrlO[0] = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_ImpAsg_InsertedFor_velCtrlO[1] = rtb_ImpAsg_InsertedFor_angRateC[1];

  // Switch: '<S9>/Switch' incorporates:
  //   Constant: '<S52>/Constant'
  //   MATLAB Function: '<S4>/Interpret RC In Cmds'
  //   RelationalOperator: '<S52>/Compare'
  //   Switch: '<S8>/Switch'

  if (flightMode == enumFlightMode::STABILIZE) {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] =
      fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[2].cmd;
  } else {
    rtb_ImpAsg_InsertedFor_velCtrlO[2] = rtb_ImpAsg_InsertedFor_angRateC[2];
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
  tlim = fcsModel_U.stateEstimate.attitude_rad[1];

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
        tlim = -1.0;
      } else {
        tlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      tlim *= rlim - 0.004245395477824045;
    } else if (std::abs(fcsModel_U.stateEstimate.attitude_rad[1]) -
               4.71238898038469 <= 0.0) {
      // 'eulerRates2bodyRates_function:26' pitch = sign(pitch)*( abs(pitch) - limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        tlim = -1.0;
      } else {
        tlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      tlim *= rlim - 0.004245395477824045;
    } else {
      // 'eulerRates2bodyRates_function:27' else
      // 'eulerRates2bodyRates_function:28' pitch = sign(pitch)*( abs(pitch) + limit); 
      if (fcsModel_U.stateEstimate.attitude_rad[1] < 0.0) {
        tlim = -1.0;
      } else {
        tlim = (fcsModel_U.stateEstimate.attitude_rad[1] > 0.0);
      }

      tlim *= rlim + 0.004245395477824045;
    }
  }

  // Construct conversion matrix
  // 'eulerRates2bodyRates_function:33' conversionMatrix = [1, 0, -sin(pitch);
  // 'eulerRates2bodyRates_function:34'     0, cos(roll), sin(roll)*cos(pitch);
  // 'eulerRates2bodyRates_function:35'     0, -sin(roll), cos(roll)*cos(pitch)]; 
  rlim = std::sin(fcsModel_U.stateEstimate.attitude_rad[0]);
  plim = std::cos(fcsModel_U.stateEstimate.attitude_rad[0]);
  tCmd = std::cos(tlim);
  rtb_Transpose[0] = 1.0;
  rtb_Transpose[3] = 0.0;
  rtb_Transpose[6] = -std::sin(tlim);
  rtb_Transpose[1] = 0.0;
  rtb_Transpose[4] = plim;
  rtb_Transpose[7] = rlim * tCmd;
  rtb_Transpose[2] = 0.0;
  rtb_Transpose[5] = -rlim;
  rtb_Transpose[8] = plim * tCmd;

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
  for (ii = 0; ii < 3; ii++) {
    // 'zeroSmallValues:11' for jj = 1:size(M,2)
    tlim = rtb_Transpose[ii];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (tlim <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      tlim = 0.0;
    }

    rtb_Transpose[ii] = tlim;
    rlim = tlim * rtb_ImpAsg_InsertedFor_velCtrlO[0];
    tlim = rtb_Transpose[ii + 3];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(tlim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      tlim = 0.0;
    }

    rtb_Transpose[ii + 3] = tlim;
    rlim += tlim * rtb_ImpAsg_InsertedFor_velCtrlO[1];
    tlim = rtb_Transpose[ii + 6];

    // 'zeroSmallValues:12' if(abs(M(ii,jj))<= abs(eps))
    if (std::abs(tlim) <= 1.0E-12) {
      // 'zeroSmallValues:13' M(ii,jj) = 0;
      tlim = 0.0;
    }

    rtb_Transpose[ii + 6] = tlim;
    rtb_ImpAsg_InsertedFor_angRateC[ii] = tlim *
      rtb_ImpAsg_InsertedFor_velCtrlO[2] + rlim;
  }

  // End of MATLAB Function: '<S8>/EulerRates2BodyRates'

  // BusCreator: '<S8>/Bus Creator' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[0].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[0].cmd = rtb_ImpAsg_InsertedFor_angRateC[0];
  rtb_VectorConcatenate[0].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[0];
  rtb_VectorConcatenate[0].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[0].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator3' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[1].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[1].cmd = rtb_ImpAsg_InsertedFor_angRateC[1];
  rtb_VectorConcatenate[1].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[1];
  rtb_VectorConcatenate[1].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
  rtb_VectorConcatenate[1].trackingCtrlCmd = 0.0;

  // BusCreator: '<S8>/Bus Creator4' incorporates:
  //   Concatenate: '<S8>/Vector Concatenate'
  //   Inport: '<Root>/stateEstimate'

  rtb_VectorConcatenate[2].feedForwardCmd = 0.0;
  rtb_VectorConcatenate[2].cmd = rtb_ImpAsg_InsertedFor_angRateC[2];
  rtb_VectorConcatenate[2].meas = fcsModel_U.stateEstimate.bodyAngRates_radps[2];
  rtb_VectorConcatenate[2].integratorReset =
    fcsModel_DW.Switch2.attCtrlInputs.ctrlInputsArray[0].integratorReset;
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
      [ForEach_itr_p], &tlim, 0.004, &fcsModel_DW.CoreSubsys[ForEach_itr_p].
      SignalConditioningBlock1);

    // End of Outputs for SubSystem: '<S10>/Signal Conditioning Block1'

    // Outputs for Atomic SubSystem: '<S10>/pidWithDebug'
    fcsModel_pidWithDebug(0.0, plim, tlim, rtb_VectorConcatenate[ForEach_itr_p].
                          integratorReset, 0.0,
                          &fcsModel_U.ctrlParams.innerLoopCtrlParams.angRateCtrlParams.ctrlParamsArray
                          [ForEach_itr_p], fcsModel_DW.CoreSubsys[ForEach_itr_p]
                          .UnitDelay_DSTATE, &rlim, &rtb_BusCreator_og, 0.004,
                          &fcsModel_DW.CoreSubsys[ForEach_itr_p].pidWithDebug);

    // End of Outputs for SubSystem: '<S10>/pidWithDebug'

    // Update for UnitDelay: '<S10>/Unit Delay'
    fcsModel_DW.CoreSubsys[ForEach_itr_p].UnitDelay_DSTATE = rlim;

    // ForEachSliceAssignment generated from: '<S10>/pidDebug'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[ForEach_itr_p]
      = rtb_BusCreator_og;

    // ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'
    rtb_ImpAsg_InsertedFor_angAccel[ForEach_itr_p] = rlim;

    // ForEachSliceAssignment generated from: '<S10>/filtMeas'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[ForEach_itr_p] =
      tlim;

    // ForEachSliceAssignment generated from: '<S10>/filtCmd'
    fcsModel_Y.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[ForEach_itr_p] =
      plim;
  }

  // End of Outputs for SubSystem: '<S7>/For Each Subsystem'
  // End of Outputs for SubSystem: '<S2>/Angular Rate Controller'

  // Product: '<S2>/Matrix Multiply' incorporates:
  //   Constant: '<S2>/Constant'
  //   ForEachSliceAssignment generated from: '<S10>/angAccelCmd_radps2'

  for (ii = 0; ii < 3; ii++) {
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

  rlim = rtb_ImpAsg_InsertedFor_angRateC[0];
  plim = rtb_ImpAsg_InsertedFor_angRateC[1];
  tCmd = rtb_ImpAsg_InsertedFor_angRateC[2];

  // RelationalOperator: '<S6>/Compare' incorporates:
  //   Constant: '<S6>/Constant'

  // Unit Conversion - from: rad/s to: rpm
  // Expression: output = (9.5493*input) + (0)
  resetIntegrator = (state == enumStateMachine::INACTIVE);
  for (ii = 0; ii < 4; ii++) {
    // Product: '<S1>/Matrix Multiply' incorporates:
    //   BusCreator: '<S2>/Bus Creator1'
    //   Constant: '<S1>/Constant'

    tlim = ((fcsModel_ConstP.Constant_Value_c[ii + 4] * rlim +
             fcsModel_ConstP.Constant_Value_c[ii] *
             fcsModel_DW.Switch2.outerLoopCmds.thrustCmd_N) +
            fcsModel_ConstP.Constant_Value_c[ii + 8] * plim) +
      fcsModel_ConstP.Constant_Value_c[ii + 12] * tCmd;

    // Saturate: '<S1>/Saturation'
    if (tlim > 839601.76328711538) {
      tlim = 839601.76328711538;
    } else if (tlim < 0.0) {
      tlim = 0.0;
    }

    // End of Saturate: '<S1>/Saturation'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn' incorporates:
    //   Sqrt: '<S1>/Sqrt'
    //   UnitConversion: '<S5>/Unit Conversion'

    tlim = 9.5492965855137211 * std::sqrt(tlim) - -0.84897325901283394 *
      fcsModel_DW.DiscreteTransferFcn_states[ii];

    // Switch: '<S1>/Switch'
    if (resetIntegrator) {
      // Outport: '<Root>/actuatorsCmds'
      fcsModel_Y.actuatorsCmds[ii] = -1.0;
    } else {
      // Outport: '<Root>/actuatorsCmds' incorporates:
      //   DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'

      fcsModel_Y.actuatorsCmds[ii] = 0.075513370493583074 * tlim +
        0.075513370493583074 * fcsModel_DW.DiscreteTransferFcn_states[ii];
    }

    // End of Switch: '<S1>/Switch'

    // DiscreteTransferFcn: '<S1>/Discrete Transfer Fcn'
    DiscreteTransferFcn_tmp[ii] = tlim;
  }

  // RateTransition: '<Root>/Rate Transition'
  if ((&fcsModel_M)->Timing.TaskCounters.TID[1] == 0) {
    fcsModel_Y.fcsDebug.outerLoopCtrlDebug = fcsModel_DW.RateTransition_Buffer0;

    // Switch: '<S3>/Switch3' incorporates:
    //   RateTransition: '<Root>/Rate Transition'
    //   Switch: '<S3>/Switch1'

    if (rtb_Compare_d) {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S155>/Bus Assignment'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    } else if (rtb_Compare_f) {
      // Switch: '<S3>/Switch1' incorporates:
      //   BusCreator: '<S3>/Bus Creator'

      rtb_BusCreator_b_frcCmd_N = rtb_Product6;
    } else {
      // BusCreator: '<S3>/Bus Creator' incorporates:
      //   BusAssignment: '<S155>/Bus Assignment'
      //   Switch: '<S3>/Switch1'

      rtb_BusCreator_b_frcCmd_N = rtb_frcCmd_N;
    }

    // End of Switch: '<S3>/Switch3'

    // BusCreator: '<S3>/Bus Creator' incorporates:
    //   BusCreator: '<S98>/Bus Creator'
    //   BusCreator: '<S99>/Bus Creator'
    //   ForEachSliceAssignment generated from: '<S103>/cmd'
    //   ForEachSliceAssignment generated from: '<S156>/filtCmd'

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
    rtb_BusCreator_b_zAccelCtrlDebu = fcsModel_rtZbusZaccelCtrlDebug;
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
    fcsModel_DW.RateTransition_Buffer0.zAccelCtrlDebug =
      rtb_BusCreator_b_zAccelCtrlDebu;
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

    // 'interpretRcInputs_function:24' throttle_is_up = false;
    // 'holdOutputAtCenter_function:6' last_input = 0;
    // SystemInitialize for Iterator SubSystem: '<S98>/NED Position Control'
    for (ForEach_itr_i = 0; ForEach_itr_i < 3; ForEach_itr_i++) {
      // SystemInitialize for Iterator SubSystem: '<S98>/NED Position Control'
      // SystemInitialize for Atomic SubSystem: '<S103>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S103>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S103>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S103>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S103>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_g[ForEach_itr_i].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S103>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S98>/NED Position Control'
    }

    // End of SystemInitialize for SubSystem: '<S98>/NED Position Control'
    // SystemInitialize for Iterator SubSystem: '<S99>/For Each Subsystem'
    for (ForEach_itr = 0; ForEach_itr < 3; ForEach_itr++) {
      // SystemInitialize for Iterator SubSystem: '<S99>/For Each Subsystem'
      // SystemInitialize for Atomic SubSystem: '<S156>/Signal Conditioning Block' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S156>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S156>/Signal Conditioning Block1' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S156>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S156>/Signal Conditioning Block2' 
      SignalConditioningBlock1_c_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        SignalConditioningBlock2);

      // End of SystemInitialize for SubSystem: '<S156>/Signal Conditioning Block2' 

      // SystemInitialize for Atomic SubSystem: '<S156>/pidWithDebug'
      fcsModel_pidWithDebug_m_Init(&fcsModel_DW.CoreSubsys_i[ForEach_itr].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S156>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S99>/For Each Subsystem'
    }

    // End of SystemInitialize for SubSystem: '<S99>/For Each Subsystem'
    // 'holdOutputAtCenter_function:6' last_input = 0;
    // SystemInitialize for Iterator SubSystem: '<S9>/Attitude Control'
    for (ForEach_itr_l = 0; ForEach_itr_l < 3; ForEach_itr_l++) {
      // SystemInitialize for Iterator SubSystem: '<S9>/Attitude Control'
      // SystemInitialize for Atomic SubSystem: '<S51>/Signal Conditioning Block' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock);

      // End of SystemInitialize for SubSystem: '<S51>/Signal Conditioning Block' 

      // SystemInitialize for Atomic SubSystem: '<S51>/Signal Conditioning Block1' 
      f_SignalConditioningBlock1_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        SignalConditioningBlock1);

      // End of SystemInitialize for SubSystem: '<S51>/Signal Conditioning Block1' 

      // SystemInitialize for Atomic SubSystem: '<S51>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys_p[ForEach_itr_l].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S51>/pidWithDebug'
      // End of SystemInitialize for SubSystem: '<S9>/Attitude Control'
    }

    // End of SystemInitialize for SubSystem: '<S9>/Attitude Control'
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

      // SystemInitialize for Atomic SubSystem: '<S10>/pidWithDebug'
      fcsModel_pidWithDebug_Init(&fcsModel_DW.CoreSubsys[ForEach_itr_p].
        pidWithDebug);

      // End of SystemInitialize for SubSystem: '<S10>/pidWithDebug'
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
