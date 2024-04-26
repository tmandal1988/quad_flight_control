//
// File: fcsModel.h
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.101
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Wed Apr 24 14:47:06 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: ARM Compatible->ARM 7
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
//    3. Traceability
// Validation result: All passed
//
#ifndef RTW_HEADER_fcsModel_h_
#define RTW_HEADER_fcsModel_h_
#include "rtwtypes.h"
#include "fcsModel_types.h"
#include <array>

// External data declarations for dependent source files
extern const busFcsDebug fcsModel_rtZbusFcsDebug;// busFcsDebug ground

// Class declaration for model fcsModel
class fcsModel final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for system '<S1>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_T {
    uint32_T Prelookup_DWORK1;         // '<S7>/Prelookup'
  };

  // Block signals and states (default storage) for system '<S14>/Discrete First Order Deriv Filter' 
  struct DW_DiscreteFirstOrderDerivFil_T {
    std::array<real_T, 2> num;
                      // '<S45>/Compute Deriv Filter Numerator And Denominator'
    std::array<real_T, 2> den;
                      // '<S45>/Compute Deriv Filter Numerator And Denominator'
    real_T DiscreteTransferFcn_states; // '<S45>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S11>/pidWithDebug'
  struct DW_pidWithDebug_fcsModel_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                   // '<S14>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S14>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S46>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S14>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S14>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S14>/Discrete-Time Integrator'
    uint8_T DiscreteTimeIntegrator_IC_LOADI;// '<S14>/Discrete-Time Integrator'
  };

  // Block signals and states (default storage) for system '<S11>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_f_T {
    std::array<real_T, 3> num;
                            // '<S31>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                            // '<S31>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S31>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S32>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S31>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S8>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_c_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S11>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S11>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S11>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S11>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S10>/Attitude Control' 
  struct DW_CoreSubsys_fcsModel_i_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S52>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S52>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S52>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S52>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S103>/holdOutputAtCenter1' 
  struct DW_holdOutputAtCenter1_fcsMod_T {
    real_T last_input;                 // '<S113>/holdOutputAtCenter'
  };

  // Block signals and states (default storage) for system '<S104>/pidWithDebug' 
  struct DW_pidWithDebug_fcsModel_i_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                  // '<S120>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S120>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S152>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S120>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S120>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S120>/Discrete-Time Integrator'
    uint8_T DiscreteTimeIntegrator_IC_LOADI;// '<S120>/Discrete-Time Integrator' 
  };

  // Block signals and states (default storage) for system '<S104>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_g_T {
    std::array<real_T, 3> num;
                           // '<S137>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                           // '<S137>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S137>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S138>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S137>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S99>/NED Position Control' 
  struct DW_CoreSubsys_fcsModel_b_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S104>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S104>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S104>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S104>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S100>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_p_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S157>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S157>/Signal Conditioning Block1' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock2;// '<S157>/Signal Conditioning Block2' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S157>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S157>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW_fcsModel_T {
    std::array<DW_CoreSubsys_fcsModel_p_T, 3> CoreSubsys_i;// '<S100>/For Each Subsystem' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S169>/pidWithDebug'
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S169>/Signal Conditioning Block' 
    std::array<DW_CoreSubsys_fcsModel_b_T, 3> CoreSubsys_g;// '<S99>/NED Position Control' 
    DW_holdOutputAtCenter1_fcsMod_T holdOutputAtCenter2;// '<S103>/holdOutputAtCenter2' 
    DW_holdOutputAtCenter1_fcsMod_T holdOutputAtCenter1;// '<S103>/holdOutputAtCenter1' 
    std::array<DW_CoreSubsys_fcsModel_i_T, 3> CoreSubsys_p;// '<S10>/Attitude Control' 
    std::array<DW_CoreSubsys_fcsModel_c_T, 3> CoreSubsys_a;// '<S8>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_T, 4> CoreSubsys;// '<S1>/For Each Subsystem' 
    busOuterLoopCtrlDebug RateTransition_Buffer0;// '<Root>/Rate Transition'
    busOuterLoopToInnerLoop Switch2;   // '<S3>/Switch2'
    busRcOutCmds rcOutCmds;            // '<S4>/Interpret RC In Cmds'
    std::array<real_T, 4> DiscreteTransferFcn_states;// '<S1>/Discrete Transfer Fcn' 
    std::array<real_T, 2> DiscreteTransferFcn_states_n;// '<S188>/Discrete Transfer Fcn' 
    std::array<real_T, 2> DiscreteTransferFcn_states_nh;// '<S189>/Discrete Transfer Fcn' 
    real_T UnitDelay_DSTATE;           // '<S169>/Unit Delay'
    real_T DelayInput2_DSTATE;         // '<S191>/Delay Input2'
    real_T DelayInput2_DSTATE_e;       // '<S190>/Delay Input2'
    real_T last_input;                 // '<S162>/holdOutputAtCenter'
    real_T last_input_c;               // '<S112>/holdOutputAtCenter'
    int32_T durationCounter_1;         // '<S4>/Chart'
    int32_T durationCounter_1_j;       // '<S4>/Chart'
    uint16_T temporalCounter_i1;       // '<S4>/Chart'
    uint8_T is_active_c1_rcInterpreter;// '<S4>/Chart'
    uint8_T is_c1_rcInterpreter;       // '<S4>/Chart'
    boolean_T throttle_is_up;          // '<S4>/Interpret RC In Cmds'
    boolean_T rcCheckFlag;             // '<S4>/Chart'
  };

  // Constant parameters (default storage)
  struct ConstP_fcsModel_T {
    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S3>/Constant'
    //    '<S156>/Constant'

    busOuterLoopToInnerLoop pooled3;

    // Expression: rpmToPwmLut(:, 1)
    //  Referenced by: '<S7>/Prelookup'

    std::array<real_T, 15> Prelookup_BreakpointsData;

    // Expression: rpmToPwmLut(:, 2)
    //  Referenced by: '<S7>/Interpolation Using Prelookup'

    std::array<real_T, 15> InterpolationUsingPrelookup_Tab;

    // Expression: allocationDataStruct.allocationMatrix
    //  Referenced by: '<S1>/Constant'

    std::array<real_T, 16> Constant_Value_c;

    // Expression: vehicleConstants.inertia_kgm2
    //  Referenced by: '<S2>/Constant'

    std::array<real_T, 9> Constant_Value_n;

    // Computed Parameter: Constant_Value_e
    //  Referenced by: '<S10>/Constant'

    std::array<uint8_T, 3> Constant_Value_e;
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_fcsModel_T {
    busRcInCmds rcCmdsIn;              // '<Root>/rcCmdsIn'
    busStateEstimate stateEstimate;    // '<Root>/stateEstimate'
    busFcsParams ctrlParams;           // '<Root>/ctrlParams'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_fcsModel_T {
    std::array<real_T, 4> actuatorsCmds;// '<Root>/actuatorsCmds'
    busFcsDebug fcsDebug;              // '<Root>/fcsDebug'
    std::array<real_T, 4> actuatorsPwmCmds;// '<Root>/actuatorsPwmCmds'
  };

  // Real-time Model Data Structure
  struct RT_MODEL_fcsModel_T {
    //
    //  Timing:
    //  The following substructure contains information regarding
    //  the timing information for the model.

    struct {
      struct {
        uint8_T TID[2];
      } TaskCounters;
    } Timing;
  };

  // Copy Constructor
  fcsModel(fcsModel const&) = delete;

  // Assignment Operator
  fcsModel& operator= (fcsModel const&) & = delete;

  // Move Constructor
  fcsModel(fcsModel &&) = delete;

  // Move Assignment Operator
  fcsModel& operator= (fcsModel &&) = delete;

  // Real-Time Model get method
  fcsModel::RT_MODEL_fcsModel_T * getRTM();

  // Root inports set method
  void setExternalInputs(const ExtU_fcsModel_T *pExtU_fcsModel_T)
  {
    fcsModel_U = *pExtU_fcsModel_T;
  }

  // Root outports get method
  const ExtY_fcsModel_T &getExternalOutputs() const
  {
    return fcsModel_Y;
  }

  // model initialize function
  void initialize();
  void ModelExternalOutputInit();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  fcsModel();

  // Destructor
  ~fcsModel();

  // private data and function members
 private:
  // External inputs
  ExtU_fcsModel_T fcsModel_U;

  // External outputs
  ExtY_fcsModel_T fcsModel_Y;

  // Block states
  DW_fcsModel_T fcsModel_DW;

  // private member function(s) for subsystem '<S14>/Discrete First Order Deriv Filter'
  static void f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
    rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
    rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW);

  // private member function(s) for subsystem '<S11>/pidWithDebug'
  static void fcsModel_pidWithDebug_Init(DW_pidWithDebug_fcsModel_T *localDW);
  static void fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc,
    const busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T
    *rty_ctrlCmd, busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_T *localDW);

  // private member function(s) for subsystem '<S30>/Compute Natural Frequency'
  static void fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps, real_T
    rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps);

  // private member function(s) for subsystem '<S30>/Compute Numerator And Denominator'
  static void ComputeNumeratorAndDenominator(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_rateNum[3], real_T rty_accelNum[3],
    real_T rty_den[3], real_T rtp_sampleTime_s);

  // private member function(s) for subsystem '<S31>/Compute Filter Numerator And Denominator'
  static void ComputeFilterNumeratorAndD_Init(real_T rty_num[3], real_T rty_den
    [3]);
  static void ComputeFilterNumeratorAndDenomi(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_num[3], real_T rty_den[3], real_T
    rtp_sampleTime_s);

  // private member function(s) for subsystem '<S11>/Signal Conditioning Block1'
  static void f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
    *localDW);
  static void fcsMod_SignalConditioningBlock1(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW);

  // private member function(s) for subsystem '<S103>/holdOutputAtCenter1'
  static void fcsModel_holdOutputAtCenter1(real_T rtu_input, real_T rtu_trigger,
    real_T *rty_output, boolean_T *rty_atCenter, DW_holdOutputAtCenter1_fcsMod_T
    *localDW);

  // private member function(s) for subsystem '<S104>/pidWithDebug'
  static void fcsModel_pidWithDebug_m_Init(DW_pidWithDebug_fcsModel_i_T *localDW);
  static void fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc,
    const busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T
    *rty_ctrlCmd, busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_i_T *localDW);

  // private member function(s) for subsystem '<S104>/Signal Conditioning Block1'
  static void SignalConditioningBlock1_c_Init(DW_SignalConditioningBlock1_g_T
    *localDW);
  static void fcsM_SignalConditioningBlock1_f(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_g_T *localDW);

  // private member function(s) for subsystem '<Root>'
  boolean_T fcsModel_checkRcCmds(const busRcInCmds
    *BusConversion_InsertedFor_Chart);

  // Real-Time Model
  RT_MODEL_fcsModel_T fcsModel_M;
};

// Constant parameters (default storage)
extern const fcsModel::ConstP_fcsModel_T fcsModel_ConstP;

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S15>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S15>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S17>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S27>/Data Type Duplicate' : Unused code path elimination
//  Block '<S27>/Data Type Propagation' : Unused code path elimination
//  Block '<S18>/Delay Input2' : Unused code path elimination
//  Block '<S18>/Difference Inputs1' : Unused code path elimination
//  Block '<S18>/Difference Inputs2' : Unused code path elimination
//  Block '<S18>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S28>/Data Type Duplicate' : Unused code path elimination
//  Block '<S28>/Data Type Propagation' : Unused code path elimination
//  Block '<S28>/LowerRelop1' : Unused code path elimination
//  Block '<S28>/Switch' : Unused code path elimination
//  Block '<S28>/Switch2' : Unused code path elimination
//  Block '<S28>/UpperRelop' : Unused code path elimination
//  Block '<S18>/delta fall limit' : Unused code path elimination
//  Block '<S18>/delta rise limit' : Unused code path elimination
//  Block '<S18>/sample time' : Unused code path elimination
//  Block '<S19>/Delay Input2' : Unused code path elimination
//  Block '<S19>/Difference Inputs1' : Unused code path elimination
//  Block '<S19>/Difference Inputs2' : Unused code path elimination
//  Block '<S19>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S29>/Data Type Duplicate' : Unused code path elimination
//  Block '<S29>/Data Type Propagation' : Unused code path elimination
//  Block '<S29>/LowerRelop1' : Unused code path elimination
//  Block '<S29>/Switch' : Unused code path elimination
//  Block '<S29>/Switch2' : Unused code path elimination
//  Block '<S29>/UpperRelop' : Unused code path elimination
//  Block '<S19>/delta fall limit' : Unused code path elimination
//  Block '<S19>/delta rise limit' : Unused code path elimination
//  Block '<S19>/sample time' : Unused code path elimination
//  Block '<S20>/Data Type Duplicate' : Unused code path elimination
//  Block '<S20>/Data Type Propagation' : Unused code path elimination
//  Block '<S21>/Data Type Duplicate' : Unused code path elimination
//  Block '<S21>/Data Type Propagation' : Unused code path elimination
//  Block '<S21>/LowerRelop1' : Unused code path elimination
//  Block '<S21>/Switch' : Unused code path elimination
//  Block '<S21>/Switch2' : Unused code path elimination
//  Block '<S21>/UpperRelop' : Unused code path elimination
//  Block '<S22>/Data Type Duplicate' : Unused code path elimination
//  Block '<S22>/Data Type Propagation' : Unused code path elimination
//  Block '<S22>/LowerRelop1' : Unused code path elimination
//  Block '<S22>/Switch' : Unused code path elimination
//  Block '<S22>/Switch2' : Unused code path elimination
//  Block '<S22>/UpperRelop' : Unused code path elimination
//  Block '<S30>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S30>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S32>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S42>/Data Type Duplicate' : Unused code path elimination
//  Block '<S42>/Data Type Propagation' : Unused code path elimination
//  Block '<S33>/Delay Input2' : Unused code path elimination
//  Block '<S33>/Difference Inputs1' : Unused code path elimination
//  Block '<S33>/Difference Inputs2' : Unused code path elimination
//  Block '<S33>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S43>/Data Type Duplicate' : Unused code path elimination
//  Block '<S43>/Data Type Propagation' : Unused code path elimination
//  Block '<S43>/LowerRelop1' : Unused code path elimination
//  Block '<S43>/Switch' : Unused code path elimination
//  Block '<S43>/Switch2' : Unused code path elimination
//  Block '<S43>/UpperRelop' : Unused code path elimination
//  Block '<S33>/delta fall limit' : Unused code path elimination
//  Block '<S33>/delta rise limit' : Unused code path elimination
//  Block '<S33>/sample time' : Unused code path elimination
//  Block '<S34>/Delay Input2' : Unused code path elimination
//  Block '<S34>/Difference Inputs1' : Unused code path elimination
//  Block '<S34>/Difference Inputs2' : Unused code path elimination
//  Block '<S34>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S44>/Data Type Duplicate' : Unused code path elimination
//  Block '<S44>/Data Type Propagation' : Unused code path elimination
//  Block '<S44>/LowerRelop1' : Unused code path elimination
//  Block '<S44>/Switch' : Unused code path elimination
//  Block '<S44>/Switch2' : Unused code path elimination
//  Block '<S44>/UpperRelop' : Unused code path elimination
//  Block '<S34>/delta fall limit' : Unused code path elimination
//  Block '<S34>/delta rise limit' : Unused code path elimination
//  Block '<S34>/sample time' : Unused code path elimination
//  Block '<S35>/Data Type Duplicate' : Unused code path elimination
//  Block '<S35>/Data Type Propagation' : Unused code path elimination
//  Block '<S36>/Data Type Duplicate' : Unused code path elimination
//  Block '<S36>/Data Type Propagation' : Unused code path elimination
//  Block '<S36>/LowerRelop1' : Unused code path elimination
//  Block '<S36>/Switch' : Unused code path elimination
//  Block '<S36>/Switch2' : Unused code path elimination
//  Block '<S36>/UpperRelop' : Unused code path elimination
//  Block '<S37>/Data Type Duplicate' : Unused code path elimination
//  Block '<S37>/Data Type Propagation' : Unused code path elimination
//  Block '<S37>/LowerRelop1' : Unused code path elimination
//  Block '<S37>/Switch' : Unused code path elimination
//  Block '<S37>/Switch2' : Unused code path elimination
//  Block '<S37>/UpperRelop' : Unused code path elimination
//  Block '<S46>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S49>/Data Type Duplicate' : Unused code path elimination
//  Block '<S49>/Data Type Propagation' : Unused code path elimination
//  Block '<S47>/Data Type Duplicate' : Unused code path elimination
//  Block '<S47>/Data Type Propagation' : Unused code path elimination
//  Block '<S59>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S59>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S61>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S71>/Data Type Duplicate' : Unused code path elimination
//  Block '<S71>/Data Type Propagation' : Unused code path elimination
//  Block '<S62>/Delay Input2' : Unused code path elimination
//  Block '<S62>/Difference Inputs1' : Unused code path elimination
//  Block '<S62>/Difference Inputs2' : Unused code path elimination
//  Block '<S62>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S72>/Data Type Duplicate' : Unused code path elimination
//  Block '<S72>/Data Type Propagation' : Unused code path elimination
//  Block '<S72>/LowerRelop1' : Unused code path elimination
//  Block '<S72>/Switch' : Unused code path elimination
//  Block '<S72>/Switch2' : Unused code path elimination
//  Block '<S72>/UpperRelop' : Unused code path elimination
//  Block '<S62>/delta fall limit' : Unused code path elimination
//  Block '<S62>/delta rise limit' : Unused code path elimination
//  Block '<S62>/sample time' : Unused code path elimination
//  Block '<S63>/Delay Input2' : Unused code path elimination
//  Block '<S63>/Difference Inputs1' : Unused code path elimination
//  Block '<S63>/Difference Inputs2' : Unused code path elimination
//  Block '<S63>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S73>/Data Type Duplicate' : Unused code path elimination
//  Block '<S73>/Data Type Propagation' : Unused code path elimination
//  Block '<S73>/LowerRelop1' : Unused code path elimination
//  Block '<S73>/Switch' : Unused code path elimination
//  Block '<S73>/Switch2' : Unused code path elimination
//  Block '<S73>/UpperRelop' : Unused code path elimination
//  Block '<S63>/delta fall limit' : Unused code path elimination
//  Block '<S63>/delta rise limit' : Unused code path elimination
//  Block '<S63>/sample time' : Unused code path elimination
//  Block '<S64>/Data Type Duplicate' : Unused code path elimination
//  Block '<S64>/Data Type Propagation' : Unused code path elimination
//  Block '<S65>/Data Type Duplicate' : Unused code path elimination
//  Block '<S65>/Data Type Propagation' : Unused code path elimination
//  Block '<S65>/LowerRelop1' : Unused code path elimination
//  Block '<S65>/Switch' : Unused code path elimination
//  Block '<S65>/Switch2' : Unused code path elimination
//  Block '<S65>/UpperRelop' : Unused code path elimination
//  Block '<S66>/Data Type Duplicate' : Unused code path elimination
//  Block '<S66>/Data Type Propagation' : Unused code path elimination
//  Block '<S66>/LowerRelop1' : Unused code path elimination
//  Block '<S66>/Switch' : Unused code path elimination
//  Block '<S66>/Switch2' : Unused code path elimination
//  Block '<S66>/UpperRelop' : Unused code path elimination
//  Block '<S74>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S74>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S76>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S86>/Data Type Duplicate' : Unused code path elimination
//  Block '<S86>/Data Type Propagation' : Unused code path elimination
//  Block '<S77>/Delay Input2' : Unused code path elimination
//  Block '<S77>/Difference Inputs1' : Unused code path elimination
//  Block '<S77>/Difference Inputs2' : Unused code path elimination
//  Block '<S77>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S87>/Data Type Duplicate' : Unused code path elimination
//  Block '<S87>/Data Type Propagation' : Unused code path elimination
//  Block '<S87>/LowerRelop1' : Unused code path elimination
//  Block '<S87>/Switch' : Unused code path elimination
//  Block '<S87>/Switch2' : Unused code path elimination
//  Block '<S87>/UpperRelop' : Unused code path elimination
//  Block '<S77>/delta fall limit' : Unused code path elimination
//  Block '<S77>/delta rise limit' : Unused code path elimination
//  Block '<S77>/sample time' : Unused code path elimination
//  Block '<S78>/Delay Input2' : Unused code path elimination
//  Block '<S78>/Difference Inputs1' : Unused code path elimination
//  Block '<S78>/Difference Inputs2' : Unused code path elimination
//  Block '<S78>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S88>/Data Type Duplicate' : Unused code path elimination
//  Block '<S88>/Data Type Propagation' : Unused code path elimination
//  Block '<S88>/LowerRelop1' : Unused code path elimination
//  Block '<S88>/Switch' : Unused code path elimination
//  Block '<S88>/Switch2' : Unused code path elimination
//  Block '<S88>/UpperRelop' : Unused code path elimination
//  Block '<S78>/delta fall limit' : Unused code path elimination
//  Block '<S78>/delta rise limit' : Unused code path elimination
//  Block '<S78>/sample time' : Unused code path elimination
//  Block '<S79>/Data Type Duplicate' : Unused code path elimination
//  Block '<S79>/Data Type Propagation' : Unused code path elimination
//  Block '<S80>/Data Type Duplicate' : Unused code path elimination
//  Block '<S80>/Data Type Propagation' : Unused code path elimination
//  Block '<S80>/LowerRelop1' : Unused code path elimination
//  Block '<S80>/Switch' : Unused code path elimination
//  Block '<S80>/Switch2' : Unused code path elimination
//  Block '<S80>/UpperRelop' : Unused code path elimination
//  Block '<S81>/Data Type Duplicate' : Unused code path elimination
//  Block '<S81>/Data Type Propagation' : Unused code path elimination
//  Block '<S81>/LowerRelop1' : Unused code path elimination
//  Block '<S81>/Switch' : Unused code path elimination
//  Block '<S81>/Switch2' : Unused code path elimination
//  Block '<S81>/UpperRelop' : Unused code path elimination
//  Block '<S90>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S93>/Data Type Duplicate' : Unused code path elimination
//  Block '<S93>/Data Type Propagation' : Unused code path elimination
//  Block '<S91>/Data Type Duplicate' : Unused code path elimination
//  Block '<S91>/Data Type Propagation' : Unused code path elimination
//  Block '<S121>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S121>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S123>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S133>/Data Type Duplicate' : Unused code path elimination
//  Block '<S133>/Data Type Propagation' : Unused code path elimination
//  Block '<S124>/Delay Input2' : Unused code path elimination
//  Block '<S124>/Difference Inputs1' : Unused code path elimination
//  Block '<S124>/Difference Inputs2' : Unused code path elimination
//  Block '<S124>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S134>/Data Type Duplicate' : Unused code path elimination
//  Block '<S134>/Data Type Propagation' : Unused code path elimination
//  Block '<S134>/LowerRelop1' : Unused code path elimination
//  Block '<S134>/Switch' : Unused code path elimination
//  Block '<S134>/Switch2' : Unused code path elimination
//  Block '<S134>/UpperRelop' : Unused code path elimination
//  Block '<S124>/delta fall limit' : Unused code path elimination
//  Block '<S124>/delta rise limit' : Unused code path elimination
//  Block '<S124>/sample time' : Unused code path elimination
//  Block '<S125>/Delay Input2' : Unused code path elimination
//  Block '<S125>/Difference Inputs1' : Unused code path elimination
//  Block '<S125>/Difference Inputs2' : Unused code path elimination
//  Block '<S125>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S135>/Data Type Duplicate' : Unused code path elimination
//  Block '<S135>/Data Type Propagation' : Unused code path elimination
//  Block '<S135>/LowerRelop1' : Unused code path elimination
//  Block '<S135>/Switch' : Unused code path elimination
//  Block '<S135>/Switch2' : Unused code path elimination
//  Block '<S135>/UpperRelop' : Unused code path elimination
//  Block '<S125>/delta fall limit' : Unused code path elimination
//  Block '<S125>/delta rise limit' : Unused code path elimination
//  Block '<S125>/sample time' : Unused code path elimination
//  Block '<S126>/Data Type Duplicate' : Unused code path elimination
//  Block '<S126>/Data Type Propagation' : Unused code path elimination
//  Block '<S127>/Data Type Duplicate' : Unused code path elimination
//  Block '<S127>/Data Type Propagation' : Unused code path elimination
//  Block '<S127>/LowerRelop1' : Unused code path elimination
//  Block '<S127>/Switch' : Unused code path elimination
//  Block '<S127>/Switch2' : Unused code path elimination
//  Block '<S127>/UpperRelop' : Unused code path elimination
//  Block '<S128>/Data Type Duplicate' : Unused code path elimination
//  Block '<S128>/Data Type Propagation' : Unused code path elimination
//  Block '<S128>/LowerRelop1' : Unused code path elimination
//  Block '<S128>/Switch' : Unused code path elimination
//  Block '<S128>/Switch2' : Unused code path elimination
//  Block '<S128>/UpperRelop' : Unused code path elimination
//  Block '<S136>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S136>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S138>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S148>/Data Type Duplicate' : Unused code path elimination
//  Block '<S148>/Data Type Propagation' : Unused code path elimination
//  Block '<S139>/Delay Input2' : Unused code path elimination
//  Block '<S139>/Difference Inputs1' : Unused code path elimination
//  Block '<S139>/Difference Inputs2' : Unused code path elimination
//  Block '<S139>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S149>/Data Type Duplicate' : Unused code path elimination
//  Block '<S149>/Data Type Propagation' : Unused code path elimination
//  Block '<S149>/LowerRelop1' : Unused code path elimination
//  Block '<S149>/Switch' : Unused code path elimination
//  Block '<S149>/Switch2' : Unused code path elimination
//  Block '<S149>/UpperRelop' : Unused code path elimination
//  Block '<S139>/delta fall limit' : Unused code path elimination
//  Block '<S139>/delta rise limit' : Unused code path elimination
//  Block '<S139>/sample time' : Unused code path elimination
//  Block '<S140>/Delay Input2' : Unused code path elimination
//  Block '<S140>/Difference Inputs1' : Unused code path elimination
//  Block '<S140>/Difference Inputs2' : Unused code path elimination
//  Block '<S140>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S150>/Data Type Duplicate' : Unused code path elimination
//  Block '<S150>/Data Type Propagation' : Unused code path elimination
//  Block '<S150>/LowerRelop1' : Unused code path elimination
//  Block '<S150>/Switch' : Unused code path elimination
//  Block '<S150>/Switch2' : Unused code path elimination
//  Block '<S150>/UpperRelop' : Unused code path elimination
//  Block '<S140>/delta fall limit' : Unused code path elimination
//  Block '<S140>/delta rise limit' : Unused code path elimination
//  Block '<S140>/sample time' : Unused code path elimination
//  Block '<S141>/Data Type Duplicate' : Unused code path elimination
//  Block '<S141>/Data Type Propagation' : Unused code path elimination
//  Block '<S142>/Data Type Duplicate' : Unused code path elimination
//  Block '<S142>/Data Type Propagation' : Unused code path elimination
//  Block '<S142>/LowerRelop1' : Unused code path elimination
//  Block '<S142>/Switch' : Unused code path elimination
//  Block '<S142>/Switch2' : Unused code path elimination
//  Block '<S142>/UpperRelop' : Unused code path elimination
//  Block '<S143>/Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Propagation' : Unused code path elimination
//  Block '<S143>/LowerRelop1' : Unused code path elimination
//  Block '<S143>/Switch' : Unused code path elimination
//  Block '<S143>/Switch2' : Unused code path elimination
//  Block '<S143>/UpperRelop' : Unused code path elimination
//  Block '<S152>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S155>/Data Type Duplicate' : Unused code path elimination
//  Block '<S155>/Data Type Propagation' : Unused code path elimination
//  Block '<S153>/Data Type Duplicate' : Unused code path elimination
//  Block '<S153>/Data Type Propagation' : Unused code path elimination
//  Block '<S173>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S173>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S175>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S185>/Data Type Duplicate' : Unused code path elimination
//  Block '<S185>/Data Type Propagation' : Unused code path elimination
//  Block '<S176>/Delay Input2' : Unused code path elimination
//  Block '<S176>/Difference Inputs1' : Unused code path elimination
//  Block '<S176>/Difference Inputs2' : Unused code path elimination
//  Block '<S176>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Propagation' : Unused code path elimination
//  Block '<S186>/LowerRelop1' : Unused code path elimination
//  Block '<S186>/Switch' : Unused code path elimination
//  Block '<S186>/Switch2' : Unused code path elimination
//  Block '<S186>/UpperRelop' : Unused code path elimination
//  Block '<S176>/delta fall limit' : Unused code path elimination
//  Block '<S176>/delta rise limit' : Unused code path elimination
//  Block '<S176>/sample time' : Unused code path elimination
//  Block '<S177>/Delay Input2' : Unused code path elimination
//  Block '<S177>/Difference Inputs1' : Unused code path elimination
//  Block '<S177>/Difference Inputs2' : Unused code path elimination
//  Block '<S177>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S187>/Data Type Duplicate' : Unused code path elimination
//  Block '<S187>/Data Type Propagation' : Unused code path elimination
//  Block '<S187>/LowerRelop1' : Unused code path elimination
//  Block '<S187>/Switch' : Unused code path elimination
//  Block '<S187>/Switch2' : Unused code path elimination
//  Block '<S187>/UpperRelop' : Unused code path elimination
//  Block '<S177>/delta fall limit' : Unused code path elimination
//  Block '<S177>/delta rise limit' : Unused code path elimination
//  Block '<S177>/sample time' : Unused code path elimination
//  Block '<S178>/Data Type Duplicate' : Unused code path elimination
//  Block '<S178>/Data Type Propagation' : Unused code path elimination
//  Block '<S179>/Data Type Duplicate' : Unused code path elimination
//  Block '<S179>/Data Type Propagation' : Unused code path elimination
//  Block '<S179>/LowerRelop1' : Unused code path elimination
//  Block '<S179>/Switch' : Unused code path elimination
//  Block '<S179>/Switch2' : Unused code path elimination
//  Block '<S179>/UpperRelop' : Unused code path elimination
//  Block '<S180>/Data Type Duplicate' : Unused code path elimination
//  Block '<S180>/Data Type Propagation' : Unused code path elimination
//  Block '<S180>/LowerRelop1' : Unused code path elimination
//  Block '<S180>/Switch' : Unused code path elimination
//  Block '<S180>/Switch2' : Unused code path elimination
//  Block '<S180>/UpperRelop' : Unused code path elimination
//  Block '<S188>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S190>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S200>/Data Type Duplicate' : Unused code path elimination
//  Block '<S200>/Data Type Propagation' : Unused code path elimination
//  Block '<S191>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S201>/Data Type Duplicate' : Unused code path elimination
//  Block '<S201>/Data Type Propagation' : Unused code path elimination
//  Block '<S192>/Delay Input2' : Unused code path elimination
//  Block '<S192>/Difference Inputs1' : Unused code path elimination
//  Block '<S192>/Difference Inputs2' : Unused code path elimination
//  Block '<S192>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S202>/Data Type Duplicate' : Unused code path elimination
//  Block '<S202>/Data Type Propagation' : Unused code path elimination
//  Block '<S202>/LowerRelop1' : Unused code path elimination
//  Block '<S202>/Switch' : Unused code path elimination
//  Block '<S202>/Switch2' : Unused code path elimination
//  Block '<S202>/UpperRelop' : Unused code path elimination
//  Block '<S192>/delta fall limit' : Unused code path elimination
//  Block '<S192>/delta rise limit' : Unused code path elimination
//  Block '<S192>/sample time' : Unused code path elimination
//  Block '<S193>/Data Type Duplicate' : Unused code path elimination
//  Block '<S193>/Data Type Propagation' : Unused code path elimination
//  Block '<S194>/Data Type Duplicate' : Unused code path elimination
//  Block '<S194>/Data Type Propagation' : Unused code path elimination
//  Block '<S195>/Data Type Duplicate' : Unused code path elimination
//  Block '<S195>/Data Type Propagation' : Unused code path elimination
//  Block '<S195>/LowerRelop1' : Unused code path elimination
//  Block '<S195>/Switch' : Unused code path elimination
//  Block '<S195>/Switch2' : Unused code path elimination
//  Block '<S195>/UpperRelop' : Unused code path elimination
//  Block '<S204>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S207>/Data Type Duplicate' : Unused code path elimination
//  Block '<S207>/Data Type Propagation' : Unused code path elimination
//  Block '<S205>/Data Type Duplicate' : Unused code path elimination
//  Block '<S205>/Data Type Propagation' : Unused code path elimination
//  Block '<S212>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S212>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S214>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S224>/Data Type Duplicate' : Unused code path elimination
//  Block '<S224>/Data Type Propagation' : Unused code path elimination
//  Block '<S215>/Delay Input2' : Unused code path elimination
//  Block '<S215>/Difference Inputs1' : Unused code path elimination
//  Block '<S215>/Difference Inputs2' : Unused code path elimination
//  Block '<S215>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S225>/Data Type Duplicate' : Unused code path elimination
//  Block '<S225>/Data Type Propagation' : Unused code path elimination
//  Block '<S225>/LowerRelop1' : Unused code path elimination
//  Block '<S225>/Switch' : Unused code path elimination
//  Block '<S225>/Switch2' : Unused code path elimination
//  Block '<S225>/UpperRelop' : Unused code path elimination
//  Block '<S215>/delta fall limit' : Unused code path elimination
//  Block '<S215>/delta rise limit' : Unused code path elimination
//  Block '<S215>/sample time' : Unused code path elimination
//  Block '<S216>/Delay Input2' : Unused code path elimination
//  Block '<S216>/Difference Inputs1' : Unused code path elimination
//  Block '<S216>/Difference Inputs2' : Unused code path elimination
//  Block '<S216>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S226>/Data Type Duplicate' : Unused code path elimination
//  Block '<S226>/Data Type Propagation' : Unused code path elimination
//  Block '<S226>/LowerRelop1' : Unused code path elimination
//  Block '<S226>/Switch' : Unused code path elimination
//  Block '<S226>/Switch2' : Unused code path elimination
//  Block '<S226>/UpperRelop' : Unused code path elimination
//  Block '<S216>/delta fall limit' : Unused code path elimination
//  Block '<S216>/delta rise limit' : Unused code path elimination
//  Block '<S216>/sample time' : Unused code path elimination
//  Block '<S217>/Data Type Duplicate' : Unused code path elimination
//  Block '<S217>/Data Type Propagation' : Unused code path elimination
//  Block '<S218>/Data Type Duplicate' : Unused code path elimination
//  Block '<S218>/Data Type Propagation' : Unused code path elimination
//  Block '<S218>/LowerRelop1' : Unused code path elimination
//  Block '<S218>/Switch' : Unused code path elimination
//  Block '<S218>/Switch2' : Unused code path elimination
//  Block '<S218>/UpperRelop' : Unused code path elimination
//  Block '<S219>/Data Type Duplicate' : Unused code path elimination
//  Block '<S219>/Data Type Propagation' : Unused code path elimination
//  Block '<S219>/LowerRelop1' : Unused code path elimination
//  Block '<S219>/Switch' : Unused code path elimination
//  Block '<S219>/Switch2' : Unused code path elimination
//  Block '<S219>/UpperRelop' : Unused code path elimination
//  Block '<S227>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S227>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S229>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S239>/Data Type Duplicate' : Unused code path elimination
//  Block '<S239>/Data Type Propagation' : Unused code path elimination
//  Block '<S230>/Delay Input2' : Unused code path elimination
//  Block '<S230>/Difference Inputs1' : Unused code path elimination
//  Block '<S230>/Difference Inputs2' : Unused code path elimination
//  Block '<S230>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S240>/Data Type Duplicate' : Unused code path elimination
//  Block '<S240>/Data Type Propagation' : Unused code path elimination
//  Block '<S240>/LowerRelop1' : Unused code path elimination
//  Block '<S240>/Switch' : Unused code path elimination
//  Block '<S240>/Switch2' : Unused code path elimination
//  Block '<S240>/UpperRelop' : Unused code path elimination
//  Block '<S230>/delta fall limit' : Unused code path elimination
//  Block '<S230>/delta rise limit' : Unused code path elimination
//  Block '<S230>/sample time' : Unused code path elimination
//  Block '<S231>/Delay Input2' : Unused code path elimination
//  Block '<S231>/Difference Inputs1' : Unused code path elimination
//  Block '<S231>/Difference Inputs2' : Unused code path elimination
//  Block '<S231>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S241>/Data Type Duplicate' : Unused code path elimination
//  Block '<S241>/Data Type Propagation' : Unused code path elimination
//  Block '<S241>/LowerRelop1' : Unused code path elimination
//  Block '<S241>/Switch' : Unused code path elimination
//  Block '<S241>/Switch2' : Unused code path elimination
//  Block '<S241>/UpperRelop' : Unused code path elimination
//  Block '<S231>/delta fall limit' : Unused code path elimination
//  Block '<S231>/delta rise limit' : Unused code path elimination
//  Block '<S231>/sample time' : Unused code path elimination
//  Block '<S232>/Data Type Duplicate' : Unused code path elimination
//  Block '<S232>/Data Type Propagation' : Unused code path elimination
//  Block '<S233>/Data Type Duplicate' : Unused code path elimination
//  Block '<S233>/Data Type Propagation' : Unused code path elimination
//  Block '<S233>/LowerRelop1' : Unused code path elimination
//  Block '<S233>/Switch' : Unused code path elimination
//  Block '<S233>/Switch2' : Unused code path elimination
//  Block '<S233>/UpperRelop' : Unused code path elimination
//  Block '<S234>/Data Type Duplicate' : Unused code path elimination
//  Block '<S234>/Data Type Propagation' : Unused code path elimination
//  Block '<S234>/LowerRelop1' : Unused code path elimination
//  Block '<S234>/Switch' : Unused code path elimination
//  Block '<S234>/Switch2' : Unused code path elimination
//  Block '<S234>/UpperRelop' : Unused code path elimination
//  Block '<S242>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S242>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S244>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S254>/Data Type Duplicate' : Unused code path elimination
//  Block '<S254>/Data Type Propagation' : Unused code path elimination
//  Block '<S245>/Delay Input2' : Unused code path elimination
//  Block '<S245>/Difference Inputs1' : Unused code path elimination
//  Block '<S245>/Difference Inputs2' : Unused code path elimination
//  Block '<S245>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S255>/Data Type Duplicate' : Unused code path elimination
//  Block '<S255>/Data Type Propagation' : Unused code path elimination
//  Block '<S255>/LowerRelop1' : Unused code path elimination
//  Block '<S255>/Switch' : Unused code path elimination
//  Block '<S255>/Switch2' : Unused code path elimination
//  Block '<S255>/UpperRelop' : Unused code path elimination
//  Block '<S245>/delta fall limit' : Unused code path elimination
//  Block '<S245>/delta rise limit' : Unused code path elimination
//  Block '<S245>/sample time' : Unused code path elimination
//  Block '<S246>/Delay Input2' : Unused code path elimination
//  Block '<S246>/Difference Inputs1' : Unused code path elimination
//  Block '<S246>/Difference Inputs2' : Unused code path elimination
//  Block '<S246>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S256>/Data Type Duplicate' : Unused code path elimination
//  Block '<S256>/Data Type Propagation' : Unused code path elimination
//  Block '<S256>/LowerRelop1' : Unused code path elimination
//  Block '<S256>/Switch' : Unused code path elimination
//  Block '<S256>/Switch2' : Unused code path elimination
//  Block '<S256>/UpperRelop' : Unused code path elimination
//  Block '<S246>/delta fall limit' : Unused code path elimination
//  Block '<S246>/delta rise limit' : Unused code path elimination
//  Block '<S246>/sample time' : Unused code path elimination
//  Block '<S247>/Data Type Duplicate' : Unused code path elimination
//  Block '<S247>/Data Type Propagation' : Unused code path elimination
//  Block '<S248>/Data Type Duplicate' : Unused code path elimination
//  Block '<S248>/Data Type Propagation' : Unused code path elimination
//  Block '<S248>/LowerRelop1' : Unused code path elimination
//  Block '<S248>/Switch' : Unused code path elimination
//  Block '<S248>/Switch2' : Unused code path elimination
//  Block '<S248>/UpperRelop' : Unused code path elimination
//  Block '<S249>/Data Type Duplicate' : Unused code path elimination
//  Block '<S249>/Data Type Propagation' : Unused code path elimination
//  Block '<S249>/LowerRelop1' : Unused code path elimination
//  Block '<S249>/Switch' : Unused code path elimination
//  Block '<S249>/Switch2' : Unused code path elimination
//  Block '<S249>/UpperRelop' : Unused code path elimination
//  Block '<S258>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S261>/Data Type Duplicate' : Unused code path elimination
//  Block '<S261>/Data Type Propagation' : Unused code path elimination
//  Block '<S259>/Data Type Duplicate' : Unused code path elimination
//  Block '<S259>/Data Type Propagation' : Unused code path elimination


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'fcsModel'
//  '<S1>'   : 'fcsModel/Allocation'
//  '<S2>'   : 'fcsModel/Inner Loop Controller'
//  '<S3>'   : 'fcsModel/Outer Loop Controller'
//  '<S4>'   : 'fcsModel/RC Interpreter'
//  '<S5>'   : 'fcsModel/Allocation/Angular Velocity Conversion'
//  '<S6>'   : 'fcsModel/Allocation/Compare To Constant'
//  '<S7>'   : 'fcsModel/Allocation/For Each Subsystem'
//  '<S8>'   : 'fcsModel/Inner Loop Controller/Angular Rate Controller'
//  '<S9>'   : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs'
//  '<S10>'  : 'fcsModel/Inner Loop Controller/Attitude Controller'
//  '<S11>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem'
//  '<S12>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S13>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S14>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug'
//  '<S15>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S16>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S17>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S18>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S19>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S20>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S21>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S22>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S23>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S24>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S25>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S26>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S27>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S28>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S29>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S30>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S31>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S32>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S33>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S34>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S35>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S36>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S37>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S38>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S39>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S40>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S41>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S42>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S43>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S44>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S45>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S46>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S47>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S48>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S49>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S50>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/Compare To Constant'
//  '<S51>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/EulerRates2BodyRates'
//  '<S52>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control'
//  '<S53>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant'
//  '<S54>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant1'
//  '<S55>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block'
//  '<S56>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1'
//  '<S57>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pickAttitudeCmdAndMeas'
//  '<S58>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug'
//  '<S59>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S60>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S61>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S62>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S63>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S64>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic'
//  '<S65>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic1'
//  '<S66>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic2'
//  '<S67>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S68>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S69>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S70>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S71>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S72>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S73>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S74>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S75>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S76>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S77>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S78>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S79>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic'
//  '<S80>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S81>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S82>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S83>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S84>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S85>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S86>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S87>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S88>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S89>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S90>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Rate Limiter Dynamic'
//  '<S91>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Saturation Dynamic'
//  '<S92>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S93>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S94>'  : 'fcsModel/Outer Loop Controller/Compare To Constant'
//  '<S95>'  : 'fcsModel/Outer Loop Controller/Compare To Constant1'
//  '<S96>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl'
//  '<S97>'  : 'fcsModel/Outer Loop Controller/assembleOuterLoopToInnerLoopBus'
//  '<S98>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs'
//  '<S99>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller'
//  '<S100>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller'
//  '<S101>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs/Compare To Constant1'
//  '<S102>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs/Compare To Constant2'
//  '<S103>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs'
//  '<S104>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control'
//  '<S105>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant'
//  '<S106>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant1'
//  '<S107>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant2'
//  '<S108>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant3'
//  '<S109>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant4'
//  '<S110>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant5'
//  '<S111>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant6'
//  '<S112>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter'
//  '<S113>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter1'
//  '<S114>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter2'
//  '<S115>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter/holdOutputAtCenter'
//  '<S116>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter1/holdOutputAtCenter'
//  '<S117>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter2/holdOutputAtCenter'
//  '<S118>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block'
//  '<S119>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1'
//  '<S120>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug'
//  '<S121>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S122>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S123>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S124>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S125>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S126>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic'
//  '<S127>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic1'
//  '<S128>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic2'
//  '<S129>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S130>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S131>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S132>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S133>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S134>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S135>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S136>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S137>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S138>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S139>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S140>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S141>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic'
//  '<S142>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S143>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S144>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S145>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S146>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S147>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S148>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S149>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S150>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S151>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S152>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic'
//  '<S153>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Saturation Dynamic'
//  '<S154>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S155>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S156>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs'
//  '<S157>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem'
//  '<S158>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant'
//  '<S159>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant1'
//  '<S160>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant2'
//  '<S161>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/accelZKiSelectorVariantSubsystem'
//  '<S162>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/holdOutputAtCenter'
//  '<S163>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/hoverThrustVariantSubsystem'
//  '<S164>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem'
//  '<S165>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/accelZKiSelectorVariantSubsystem/accelZCtrlKiPassThrough'
//  '<S166>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/holdOutputAtCenter/holdOutputAtCenter'
//  '<S167>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/hoverThrustVariantSubsystem/constantHoverThrust'
//  '<S168>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl'
//  '<S169>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller'
//  '<S170>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block'
//  '<S171>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1'
//  '<S172>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug'
//  '<S173>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S174>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S175>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S176>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S177>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S178>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic'
//  '<S179>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic1'
//  '<S180>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic2'
//  '<S181>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S182>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S183>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S184>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S185>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S186>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S187>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S188>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S189>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S190>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S191>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S192>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S193>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic'
//  '<S194>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S195>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S196>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S197>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S198>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S199>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S200>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S201>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S202>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S203>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S204>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Rate Limiter Dynamic'
//  '<S205>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Saturation Dynamic'
//  '<S206>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S207>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S208>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S209>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S210>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2'
//  '<S211>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug'
//  '<S212>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S213>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S214>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S215>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S216>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S217>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S218>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S219>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S220>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S221>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S222>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S223>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S224>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S225>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S226>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S227>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S228>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S229>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S230>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S231>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S232>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S233>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S234>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S235>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S236>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S237>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S238>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S239>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S240>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S241>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S242>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter'
//  '<S243>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter'
//  '<S244>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic'
//  '<S245>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic1'
//  '<S246>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic2'
//  '<S247>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic'
//  '<S248>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic1'
//  '<S249>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic2'
//  '<S250>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S251>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S252>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S253>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S254>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S255>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S256>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S257>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S258>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S259>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S260>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S261>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S262>' : 'fcsModel/RC Interpreter/Chart'
//  '<S263>' : 'fcsModel/RC Interpreter/Interpret RC In Cmds'


//-
//  Requirements for '<Root>': fcsModel

#endif                                 // RTW_HEADER_fcsModel_h_

//
// File trailer for generated code.
//
// [EOF]
//
