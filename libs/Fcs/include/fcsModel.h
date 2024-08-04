//
// File: fcsModel.h
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
#ifndef RTW_HEADER_fcsModel_h_
#define RTW_HEADER_fcsModel_h_
#include "rtwtypes.h"
#include "fcsModel_types.h"
#include <array>

// External data declarations for dependent source files
extern const busXyBodyAccelCtrIDebug fcsModel_rtZbusXyBodyAccelCtrIDebug;// busXyBodyAccelCtrIDebug ground 
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

  // Block signals and states (default storage) for system '<S24>/Discrete First Order Deriv Filter' 
  struct DW_DiscreteFirstOrderDerivFil_T {
    std::array<real_T, 2> num;
                      // '<S55>/Compute Deriv Filter Numerator And Denominator'
    std::array<real_T, 2> den;
                      // '<S55>/Compute Deriv Filter Numerator And Denominator'
    real_T DiscreteTransferFcn_states; // '<S55>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S21>/pidWithDebug'
  struct DW_pidWithDebug_fcsModel_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                   // '<S24>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S24>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S56>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S24>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S24>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S24>/Discrete-Time Integrator'
    uint8_T DiscreteTimeIntegrator_IC_LOADI;// '<S24>/Discrete-Time Integrator'
  };

  // Block signals and states (default storage) for system '<S21>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_f_T {
    std::array<real_T, 3> num;
                            // '<S41>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                            // '<S41>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S41>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S42>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S41>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S18>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_c_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S21>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S21>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S21>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S21>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S20>/Attitude Control' 
  struct DW_CoreSubsys_fcsModel_i_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S62>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S62>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S62>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S62>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S113>/holdOutputAtCenter1' 
  struct DW_holdOutputAtCenter1_fcsMod_T {
    real_T last_input;                 // '<S123>/holdOutputAtCenter'
  };

  // Block signals and states (default storage) for system '<S114>/pidWithDebug' 
  struct DW_pidWithDebug_fcsModel_i_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                  // '<S130>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S130>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S162>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S130>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S130>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S130>/Discrete-Time Integrator'
    uint8_T DiscreteTimeIntegrator_IC_LOADI;// '<S130>/Discrete-Time Integrator' 
  };

  // Block signals and states (default storage) for system '<S114>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_g_T {
    std::array<real_T, 3> num;
                           // '<S147>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                           // '<S147>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S147>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S148>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S147>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S109>/NED Position Control' 
  struct DW_CoreSubsys_fcsModel_b_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S114>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S114>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S114>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S114>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S110>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_p_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S167>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S167>/Signal Conditioning Block1' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock2;// '<S167>/Signal Conditioning Block2' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S167>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S167>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW_fcsModel_T {
    std::array<DW_CoreSubsys_fcsModel_p_T, 3> CoreSubsys_i;// '<S110>/For Each Subsystem' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S182>/pidWithDebug'
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S182>/Signal Conditioning Block' 
    std::array<DW_CoreSubsys_fcsModel_b_T, 3> CoreSubsys_g;// '<S109>/NED Position Control' 
    DW_holdOutputAtCenter1_fcsMod_T holdOutputAtCenter2;// '<S113>/holdOutputAtCenter2' 
    DW_holdOutputAtCenter1_fcsMod_T holdOutputAtCenter1;// '<S113>/holdOutputAtCenter1' 
    std::array<DW_CoreSubsys_fcsModel_i_T, 3> CoreSubsys_p;// '<S20>/Attitude Control' 
    std::array<DW_CoreSubsys_fcsModel_c_T, 3> CoreSubsys_a;// '<S18>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_T, 4> CoreSubsys;// '<S1>/For Each Subsystem' 
    busOuterLoopCtrlDebug RateTransition_Buffer0;// '<Root>/Rate Transition'
    busOuterLoopToInnerLoop Switch2;   // '<S3>/Switch2'
    busRcOutCmds rcOutCmds;            // '<S4>/Interpret RC In Cmds'
    std::array<real_T, 4> DiscreteTransferFcn_states_d;// '<S1>/Discrete Transfer Fcn' 
    std::array<real_T, 2> DiscreteTransferFcn_states_n;// '<S201>/Discrete Transfer Fcn' 
    std::array<real_T, 2> DiscreteTransferFcn_states_nh;// '<S202>/Discrete Transfer Fcn' 
    real_T DiscreteTransferFcn_states; // '<S179>/Discrete Transfer Fcn'
    real_T UnitDelay_DSTATE;           // '<S182>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S15>/Unit Delay1'
    real_T UnitDelay_DSTATE_j;         // '<S15>/Unit Delay'
    real_T DiscreteTransferFcn_states_c;// '<S15>/Discrete Transfer Fcn'
    real_T DelayInput2_DSTATE;         // '<S204>/Delay Input2'
    real_T DelayInput2_DSTATE_e;       // '<S203>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S179>/Discrete Transfer Fcn'
    real_T NextOutput;                 // '<S16>/White Noise'
    real_T last_input;                 // '<S172>/holdOutputAtCenter'
    real_T last_input_c;               // '<S122>/holdOutputAtCenter'
    int32_T durationCounter_1;         // '<S4>/Chart'
    int32_T durationCounter_1_j;       // '<S4>/Chart'
    uint32_T RandSeed;                 // '<S16>/White Noise'
    enumChirpTrigger UnitDelay_DSTATE_g;// '<S4>/Unit Delay'
    uint16_T temporalCounter_i1;       // '<S4>/Chart'
    uint8_T chirpCount_;               // '<S4>/Interpret RC In Cmds'
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
    //    '<S166>/Constant'

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
    //  Referenced by: '<S20>/Constant'

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

  // private member function(s) for subsystem '<S24>/Discrete First Order Deriv Filter'
  static void f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
    rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
    rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW);

  // private member function(s) for subsystem '<S21>/pidWithDebug'
  static void fcsModel_pidWithDebug_Init(DW_pidWithDebug_fcsModel_T *localDW);
  static void fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc,
    const busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T
    *rty_ctrlCmd, busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_T *localDW);

  // private member function(s) for subsystem '<S40>/Compute Natural Frequency'
  static void fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps, real_T
    rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps);

  // private member function(s) for subsystem '<S40>/Compute Numerator And Denominator'
  static void ComputeNumeratorAndDenominator(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_rateNum[3], real_T rty_accelNum[3],
    real_T rty_den[3], real_T rtp_sampleTime_s);

  // private member function(s) for subsystem '<S41>/Compute Filter Numerator And Denominator'
  static void ComputeFilterNumeratorAndD_Init(real_T rty_num[3], real_T rty_den
    [3]);
  static void ComputeFilterNumeratorAndDenomi(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_num[3], real_T rty_den[3], real_T
    rtp_sampleTime_s);

  // private member function(s) for subsystem '<S21>/Signal Conditioning Block1'
  static void f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
    *localDW);
  static void fcsMod_SignalConditioningBlock1(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW);

  // private member function(s) for subsystem '<S113>/holdOutputAtCenter1'
  static void fcsModel_holdOutputAtCenter1(real_T rtu_input, real_T rtu_trigger,
    real_T *rty_output, boolean_T *rty_atCenter, DW_holdOutputAtCenter1_fcsMod_T
    *localDW);

  // private member function(s) for subsystem '<S114>/pidWithDebug'
  static void fcsModel_pidWithDebug_m_Init(DW_pidWithDebug_fcsModel_i_T *localDW);
  static void fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, real_T rtu_integratorIc,
    const busPidParams *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T
    *rty_ctrlCmd, busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_i_T *localDW);

  // private member function(s) for subsystem '<S114>/Signal Conditioning Block1'
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
//  Block '<S25>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S25>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S27>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S37>/Data Type Duplicate' : Unused code path elimination
//  Block '<S37>/Data Type Propagation' : Unused code path elimination
//  Block '<S28>/Delay Input2' : Unused code path elimination
//  Block '<S28>/Difference Inputs1' : Unused code path elimination
//  Block '<S28>/Difference Inputs2' : Unused code path elimination
//  Block '<S28>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S38>/Data Type Duplicate' : Unused code path elimination
//  Block '<S38>/Data Type Propagation' : Unused code path elimination
//  Block '<S38>/LowerRelop1' : Unused code path elimination
//  Block '<S38>/Switch' : Unused code path elimination
//  Block '<S38>/Switch2' : Unused code path elimination
//  Block '<S38>/UpperRelop' : Unused code path elimination
//  Block '<S28>/delta fall limit' : Unused code path elimination
//  Block '<S28>/delta rise limit' : Unused code path elimination
//  Block '<S28>/sample time' : Unused code path elimination
//  Block '<S29>/Delay Input2' : Unused code path elimination
//  Block '<S29>/Difference Inputs1' : Unused code path elimination
//  Block '<S29>/Difference Inputs2' : Unused code path elimination
//  Block '<S29>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S39>/Data Type Duplicate' : Unused code path elimination
//  Block '<S39>/Data Type Propagation' : Unused code path elimination
//  Block '<S39>/LowerRelop1' : Unused code path elimination
//  Block '<S39>/Switch' : Unused code path elimination
//  Block '<S39>/Switch2' : Unused code path elimination
//  Block '<S39>/UpperRelop' : Unused code path elimination
//  Block '<S29>/delta fall limit' : Unused code path elimination
//  Block '<S29>/delta rise limit' : Unused code path elimination
//  Block '<S29>/sample time' : Unused code path elimination
//  Block '<S30>/Data Type Duplicate' : Unused code path elimination
//  Block '<S30>/Data Type Propagation' : Unused code path elimination
//  Block '<S31>/Data Type Duplicate' : Unused code path elimination
//  Block '<S31>/Data Type Propagation' : Unused code path elimination
//  Block '<S31>/LowerRelop1' : Unused code path elimination
//  Block '<S31>/Switch' : Unused code path elimination
//  Block '<S31>/Switch2' : Unused code path elimination
//  Block '<S31>/UpperRelop' : Unused code path elimination
//  Block '<S32>/Data Type Duplicate' : Unused code path elimination
//  Block '<S32>/Data Type Propagation' : Unused code path elimination
//  Block '<S32>/LowerRelop1' : Unused code path elimination
//  Block '<S32>/Switch' : Unused code path elimination
//  Block '<S32>/Switch2' : Unused code path elimination
//  Block '<S32>/UpperRelop' : Unused code path elimination
//  Block '<S40>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S40>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S42>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S52>/Data Type Duplicate' : Unused code path elimination
//  Block '<S52>/Data Type Propagation' : Unused code path elimination
//  Block '<S43>/Delay Input2' : Unused code path elimination
//  Block '<S43>/Difference Inputs1' : Unused code path elimination
//  Block '<S43>/Difference Inputs2' : Unused code path elimination
//  Block '<S43>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Propagation' : Unused code path elimination
//  Block '<S53>/LowerRelop1' : Unused code path elimination
//  Block '<S53>/Switch' : Unused code path elimination
//  Block '<S53>/Switch2' : Unused code path elimination
//  Block '<S53>/UpperRelop' : Unused code path elimination
//  Block '<S43>/delta fall limit' : Unused code path elimination
//  Block '<S43>/delta rise limit' : Unused code path elimination
//  Block '<S43>/sample time' : Unused code path elimination
//  Block '<S44>/Delay Input2' : Unused code path elimination
//  Block '<S44>/Difference Inputs1' : Unused code path elimination
//  Block '<S44>/Difference Inputs2' : Unused code path elimination
//  Block '<S44>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S54>/Data Type Duplicate' : Unused code path elimination
//  Block '<S54>/Data Type Propagation' : Unused code path elimination
//  Block '<S54>/LowerRelop1' : Unused code path elimination
//  Block '<S54>/Switch' : Unused code path elimination
//  Block '<S54>/Switch2' : Unused code path elimination
//  Block '<S54>/UpperRelop' : Unused code path elimination
//  Block '<S44>/delta fall limit' : Unused code path elimination
//  Block '<S44>/delta rise limit' : Unused code path elimination
//  Block '<S44>/sample time' : Unused code path elimination
//  Block '<S45>/Data Type Duplicate' : Unused code path elimination
//  Block '<S45>/Data Type Propagation' : Unused code path elimination
//  Block '<S46>/Data Type Duplicate' : Unused code path elimination
//  Block '<S46>/Data Type Propagation' : Unused code path elimination
//  Block '<S46>/LowerRelop1' : Unused code path elimination
//  Block '<S46>/Switch' : Unused code path elimination
//  Block '<S46>/Switch2' : Unused code path elimination
//  Block '<S46>/UpperRelop' : Unused code path elimination
//  Block '<S47>/Data Type Duplicate' : Unused code path elimination
//  Block '<S47>/Data Type Propagation' : Unused code path elimination
//  Block '<S47>/LowerRelop1' : Unused code path elimination
//  Block '<S47>/Switch' : Unused code path elimination
//  Block '<S47>/Switch2' : Unused code path elimination
//  Block '<S47>/UpperRelop' : Unused code path elimination
//  Block '<S56>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S59>/Data Type Duplicate' : Unused code path elimination
//  Block '<S59>/Data Type Propagation' : Unused code path elimination
//  Block '<S57>/Data Type Duplicate' : Unused code path elimination
//  Block '<S57>/Data Type Propagation' : Unused code path elimination
//  Block '<S69>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S69>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S71>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S81>/Data Type Duplicate' : Unused code path elimination
//  Block '<S81>/Data Type Propagation' : Unused code path elimination
//  Block '<S72>/Delay Input2' : Unused code path elimination
//  Block '<S72>/Difference Inputs1' : Unused code path elimination
//  Block '<S72>/Difference Inputs2' : Unused code path elimination
//  Block '<S72>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S82>/Data Type Duplicate' : Unused code path elimination
//  Block '<S82>/Data Type Propagation' : Unused code path elimination
//  Block '<S82>/LowerRelop1' : Unused code path elimination
//  Block '<S82>/Switch' : Unused code path elimination
//  Block '<S82>/Switch2' : Unused code path elimination
//  Block '<S82>/UpperRelop' : Unused code path elimination
//  Block '<S72>/delta fall limit' : Unused code path elimination
//  Block '<S72>/delta rise limit' : Unused code path elimination
//  Block '<S72>/sample time' : Unused code path elimination
//  Block '<S73>/Delay Input2' : Unused code path elimination
//  Block '<S73>/Difference Inputs1' : Unused code path elimination
//  Block '<S73>/Difference Inputs2' : Unused code path elimination
//  Block '<S73>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S83>/Data Type Duplicate' : Unused code path elimination
//  Block '<S83>/Data Type Propagation' : Unused code path elimination
//  Block '<S83>/LowerRelop1' : Unused code path elimination
//  Block '<S83>/Switch' : Unused code path elimination
//  Block '<S83>/Switch2' : Unused code path elimination
//  Block '<S83>/UpperRelop' : Unused code path elimination
//  Block '<S73>/delta fall limit' : Unused code path elimination
//  Block '<S73>/delta rise limit' : Unused code path elimination
//  Block '<S73>/sample time' : Unused code path elimination
//  Block '<S74>/Data Type Duplicate' : Unused code path elimination
//  Block '<S74>/Data Type Propagation' : Unused code path elimination
//  Block '<S75>/Data Type Duplicate' : Unused code path elimination
//  Block '<S75>/Data Type Propagation' : Unused code path elimination
//  Block '<S75>/LowerRelop1' : Unused code path elimination
//  Block '<S75>/Switch' : Unused code path elimination
//  Block '<S75>/Switch2' : Unused code path elimination
//  Block '<S75>/UpperRelop' : Unused code path elimination
//  Block '<S76>/Data Type Duplicate' : Unused code path elimination
//  Block '<S76>/Data Type Propagation' : Unused code path elimination
//  Block '<S76>/LowerRelop1' : Unused code path elimination
//  Block '<S76>/Switch' : Unused code path elimination
//  Block '<S76>/Switch2' : Unused code path elimination
//  Block '<S76>/UpperRelop' : Unused code path elimination
//  Block '<S84>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S84>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S86>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S96>/Data Type Duplicate' : Unused code path elimination
//  Block '<S96>/Data Type Propagation' : Unused code path elimination
//  Block '<S87>/Delay Input2' : Unused code path elimination
//  Block '<S87>/Difference Inputs1' : Unused code path elimination
//  Block '<S87>/Difference Inputs2' : Unused code path elimination
//  Block '<S87>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S97>/Data Type Duplicate' : Unused code path elimination
//  Block '<S97>/Data Type Propagation' : Unused code path elimination
//  Block '<S97>/LowerRelop1' : Unused code path elimination
//  Block '<S97>/Switch' : Unused code path elimination
//  Block '<S97>/Switch2' : Unused code path elimination
//  Block '<S97>/UpperRelop' : Unused code path elimination
//  Block '<S87>/delta fall limit' : Unused code path elimination
//  Block '<S87>/delta rise limit' : Unused code path elimination
//  Block '<S87>/sample time' : Unused code path elimination
//  Block '<S88>/Delay Input2' : Unused code path elimination
//  Block '<S88>/Difference Inputs1' : Unused code path elimination
//  Block '<S88>/Difference Inputs2' : Unused code path elimination
//  Block '<S88>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S98>/Data Type Duplicate' : Unused code path elimination
//  Block '<S98>/Data Type Propagation' : Unused code path elimination
//  Block '<S98>/LowerRelop1' : Unused code path elimination
//  Block '<S98>/Switch' : Unused code path elimination
//  Block '<S98>/Switch2' : Unused code path elimination
//  Block '<S98>/UpperRelop' : Unused code path elimination
//  Block '<S88>/delta fall limit' : Unused code path elimination
//  Block '<S88>/delta rise limit' : Unused code path elimination
//  Block '<S88>/sample time' : Unused code path elimination
//  Block '<S89>/Data Type Duplicate' : Unused code path elimination
//  Block '<S89>/Data Type Propagation' : Unused code path elimination
//  Block '<S90>/Data Type Duplicate' : Unused code path elimination
//  Block '<S90>/Data Type Propagation' : Unused code path elimination
//  Block '<S90>/LowerRelop1' : Unused code path elimination
//  Block '<S90>/Switch' : Unused code path elimination
//  Block '<S90>/Switch2' : Unused code path elimination
//  Block '<S90>/UpperRelop' : Unused code path elimination
//  Block '<S91>/Data Type Duplicate' : Unused code path elimination
//  Block '<S91>/Data Type Propagation' : Unused code path elimination
//  Block '<S91>/LowerRelop1' : Unused code path elimination
//  Block '<S91>/Switch' : Unused code path elimination
//  Block '<S91>/Switch2' : Unused code path elimination
//  Block '<S91>/UpperRelop' : Unused code path elimination
//  Block '<S100>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S103>/Data Type Duplicate' : Unused code path elimination
//  Block '<S103>/Data Type Propagation' : Unused code path elimination
//  Block '<S101>/Data Type Duplicate' : Unused code path elimination
//  Block '<S101>/Data Type Propagation' : Unused code path elimination
//  Block '<S131>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S131>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S133>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Propagation' : Unused code path elimination
//  Block '<S134>/Delay Input2' : Unused code path elimination
//  Block '<S134>/Difference Inputs1' : Unused code path elimination
//  Block '<S134>/Difference Inputs2' : Unused code path elimination
//  Block '<S134>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S144>/Data Type Duplicate' : Unused code path elimination
//  Block '<S144>/Data Type Propagation' : Unused code path elimination
//  Block '<S144>/LowerRelop1' : Unused code path elimination
//  Block '<S144>/Switch' : Unused code path elimination
//  Block '<S144>/Switch2' : Unused code path elimination
//  Block '<S144>/UpperRelop' : Unused code path elimination
//  Block '<S134>/delta fall limit' : Unused code path elimination
//  Block '<S134>/delta rise limit' : Unused code path elimination
//  Block '<S134>/sample time' : Unused code path elimination
//  Block '<S135>/Delay Input2' : Unused code path elimination
//  Block '<S135>/Difference Inputs1' : Unused code path elimination
//  Block '<S135>/Difference Inputs2' : Unused code path elimination
//  Block '<S135>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S145>/Data Type Duplicate' : Unused code path elimination
//  Block '<S145>/Data Type Propagation' : Unused code path elimination
//  Block '<S145>/LowerRelop1' : Unused code path elimination
//  Block '<S145>/Switch' : Unused code path elimination
//  Block '<S145>/Switch2' : Unused code path elimination
//  Block '<S145>/UpperRelop' : Unused code path elimination
//  Block '<S135>/delta fall limit' : Unused code path elimination
//  Block '<S135>/delta rise limit' : Unused code path elimination
//  Block '<S135>/sample time' : Unused code path elimination
//  Block '<S136>/Data Type Duplicate' : Unused code path elimination
//  Block '<S136>/Data Type Propagation' : Unused code path elimination
//  Block '<S137>/Data Type Duplicate' : Unused code path elimination
//  Block '<S137>/Data Type Propagation' : Unused code path elimination
//  Block '<S137>/LowerRelop1' : Unused code path elimination
//  Block '<S137>/Switch' : Unused code path elimination
//  Block '<S137>/Switch2' : Unused code path elimination
//  Block '<S137>/UpperRelop' : Unused code path elimination
//  Block '<S138>/Data Type Duplicate' : Unused code path elimination
//  Block '<S138>/Data Type Propagation' : Unused code path elimination
//  Block '<S138>/LowerRelop1' : Unused code path elimination
//  Block '<S138>/Switch' : Unused code path elimination
//  Block '<S138>/Switch2' : Unused code path elimination
//  Block '<S138>/UpperRelop' : Unused code path elimination
//  Block '<S146>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S146>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S148>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S158>/Data Type Duplicate' : Unused code path elimination
//  Block '<S158>/Data Type Propagation' : Unused code path elimination
//  Block '<S149>/Delay Input2' : Unused code path elimination
//  Block '<S149>/Difference Inputs1' : Unused code path elimination
//  Block '<S149>/Difference Inputs2' : Unused code path elimination
//  Block '<S149>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Propagation' : Unused code path elimination
//  Block '<S159>/LowerRelop1' : Unused code path elimination
//  Block '<S159>/Switch' : Unused code path elimination
//  Block '<S159>/Switch2' : Unused code path elimination
//  Block '<S159>/UpperRelop' : Unused code path elimination
//  Block '<S149>/delta fall limit' : Unused code path elimination
//  Block '<S149>/delta rise limit' : Unused code path elimination
//  Block '<S149>/sample time' : Unused code path elimination
//  Block '<S150>/Delay Input2' : Unused code path elimination
//  Block '<S150>/Difference Inputs1' : Unused code path elimination
//  Block '<S150>/Difference Inputs2' : Unused code path elimination
//  Block '<S150>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Propagation' : Unused code path elimination
//  Block '<S160>/LowerRelop1' : Unused code path elimination
//  Block '<S160>/Switch' : Unused code path elimination
//  Block '<S160>/Switch2' : Unused code path elimination
//  Block '<S160>/UpperRelop' : Unused code path elimination
//  Block '<S150>/delta fall limit' : Unused code path elimination
//  Block '<S150>/delta rise limit' : Unused code path elimination
//  Block '<S150>/sample time' : Unused code path elimination
//  Block '<S151>/Data Type Duplicate' : Unused code path elimination
//  Block '<S151>/Data Type Propagation' : Unused code path elimination
//  Block '<S152>/Data Type Duplicate' : Unused code path elimination
//  Block '<S152>/Data Type Propagation' : Unused code path elimination
//  Block '<S152>/LowerRelop1' : Unused code path elimination
//  Block '<S152>/Switch' : Unused code path elimination
//  Block '<S152>/Switch2' : Unused code path elimination
//  Block '<S152>/UpperRelop' : Unused code path elimination
//  Block '<S153>/Data Type Duplicate' : Unused code path elimination
//  Block '<S153>/Data Type Propagation' : Unused code path elimination
//  Block '<S153>/LowerRelop1' : Unused code path elimination
//  Block '<S153>/Switch' : Unused code path elimination
//  Block '<S153>/Switch2' : Unused code path elimination
//  Block '<S153>/UpperRelop' : Unused code path elimination
//  Block '<S162>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S165>/Data Type Duplicate' : Unused code path elimination
//  Block '<S165>/Data Type Propagation' : Unused code path elimination
//  Block '<S163>/Data Type Duplicate' : Unused code path elimination
//  Block '<S163>/Data Type Propagation' : Unused code path elimination
//  Block '<S186>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S186>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S188>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S198>/Data Type Duplicate' : Unused code path elimination
//  Block '<S198>/Data Type Propagation' : Unused code path elimination
//  Block '<S189>/Delay Input2' : Unused code path elimination
//  Block '<S189>/Difference Inputs1' : Unused code path elimination
//  Block '<S189>/Difference Inputs2' : Unused code path elimination
//  Block '<S189>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S199>/Data Type Duplicate' : Unused code path elimination
//  Block '<S199>/Data Type Propagation' : Unused code path elimination
//  Block '<S199>/LowerRelop1' : Unused code path elimination
//  Block '<S199>/Switch' : Unused code path elimination
//  Block '<S199>/Switch2' : Unused code path elimination
//  Block '<S199>/UpperRelop' : Unused code path elimination
//  Block '<S189>/delta fall limit' : Unused code path elimination
//  Block '<S189>/delta rise limit' : Unused code path elimination
//  Block '<S189>/sample time' : Unused code path elimination
//  Block '<S190>/Delay Input2' : Unused code path elimination
//  Block '<S190>/Difference Inputs1' : Unused code path elimination
//  Block '<S190>/Difference Inputs2' : Unused code path elimination
//  Block '<S190>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S200>/Data Type Duplicate' : Unused code path elimination
//  Block '<S200>/Data Type Propagation' : Unused code path elimination
//  Block '<S200>/LowerRelop1' : Unused code path elimination
//  Block '<S200>/Switch' : Unused code path elimination
//  Block '<S200>/Switch2' : Unused code path elimination
//  Block '<S200>/UpperRelop' : Unused code path elimination
//  Block '<S190>/delta fall limit' : Unused code path elimination
//  Block '<S190>/delta rise limit' : Unused code path elimination
//  Block '<S190>/sample time' : Unused code path elimination
//  Block '<S191>/Data Type Duplicate' : Unused code path elimination
//  Block '<S191>/Data Type Propagation' : Unused code path elimination
//  Block '<S192>/Data Type Duplicate' : Unused code path elimination
//  Block '<S192>/Data Type Propagation' : Unused code path elimination
//  Block '<S192>/LowerRelop1' : Unused code path elimination
//  Block '<S192>/Switch' : Unused code path elimination
//  Block '<S192>/Switch2' : Unused code path elimination
//  Block '<S192>/UpperRelop' : Unused code path elimination
//  Block '<S193>/Data Type Duplicate' : Unused code path elimination
//  Block '<S193>/Data Type Propagation' : Unused code path elimination
//  Block '<S193>/LowerRelop1' : Unused code path elimination
//  Block '<S193>/Switch' : Unused code path elimination
//  Block '<S193>/Switch2' : Unused code path elimination
//  Block '<S193>/UpperRelop' : Unused code path elimination
//  Block '<S201>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S203>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S213>/Data Type Duplicate' : Unused code path elimination
//  Block '<S213>/Data Type Propagation' : Unused code path elimination
//  Block '<S204>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S214>/Data Type Duplicate' : Unused code path elimination
//  Block '<S214>/Data Type Propagation' : Unused code path elimination
//  Block '<S205>/Delay Input2' : Unused code path elimination
//  Block '<S205>/Difference Inputs1' : Unused code path elimination
//  Block '<S205>/Difference Inputs2' : Unused code path elimination
//  Block '<S205>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S215>/Data Type Duplicate' : Unused code path elimination
//  Block '<S215>/Data Type Propagation' : Unused code path elimination
//  Block '<S215>/LowerRelop1' : Unused code path elimination
//  Block '<S215>/Switch' : Unused code path elimination
//  Block '<S215>/Switch2' : Unused code path elimination
//  Block '<S215>/UpperRelop' : Unused code path elimination
//  Block '<S205>/delta fall limit' : Unused code path elimination
//  Block '<S205>/delta rise limit' : Unused code path elimination
//  Block '<S205>/sample time' : Unused code path elimination
//  Block '<S206>/Data Type Duplicate' : Unused code path elimination
//  Block '<S206>/Data Type Propagation' : Unused code path elimination
//  Block '<S207>/Data Type Duplicate' : Unused code path elimination
//  Block '<S207>/Data Type Propagation' : Unused code path elimination
//  Block '<S208>/Data Type Duplicate' : Unused code path elimination
//  Block '<S208>/Data Type Propagation' : Unused code path elimination
//  Block '<S208>/LowerRelop1' : Unused code path elimination
//  Block '<S208>/Switch' : Unused code path elimination
//  Block '<S208>/Switch2' : Unused code path elimination
//  Block '<S208>/UpperRelop' : Unused code path elimination
//  Block '<S217>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S220>/Data Type Duplicate' : Unused code path elimination
//  Block '<S220>/Data Type Propagation' : Unused code path elimination
//  Block '<S218>/Data Type Duplicate' : Unused code path elimination
//  Block '<S218>/Data Type Propagation' : Unused code path elimination
//  Block '<S225>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S225>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S227>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S237>/Data Type Duplicate' : Unused code path elimination
//  Block '<S237>/Data Type Propagation' : Unused code path elimination
//  Block '<S228>/Delay Input2' : Unused code path elimination
//  Block '<S228>/Difference Inputs1' : Unused code path elimination
//  Block '<S228>/Difference Inputs2' : Unused code path elimination
//  Block '<S228>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S238>/Data Type Duplicate' : Unused code path elimination
//  Block '<S238>/Data Type Propagation' : Unused code path elimination
//  Block '<S238>/LowerRelop1' : Unused code path elimination
//  Block '<S238>/Switch' : Unused code path elimination
//  Block '<S238>/Switch2' : Unused code path elimination
//  Block '<S238>/UpperRelop' : Unused code path elimination
//  Block '<S228>/delta fall limit' : Unused code path elimination
//  Block '<S228>/delta rise limit' : Unused code path elimination
//  Block '<S228>/sample time' : Unused code path elimination
//  Block '<S229>/Delay Input2' : Unused code path elimination
//  Block '<S229>/Difference Inputs1' : Unused code path elimination
//  Block '<S229>/Difference Inputs2' : Unused code path elimination
//  Block '<S229>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S239>/Data Type Duplicate' : Unused code path elimination
//  Block '<S239>/Data Type Propagation' : Unused code path elimination
//  Block '<S239>/LowerRelop1' : Unused code path elimination
//  Block '<S239>/Switch' : Unused code path elimination
//  Block '<S239>/Switch2' : Unused code path elimination
//  Block '<S239>/UpperRelop' : Unused code path elimination
//  Block '<S229>/delta fall limit' : Unused code path elimination
//  Block '<S229>/delta rise limit' : Unused code path elimination
//  Block '<S229>/sample time' : Unused code path elimination
//  Block '<S230>/Data Type Duplicate' : Unused code path elimination
//  Block '<S230>/Data Type Propagation' : Unused code path elimination
//  Block '<S231>/Data Type Duplicate' : Unused code path elimination
//  Block '<S231>/Data Type Propagation' : Unused code path elimination
//  Block '<S231>/LowerRelop1' : Unused code path elimination
//  Block '<S231>/Switch' : Unused code path elimination
//  Block '<S231>/Switch2' : Unused code path elimination
//  Block '<S231>/UpperRelop' : Unused code path elimination
//  Block '<S232>/Data Type Duplicate' : Unused code path elimination
//  Block '<S232>/Data Type Propagation' : Unused code path elimination
//  Block '<S232>/LowerRelop1' : Unused code path elimination
//  Block '<S232>/Switch' : Unused code path elimination
//  Block '<S232>/Switch2' : Unused code path elimination
//  Block '<S232>/UpperRelop' : Unused code path elimination
//  Block '<S240>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S240>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S242>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S252>/Data Type Duplicate' : Unused code path elimination
//  Block '<S252>/Data Type Propagation' : Unused code path elimination
//  Block '<S243>/Delay Input2' : Unused code path elimination
//  Block '<S243>/Difference Inputs1' : Unused code path elimination
//  Block '<S243>/Difference Inputs2' : Unused code path elimination
//  Block '<S243>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S253>/Data Type Duplicate' : Unused code path elimination
//  Block '<S253>/Data Type Propagation' : Unused code path elimination
//  Block '<S253>/LowerRelop1' : Unused code path elimination
//  Block '<S253>/Switch' : Unused code path elimination
//  Block '<S253>/Switch2' : Unused code path elimination
//  Block '<S253>/UpperRelop' : Unused code path elimination
//  Block '<S243>/delta fall limit' : Unused code path elimination
//  Block '<S243>/delta rise limit' : Unused code path elimination
//  Block '<S243>/sample time' : Unused code path elimination
//  Block '<S244>/Delay Input2' : Unused code path elimination
//  Block '<S244>/Difference Inputs1' : Unused code path elimination
//  Block '<S244>/Difference Inputs2' : Unused code path elimination
//  Block '<S244>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S254>/Data Type Duplicate' : Unused code path elimination
//  Block '<S254>/Data Type Propagation' : Unused code path elimination
//  Block '<S254>/LowerRelop1' : Unused code path elimination
//  Block '<S254>/Switch' : Unused code path elimination
//  Block '<S254>/Switch2' : Unused code path elimination
//  Block '<S254>/UpperRelop' : Unused code path elimination
//  Block '<S244>/delta fall limit' : Unused code path elimination
//  Block '<S244>/delta rise limit' : Unused code path elimination
//  Block '<S244>/sample time' : Unused code path elimination
//  Block '<S245>/Data Type Duplicate' : Unused code path elimination
//  Block '<S245>/Data Type Propagation' : Unused code path elimination
//  Block '<S246>/Data Type Duplicate' : Unused code path elimination
//  Block '<S246>/Data Type Propagation' : Unused code path elimination
//  Block '<S246>/LowerRelop1' : Unused code path elimination
//  Block '<S246>/Switch' : Unused code path elimination
//  Block '<S246>/Switch2' : Unused code path elimination
//  Block '<S246>/UpperRelop' : Unused code path elimination
//  Block '<S247>/Data Type Duplicate' : Unused code path elimination
//  Block '<S247>/Data Type Propagation' : Unused code path elimination
//  Block '<S247>/LowerRelop1' : Unused code path elimination
//  Block '<S247>/Switch' : Unused code path elimination
//  Block '<S247>/Switch2' : Unused code path elimination
//  Block '<S247>/UpperRelop' : Unused code path elimination
//  Block '<S255>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S255>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S257>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S267>/Data Type Duplicate' : Unused code path elimination
//  Block '<S267>/Data Type Propagation' : Unused code path elimination
//  Block '<S258>/Delay Input2' : Unused code path elimination
//  Block '<S258>/Difference Inputs1' : Unused code path elimination
//  Block '<S258>/Difference Inputs2' : Unused code path elimination
//  Block '<S258>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S268>/Data Type Duplicate' : Unused code path elimination
//  Block '<S268>/Data Type Propagation' : Unused code path elimination
//  Block '<S268>/LowerRelop1' : Unused code path elimination
//  Block '<S268>/Switch' : Unused code path elimination
//  Block '<S268>/Switch2' : Unused code path elimination
//  Block '<S268>/UpperRelop' : Unused code path elimination
//  Block '<S258>/delta fall limit' : Unused code path elimination
//  Block '<S258>/delta rise limit' : Unused code path elimination
//  Block '<S258>/sample time' : Unused code path elimination
//  Block '<S259>/Delay Input2' : Unused code path elimination
//  Block '<S259>/Difference Inputs1' : Unused code path elimination
//  Block '<S259>/Difference Inputs2' : Unused code path elimination
//  Block '<S259>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S269>/Data Type Duplicate' : Unused code path elimination
//  Block '<S269>/Data Type Propagation' : Unused code path elimination
//  Block '<S269>/LowerRelop1' : Unused code path elimination
//  Block '<S269>/Switch' : Unused code path elimination
//  Block '<S269>/Switch2' : Unused code path elimination
//  Block '<S269>/UpperRelop' : Unused code path elimination
//  Block '<S259>/delta fall limit' : Unused code path elimination
//  Block '<S259>/delta rise limit' : Unused code path elimination
//  Block '<S259>/sample time' : Unused code path elimination
//  Block '<S260>/Data Type Duplicate' : Unused code path elimination
//  Block '<S260>/Data Type Propagation' : Unused code path elimination
//  Block '<S261>/Data Type Duplicate' : Unused code path elimination
//  Block '<S261>/Data Type Propagation' : Unused code path elimination
//  Block '<S261>/LowerRelop1' : Unused code path elimination
//  Block '<S261>/Switch' : Unused code path elimination
//  Block '<S261>/Switch2' : Unused code path elimination
//  Block '<S261>/UpperRelop' : Unused code path elimination
//  Block '<S262>/Data Type Duplicate' : Unused code path elimination
//  Block '<S262>/Data Type Propagation' : Unused code path elimination
//  Block '<S262>/LowerRelop1' : Unused code path elimination
//  Block '<S262>/Switch' : Unused code path elimination
//  Block '<S262>/Switch2' : Unused code path elimination
//  Block '<S262>/UpperRelop' : Unused code path elimination
//  Block '<S271>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S274>/Data Type Duplicate' : Unused code path elimination
//  Block '<S274>/Data Type Propagation' : Unused code path elimination
//  Block '<S272>/Data Type Duplicate' : Unused code path elimination
//  Block '<S272>/Data Type Propagation' : Unused code path elimination


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
//  '<S8>'   : 'fcsModel/Allocation/sysIdInjection'
//  '<S9>'   : 'fcsModel/Allocation/sysIdInjection/Compare To Constant'
//  '<S10>'  : 'fcsModel/Allocation/sysIdInjection/Compare To Constant1'
//  '<S11>'  : 'fcsModel/Allocation/sysIdInjection/Compare To Constant2'
//  '<S12>'  : 'fcsModel/Allocation/sysIdInjection/Compare To Constant3'
//  '<S13>'  : 'fcsModel/Allocation/sysIdInjection/Compare To Constant4'
//  '<S14>'  : 'fcsModel/Allocation/sysIdInjection/sysIdInputGeneration'
//  '<S15>'  : 'fcsModel/Allocation/sysIdInjection/sysIdInputGeneration/chirpInjection'
//  '<S16>'  : 'fcsModel/Allocation/sysIdInjection/sysIdInputGeneration/chirpInjection/Band-Limited White Noise'
//  '<S17>'  : 'fcsModel/Allocation/sysIdInjection/sysIdInputGeneration/chirpInjection/Generate Chirp'
//  '<S18>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller'
//  '<S19>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs'
//  '<S20>'  : 'fcsModel/Inner Loop Controller/Attitude Controller'
//  '<S21>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem'
//  '<S22>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S23>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S24>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug'
//  '<S25>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S26>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S27>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S28>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S29>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S30>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S31>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S32>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S33>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S34>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S35>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S36>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S37>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S38>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S39>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S40>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S41>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S42>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S43>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S44>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S45>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S46>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S47>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S48>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S49>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S50>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S51>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S52>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S53>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S54>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S55>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S56>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S57>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S58>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S59>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S60>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/Compare To Constant'
//  '<S61>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/EulerRates2BodyRates'
//  '<S62>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control'
//  '<S63>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant'
//  '<S64>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant1'
//  '<S65>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block'
//  '<S66>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1'
//  '<S67>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pickAttitudeCmdAndMeas'
//  '<S68>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug'
//  '<S69>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S70>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S71>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S72>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S73>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S74>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic'
//  '<S75>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic1'
//  '<S76>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Saturation Dynamic2'
//  '<S77>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S78>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S79>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S80>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S81>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S82>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S83>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S84>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S85>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S86>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S87>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S88>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S89>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic'
//  '<S90>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S91>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S92>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S93>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S94>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S95>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S96>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S97>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S98>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S99>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S100>' : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Rate Limiter Dynamic'
//  '<S101>' : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Saturation Dynamic'
//  '<S102>' : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S103>' : 'fcsModel/Inner Loop Controller/Attitude Controller/Attitude Control/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S104>' : 'fcsModel/Outer Loop Controller/Compare To Constant'
//  '<S105>' : 'fcsModel/Outer Loop Controller/Compare To Constant1'
//  '<S106>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl'
//  '<S107>' : 'fcsModel/Outer Loop Controller/assembleOuterLoopToInnerLoopBus'
//  '<S108>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs'
//  '<S109>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller'
//  '<S110>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller'
//  '<S111>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs/Compare To Constant1'
//  '<S112>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs/Compare To Constant2'
//  '<S113>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs'
//  '<S114>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control'
//  '<S115>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant'
//  '<S116>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant1'
//  '<S117>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant2'
//  '<S118>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant3'
//  '<S119>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant4'
//  '<S120>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant5'
//  '<S121>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant6'
//  '<S122>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter'
//  '<S123>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter1'
//  '<S124>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter2'
//  '<S125>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter/holdOutputAtCenter'
//  '<S126>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter1/holdOutputAtCenter'
//  '<S127>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter2/holdOutputAtCenter'
//  '<S128>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block'
//  '<S129>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1'
//  '<S130>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug'
//  '<S131>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S132>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S133>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S134>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S135>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S136>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic'
//  '<S137>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic1'
//  '<S138>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic2'
//  '<S139>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S140>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S141>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S142>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S143>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S144>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S145>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S146>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S147>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S148>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S149>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S150>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S151>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic'
//  '<S152>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S153>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S154>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S155>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S156>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S157>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S158>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S159>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S160>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S161>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S162>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic'
//  '<S163>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Saturation Dynamic'
//  '<S164>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S165>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S166>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs'
//  '<S167>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem'
//  '<S168>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant'
//  '<S169>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant1'
//  '<S170>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/Compare To Constant2'
//  '<S171>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/accelZKiSelectorVariantSubsystem'
//  '<S172>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/holdOutputAtCenter'
//  '<S173>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/hoverThrustVariantSubsystem'
//  '<S174>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/nedAccelToRollPitchCmd'
//  '<S175>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem'
//  '<S176>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/accelZKiSelectorVariantSubsystem/accelZCtrlKiPassThrough'
//  '<S177>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/holdOutputAtCenter/holdOutputAtCenter'
//  '<S178>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/hoverThrustVariantSubsystem/constantHoverThrust'
//  '<S179>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/nedAccelToRollPitchCmd/kinematicInversion'
//  '<S180>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/nedAccelToRollPitchCmd/kinematicInversion/NE Accel Cmds To Roll Pitch Cmds'
//  '<S181>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl'
//  '<S182>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller'
//  '<S183>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block'
//  '<S184>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1'
//  '<S185>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug'
//  '<S186>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S187>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S188>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S189>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S190>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S191>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic'
//  '<S192>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic1'
//  '<S193>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Saturation Dynamic2'
//  '<S194>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S195>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S196>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S197>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S198>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S199>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S200>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S201>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S202>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S203>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S204>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S205>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S206>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic'
//  '<S207>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S208>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S209>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S210>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S211>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S212>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S213>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S214>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S215>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S216>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S217>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Rate Limiter Dynamic'
//  '<S218>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Saturation Dynamic'
//  '<S219>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S220>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs/zAccelCtrlVariantSubsystem/zAccelCtrl/Z Accel Controller/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S221>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S222>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S223>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2'
//  '<S224>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug'
//  '<S225>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S226>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S227>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S228>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S229>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S230>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S231>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S232>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S233>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S234>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S235>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S236>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S237>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S238>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S239>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S240>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S241>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S242>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S243>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S244>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S245>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S246>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S247>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S248>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S249>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S250>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S251>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S252>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S253>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S254>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S255>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter'
//  '<S256>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter'
//  '<S257>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic'
//  '<S258>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic1'
//  '<S259>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic2'
//  '<S260>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic'
//  '<S261>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic1'
//  '<S262>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Saturation Dynamic2'
//  '<S263>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S264>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S265>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S266>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S267>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S268>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S269>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block2/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S270>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S271>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S272>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S273>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S274>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S275>' : 'fcsModel/RC Interpreter/Chart'
//  '<S276>' : 'fcsModel/RC Interpreter/Interpret RC In Cmds'


//-
//  Requirements for '<Root>': fcsModel

#endif                                 // RTW_HEADER_fcsModel_h_

//
// File trailer for generated code.
//
// [EOF]
//
