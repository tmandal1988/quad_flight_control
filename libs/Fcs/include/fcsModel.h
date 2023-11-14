//
// File: fcsModel.h
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.87
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Nov 14 14:08:01 2023
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
extern const busOuterLoopToInnerLoop fcsModel_rtZbusOuterLoopToInnerLoop;// busOuterLoopToInnerLoop ground 

// Class declaration for model fcsModel
class fcsModel final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for system '<S13>/Discrete First Order Deriv Filter' 
  struct DW_DiscreteFirstOrderDerivFil_T {
    std::array<real_T, 2> num;
                      // '<S44>/Compute Deriv Filter Numerator And Denominator'
    std::array<real_T, 2> den;
                      // '<S44>/Compute Deriv Filter Numerator And Denominator'
    real_T DiscreteTransferFcn_states; // '<S44>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S10>/pidWithDebug'
  struct DW_pidWithDebug_fcsModel_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                   // '<S13>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S13>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S45>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S13>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S13>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S13>/Discrete-Time Integrator'
  };

  // Block signals and states (default storage) for system '<S10>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_f_T {
    std::array<real_T, 3> num;
                            // '<S30>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                            // '<S30>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S30>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S31>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S30>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S7>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S10>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S10>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S10>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S10>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S9>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_i_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S54>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S54>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S54>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S54>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S101>/pidWithDebug' 
  struct DW_pidWithDebug_fcsModel_i_T {
    DW_DiscreteFirstOrderDerivFil_T DiscreteFirstOrderDerivFilter;
                                  // '<S113>/Discrete First Order Deriv Filter'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S113>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S145>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S113>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S113>/Unit Delay1'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S113>/Discrete-Time Integrator'
  };

  // Block signals and states (default storage) for system '<S101>/Signal Conditioning Block1' 
  struct DW_SignalConditioningBlock1_g_T {
    std::array<real_T, 3> num;
                           // '<S130>/Compute Filter Numerator And Denominator'
    std::array<real_T, 3> den;
                           // '<S130>/Compute Filter Numerator And Denominator'
    std::array<real_T, 2> DiscreteTransferFcn_states;// '<S130>/Discrete Transfer Fcn' 
    real_T DelayInput2_DSTATE;         // '<S131>/Delay Input2'
    real_T DiscreteTransferFcn_tmp;    // '<S130>/Discrete Transfer Fcn'
  };

  // Block signals and states (default storage) for system '<S97>/NED Position Control' 
  struct DW_CoreSubsys_fcsModel_b_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S101>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S101>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S101>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S101>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S98>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_p_T {
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock;// '<S150>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_g_T SignalConditioningBlock1;// '<S150>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_i_T pidWithDebug;// '<S150>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S150>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW_fcsModel_T {
    std::array<DW_CoreSubsys_fcsModel_p_T, 3> CoreSubsys_i;// '<S98>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_b_T, 3> CoreSubsys_g;// '<S97>/NED Position Control' 
    std::array<DW_CoreSubsys_fcsModel_i_T, 3> CoreSubsys_p;// '<S9>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_T, 3> CoreSubsys;// '<S7>/For Each Subsystem' 
    busOuterLoopCtrlDebug RateTransition_Buffer0;// '<Root>/Rate Transition'
    busOuterLoopToInnerLoop Switch;    // '<S3>/Switch'
    busRcOutCmds rcOutCmds;            // '<S4>/Interpret RC In Cmds'
    std::array<real_T, 4> DiscreteTransferFcn_states;// '<S1>/Discrete Transfer Fcn' 
    std::array<real_T, 2> DiscreteTransferFcn_states_e;// '<S106>/Discrete Transfer Fcn' 
    real_T last_input;                 // '<S105>/holdOutputAtCenter'
    real_T pers_base_alt_point_m;      // '<S104>/computeAltFromBaro'
    int32_T durationCounter_1;         // '<S4>/Chart'
    int32_T durationCounter_1_j;       // '<S4>/Chart'
    uint16_T temporalCounter_i1;       // '<S4>/Chart'
    uint8_T is_active_c1_rcInterpreter;// '<S4>/Chart'
    uint8_T is_c1_rcInterpreter;       // '<S4>/Chart'
    boolean_T rcCheckFlag;             // '<S4>/Chart'
  };

  // Constant parameters (default storage)
  struct ConstP_fcsModel_T {
    // Expression: allocationDataStruct.allocationMatrix
    //  Referenced by: '<S1>/Constant'

    std::array<real_T, 16> Constant_Value_c;

    // Expression: vehicleConstants.inertia_kgm2
    //  Referenced by: '<S2>/Constant'

    std::array<real_T, 9> Constant_Value_n;
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

  // private member function(s) for subsystem '<S13>/Discrete First Order Deriv Filter'
  static void f_DiscreteFirstOrderDerivFilter(real_T rtu_input, real_T
    rtu_filterBandwidth_radps, real_T *rty_filteredInputRate, real_T
    rtp_sampleTime_s, DW_DiscreteFirstOrderDerivFil_T *localDW);

  // private member function(s) for subsystem '<S10>/pidWithDebug'
  static void fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
    *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
    busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_T *localDW);

  // private member function(s) for subsystem '<S29>/Compute Natural Frequency'
  static void fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps, real_T
    rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps);

  // private member function(s) for subsystem '<S29>/Compute Numerator And Denominator'
  static void ComputeNumeratorAndDenominator(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_rateNum[3], real_T rty_accelNum[3],
    real_T rty_den[3], real_T rtp_sampleTime_s);

  // private member function(s) for subsystem '<S30>/Compute Filter Numerator And Denominator'
  static void ComputeFilterNumeratorAndD_Init(real_T rty_num[3], real_T rty_den
    [3]);
  static void ComputeFilterNumeratorAndDenomi(real_T rtu_naturalFrequency_radps,
    real_T rtu_dampingRatio_nd, real_T rty_num[3], real_T rty_den[3], real_T
    rtp_sampleTime_s);

  // private member function(s) for subsystem '<S10>/Signal Conditioning Block1'
  static void f_SignalConditioningBlock1_Init(DW_SignalConditioningBlock1_f_T
    *localDW);
  static void fcsMod_SignalConditioningBlock1(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW);

  // private member function(s) for subsystem '<S101>/pidWithDebug'
  static void fcsModel_pidWithDebug_j(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
    *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
    busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_i_T *localDW);

  // private member function(s) for subsystem '<S101>/Signal Conditioning Block1'
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
//  Block '<S14>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S14>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S16>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S26>/Data Type Duplicate' : Unused code path elimination
//  Block '<S26>/Data Type Propagation' : Unused code path elimination
//  Block '<S17>/Delay Input2' : Unused code path elimination
//  Block '<S17>/Difference Inputs1' : Unused code path elimination
//  Block '<S17>/Difference Inputs2' : Unused code path elimination
//  Block '<S17>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S27>/Data Type Duplicate' : Unused code path elimination
//  Block '<S27>/Data Type Propagation' : Unused code path elimination
//  Block '<S27>/LowerRelop1' : Unused code path elimination
//  Block '<S27>/Switch' : Unused code path elimination
//  Block '<S27>/Switch2' : Unused code path elimination
//  Block '<S27>/UpperRelop' : Unused code path elimination
//  Block '<S17>/delta fall limit' : Unused code path elimination
//  Block '<S17>/delta rise limit' : Unused code path elimination
//  Block '<S17>/sample time' : Unused code path elimination
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
//  Block '<S19>/Data Type Duplicate' : Unused code path elimination
//  Block '<S19>/Data Type Propagation' : Unused code path elimination
//  Block '<S20>/Data Type Duplicate' : Unused code path elimination
//  Block '<S20>/Data Type Propagation' : Unused code path elimination
//  Block '<S20>/LowerRelop1' : Unused code path elimination
//  Block '<S20>/Switch' : Unused code path elimination
//  Block '<S20>/Switch2' : Unused code path elimination
//  Block '<S20>/UpperRelop' : Unused code path elimination
//  Block '<S21>/Data Type Duplicate' : Unused code path elimination
//  Block '<S21>/Data Type Propagation' : Unused code path elimination
//  Block '<S21>/LowerRelop1' : Unused code path elimination
//  Block '<S21>/Switch' : Unused code path elimination
//  Block '<S21>/Switch2' : Unused code path elimination
//  Block '<S21>/UpperRelop' : Unused code path elimination
//  Block '<S29>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S29>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S31>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S41>/Data Type Duplicate' : Unused code path elimination
//  Block '<S41>/Data Type Propagation' : Unused code path elimination
//  Block '<S32>/Delay Input2' : Unused code path elimination
//  Block '<S32>/Difference Inputs1' : Unused code path elimination
//  Block '<S32>/Difference Inputs2' : Unused code path elimination
//  Block '<S32>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S42>/Data Type Duplicate' : Unused code path elimination
//  Block '<S42>/Data Type Propagation' : Unused code path elimination
//  Block '<S42>/LowerRelop1' : Unused code path elimination
//  Block '<S42>/Switch' : Unused code path elimination
//  Block '<S42>/Switch2' : Unused code path elimination
//  Block '<S42>/UpperRelop' : Unused code path elimination
//  Block '<S32>/delta fall limit' : Unused code path elimination
//  Block '<S32>/delta rise limit' : Unused code path elimination
//  Block '<S32>/sample time' : Unused code path elimination
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
//  Block '<S34>/Data Type Duplicate' : Unused code path elimination
//  Block '<S34>/Data Type Propagation' : Unused code path elimination
//  Block '<S35>/Data Type Duplicate' : Unused code path elimination
//  Block '<S35>/Data Type Propagation' : Unused code path elimination
//  Block '<S35>/LowerRelop1' : Unused code path elimination
//  Block '<S35>/Switch' : Unused code path elimination
//  Block '<S35>/Switch2' : Unused code path elimination
//  Block '<S35>/UpperRelop' : Unused code path elimination
//  Block '<S36>/Data Type Duplicate' : Unused code path elimination
//  Block '<S36>/Data Type Propagation' : Unused code path elimination
//  Block '<S36>/LowerRelop1' : Unused code path elimination
//  Block '<S36>/Switch' : Unused code path elimination
//  Block '<S36>/Switch2' : Unused code path elimination
//  Block '<S36>/UpperRelop' : Unused code path elimination
//  Block '<S45>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S48>/Data Type Duplicate' : Unused code path elimination
//  Block '<S48>/Data Type Propagation' : Unused code path elimination
//  Block '<S46>/Data Type Duplicate' : Unused code path elimination
//  Block '<S46>/Data Type Propagation' : Unused code path elimination
//  Block '<S58>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S58>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S60>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S70>/Data Type Duplicate' : Unused code path elimination
//  Block '<S70>/Data Type Propagation' : Unused code path elimination
//  Block '<S61>/Delay Input2' : Unused code path elimination
//  Block '<S61>/Difference Inputs1' : Unused code path elimination
//  Block '<S61>/Difference Inputs2' : Unused code path elimination
//  Block '<S61>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S71>/Data Type Duplicate' : Unused code path elimination
//  Block '<S71>/Data Type Propagation' : Unused code path elimination
//  Block '<S71>/LowerRelop1' : Unused code path elimination
//  Block '<S71>/Switch' : Unused code path elimination
//  Block '<S71>/Switch2' : Unused code path elimination
//  Block '<S71>/UpperRelop' : Unused code path elimination
//  Block '<S61>/delta fall limit' : Unused code path elimination
//  Block '<S61>/delta rise limit' : Unused code path elimination
//  Block '<S61>/sample time' : Unused code path elimination
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
//  Block '<S63>/Data Type Duplicate' : Unused code path elimination
//  Block '<S63>/Data Type Propagation' : Unused code path elimination
//  Block '<S64>/Data Type Duplicate' : Unused code path elimination
//  Block '<S64>/Data Type Propagation' : Unused code path elimination
//  Block '<S64>/LowerRelop1' : Unused code path elimination
//  Block '<S64>/Switch' : Unused code path elimination
//  Block '<S64>/Switch2' : Unused code path elimination
//  Block '<S64>/UpperRelop' : Unused code path elimination
//  Block '<S65>/Data Type Duplicate' : Unused code path elimination
//  Block '<S65>/Data Type Propagation' : Unused code path elimination
//  Block '<S65>/LowerRelop1' : Unused code path elimination
//  Block '<S65>/Switch' : Unused code path elimination
//  Block '<S65>/Switch2' : Unused code path elimination
//  Block '<S65>/UpperRelop' : Unused code path elimination
//  Block '<S73>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S73>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S75>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S85>/Data Type Duplicate' : Unused code path elimination
//  Block '<S85>/Data Type Propagation' : Unused code path elimination
//  Block '<S76>/Delay Input2' : Unused code path elimination
//  Block '<S76>/Difference Inputs1' : Unused code path elimination
//  Block '<S76>/Difference Inputs2' : Unused code path elimination
//  Block '<S76>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S86>/Data Type Duplicate' : Unused code path elimination
//  Block '<S86>/Data Type Propagation' : Unused code path elimination
//  Block '<S86>/LowerRelop1' : Unused code path elimination
//  Block '<S86>/Switch' : Unused code path elimination
//  Block '<S86>/Switch2' : Unused code path elimination
//  Block '<S86>/UpperRelop' : Unused code path elimination
//  Block '<S76>/delta fall limit' : Unused code path elimination
//  Block '<S76>/delta rise limit' : Unused code path elimination
//  Block '<S76>/sample time' : Unused code path elimination
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
//  Block '<S78>/Data Type Duplicate' : Unused code path elimination
//  Block '<S78>/Data Type Propagation' : Unused code path elimination
//  Block '<S79>/Data Type Duplicate' : Unused code path elimination
//  Block '<S79>/Data Type Propagation' : Unused code path elimination
//  Block '<S79>/LowerRelop1' : Unused code path elimination
//  Block '<S79>/Switch' : Unused code path elimination
//  Block '<S79>/Switch2' : Unused code path elimination
//  Block '<S79>/UpperRelop' : Unused code path elimination
//  Block '<S80>/Data Type Duplicate' : Unused code path elimination
//  Block '<S80>/Data Type Propagation' : Unused code path elimination
//  Block '<S80>/LowerRelop1' : Unused code path elimination
//  Block '<S80>/Switch' : Unused code path elimination
//  Block '<S80>/Switch2' : Unused code path elimination
//  Block '<S80>/UpperRelop' : Unused code path elimination
//  Block '<S89>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S92>/Data Type Duplicate' : Unused code path elimination
//  Block '<S92>/Data Type Propagation' : Unused code path elimination
//  Block '<S90>/Data Type Duplicate' : Unused code path elimination
//  Block '<S90>/Data Type Propagation' : Unused code path elimination
//  Block '<S114>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S114>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S116>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S126>/Data Type Duplicate' : Unused code path elimination
//  Block '<S126>/Data Type Propagation' : Unused code path elimination
//  Block '<S117>/Delay Input2' : Unused code path elimination
//  Block '<S117>/Difference Inputs1' : Unused code path elimination
//  Block '<S117>/Difference Inputs2' : Unused code path elimination
//  Block '<S117>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S127>/Data Type Duplicate' : Unused code path elimination
//  Block '<S127>/Data Type Propagation' : Unused code path elimination
//  Block '<S127>/LowerRelop1' : Unused code path elimination
//  Block '<S127>/Switch' : Unused code path elimination
//  Block '<S127>/Switch2' : Unused code path elimination
//  Block '<S127>/UpperRelop' : Unused code path elimination
//  Block '<S117>/delta fall limit' : Unused code path elimination
//  Block '<S117>/delta rise limit' : Unused code path elimination
//  Block '<S117>/sample time' : Unused code path elimination
//  Block '<S118>/Delay Input2' : Unused code path elimination
//  Block '<S118>/Difference Inputs1' : Unused code path elimination
//  Block '<S118>/Difference Inputs2' : Unused code path elimination
//  Block '<S118>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S128>/Data Type Duplicate' : Unused code path elimination
//  Block '<S128>/Data Type Propagation' : Unused code path elimination
//  Block '<S128>/LowerRelop1' : Unused code path elimination
//  Block '<S128>/Switch' : Unused code path elimination
//  Block '<S128>/Switch2' : Unused code path elimination
//  Block '<S128>/UpperRelop' : Unused code path elimination
//  Block '<S118>/delta fall limit' : Unused code path elimination
//  Block '<S118>/delta rise limit' : Unused code path elimination
//  Block '<S118>/sample time' : Unused code path elimination
//  Block '<S119>/Data Type Duplicate' : Unused code path elimination
//  Block '<S119>/Data Type Propagation' : Unused code path elimination
//  Block '<S120>/Data Type Duplicate' : Unused code path elimination
//  Block '<S120>/Data Type Propagation' : Unused code path elimination
//  Block '<S120>/LowerRelop1' : Unused code path elimination
//  Block '<S120>/Switch' : Unused code path elimination
//  Block '<S120>/Switch2' : Unused code path elimination
//  Block '<S120>/UpperRelop' : Unused code path elimination
//  Block '<S121>/Data Type Duplicate' : Unused code path elimination
//  Block '<S121>/Data Type Propagation' : Unused code path elimination
//  Block '<S121>/LowerRelop1' : Unused code path elimination
//  Block '<S121>/Switch' : Unused code path elimination
//  Block '<S121>/Switch2' : Unused code path elimination
//  Block '<S121>/UpperRelop' : Unused code path elimination
//  Block '<S129>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S129>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S131>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S141>/Data Type Duplicate' : Unused code path elimination
//  Block '<S141>/Data Type Propagation' : Unused code path elimination
//  Block '<S132>/Delay Input2' : Unused code path elimination
//  Block '<S132>/Difference Inputs1' : Unused code path elimination
//  Block '<S132>/Difference Inputs2' : Unused code path elimination
//  Block '<S132>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S142>/Data Type Duplicate' : Unused code path elimination
//  Block '<S142>/Data Type Propagation' : Unused code path elimination
//  Block '<S142>/LowerRelop1' : Unused code path elimination
//  Block '<S142>/Switch' : Unused code path elimination
//  Block '<S142>/Switch2' : Unused code path elimination
//  Block '<S142>/UpperRelop' : Unused code path elimination
//  Block '<S132>/delta fall limit' : Unused code path elimination
//  Block '<S132>/delta rise limit' : Unused code path elimination
//  Block '<S132>/sample time' : Unused code path elimination
//  Block '<S133>/Delay Input2' : Unused code path elimination
//  Block '<S133>/Difference Inputs1' : Unused code path elimination
//  Block '<S133>/Difference Inputs2' : Unused code path elimination
//  Block '<S133>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Propagation' : Unused code path elimination
//  Block '<S143>/LowerRelop1' : Unused code path elimination
//  Block '<S143>/Switch' : Unused code path elimination
//  Block '<S143>/Switch2' : Unused code path elimination
//  Block '<S143>/UpperRelop' : Unused code path elimination
//  Block '<S133>/delta fall limit' : Unused code path elimination
//  Block '<S133>/delta rise limit' : Unused code path elimination
//  Block '<S133>/sample time' : Unused code path elimination
//  Block '<S134>/Data Type Duplicate' : Unused code path elimination
//  Block '<S134>/Data Type Propagation' : Unused code path elimination
//  Block '<S135>/Data Type Duplicate' : Unused code path elimination
//  Block '<S135>/Data Type Propagation' : Unused code path elimination
//  Block '<S135>/LowerRelop1' : Unused code path elimination
//  Block '<S135>/Switch' : Unused code path elimination
//  Block '<S135>/Switch2' : Unused code path elimination
//  Block '<S135>/UpperRelop' : Unused code path elimination
//  Block '<S136>/Data Type Duplicate' : Unused code path elimination
//  Block '<S136>/Data Type Propagation' : Unused code path elimination
//  Block '<S136>/LowerRelop1' : Unused code path elimination
//  Block '<S136>/Switch' : Unused code path elimination
//  Block '<S136>/Switch2' : Unused code path elimination
//  Block '<S136>/UpperRelop' : Unused code path elimination
//  Block '<S145>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S148>/Data Type Duplicate' : Unused code path elimination
//  Block '<S148>/Data Type Propagation' : Unused code path elimination
//  Block '<S146>/Data Type Duplicate' : Unused code path elimination
//  Block '<S146>/Data Type Propagation' : Unused code path elimination
//  Block '<S154>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S154>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S156>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S166>/Data Type Duplicate' : Unused code path elimination
//  Block '<S166>/Data Type Propagation' : Unused code path elimination
//  Block '<S157>/Delay Input2' : Unused code path elimination
//  Block '<S157>/Difference Inputs1' : Unused code path elimination
//  Block '<S157>/Difference Inputs2' : Unused code path elimination
//  Block '<S157>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S167>/Data Type Duplicate' : Unused code path elimination
//  Block '<S167>/Data Type Propagation' : Unused code path elimination
//  Block '<S167>/LowerRelop1' : Unused code path elimination
//  Block '<S167>/Switch' : Unused code path elimination
//  Block '<S167>/Switch2' : Unused code path elimination
//  Block '<S167>/UpperRelop' : Unused code path elimination
//  Block '<S157>/delta fall limit' : Unused code path elimination
//  Block '<S157>/delta rise limit' : Unused code path elimination
//  Block '<S157>/sample time' : Unused code path elimination
//  Block '<S158>/Delay Input2' : Unused code path elimination
//  Block '<S158>/Difference Inputs1' : Unused code path elimination
//  Block '<S158>/Difference Inputs2' : Unused code path elimination
//  Block '<S158>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S168>/Data Type Duplicate' : Unused code path elimination
//  Block '<S168>/Data Type Propagation' : Unused code path elimination
//  Block '<S168>/LowerRelop1' : Unused code path elimination
//  Block '<S168>/Switch' : Unused code path elimination
//  Block '<S168>/Switch2' : Unused code path elimination
//  Block '<S168>/UpperRelop' : Unused code path elimination
//  Block '<S158>/delta fall limit' : Unused code path elimination
//  Block '<S158>/delta rise limit' : Unused code path elimination
//  Block '<S158>/sample time' : Unused code path elimination
//  Block '<S159>/Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Propagation' : Unused code path elimination
//  Block '<S160>/Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Propagation' : Unused code path elimination
//  Block '<S160>/LowerRelop1' : Unused code path elimination
//  Block '<S160>/Switch' : Unused code path elimination
//  Block '<S160>/Switch2' : Unused code path elimination
//  Block '<S160>/UpperRelop' : Unused code path elimination
//  Block '<S161>/Data Type Duplicate' : Unused code path elimination
//  Block '<S161>/Data Type Propagation' : Unused code path elimination
//  Block '<S161>/LowerRelop1' : Unused code path elimination
//  Block '<S161>/Switch' : Unused code path elimination
//  Block '<S161>/Switch2' : Unused code path elimination
//  Block '<S161>/UpperRelop' : Unused code path elimination
//  Block '<S169>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S169>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S171>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S181>/Data Type Duplicate' : Unused code path elimination
//  Block '<S181>/Data Type Propagation' : Unused code path elimination
//  Block '<S172>/Delay Input2' : Unused code path elimination
//  Block '<S172>/Difference Inputs1' : Unused code path elimination
//  Block '<S172>/Difference Inputs2' : Unused code path elimination
//  Block '<S172>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S182>/Data Type Duplicate' : Unused code path elimination
//  Block '<S182>/Data Type Propagation' : Unused code path elimination
//  Block '<S182>/LowerRelop1' : Unused code path elimination
//  Block '<S182>/Switch' : Unused code path elimination
//  Block '<S182>/Switch2' : Unused code path elimination
//  Block '<S182>/UpperRelop' : Unused code path elimination
//  Block '<S172>/delta fall limit' : Unused code path elimination
//  Block '<S172>/delta rise limit' : Unused code path elimination
//  Block '<S172>/sample time' : Unused code path elimination
//  Block '<S173>/Delay Input2' : Unused code path elimination
//  Block '<S173>/Difference Inputs1' : Unused code path elimination
//  Block '<S173>/Difference Inputs2' : Unused code path elimination
//  Block '<S173>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S183>/Data Type Duplicate' : Unused code path elimination
//  Block '<S183>/Data Type Propagation' : Unused code path elimination
//  Block '<S183>/LowerRelop1' : Unused code path elimination
//  Block '<S183>/Switch' : Unused code path elimination
//  Block '<S183>/Switch2' : Unused code path elimination
//  Block '<S183>/UpperRelop' : Unused code path elimination
//  Block '<S173>/delta fall limit' : Unused code path elimination
//  Block '<S173>/delta rise limit' : Unused code path elimination
//  Block '<S173>/sample time' : Unused code path elimination
//  Block '<S174>/Data Type Duplicate' : Unused code path elimination
//  Block '<S174>/Data Type Propagation' : Unused code path elimination
//  Block '<S175>/Data Type Duplicate' : Unused code path elimination
//  Block '<S175>/Data Type Propagation' : Unused code path elimination
//  Block '<S175>/LowerRelop1' : Unused code path elimination
//  Block '<S175>/Switch' : Unused code path elimination
//  Block '<S175>/Switch2' : Unused code path elimination
//  Block '<S175>/UpperRelop' : Unused code path elimination
//  Block '<S176>/Data Type Duplicate' : Unused code path elimination
//  Block '<S176>/Data Type Propagation' : Unused code path elimination
//  Block '<S176>/LowerRelop1' : Unused code path elimination
//  Block '<S176>/Switch' : Unused code path elimination
//  Block '<S176>/Switch2' : Unused code path elimination
//  Block '<S176>/UpperRelop' : Unused code path elimination
//  Block '<S185>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S188>/Data Type Duplicate' : Unused code path elimination
//  Block '<S188>/Data Type Propagation' : Unused code path elimination
//  Block '<S186>/Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Propagation' : Unused code path elimination


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
//  '<S7>'   : 'fcsModel/Inner Loop Controller/Angular Rate Controller'
//  '<S8>'   : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs'
//  '<S9>'   : 'fcsModel/Inner Loop Controller/Attitude Controller'
//  '<S10>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem'
//  '<S11>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S12>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S13>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug'
//  '<S14>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S15>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S16>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S17>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S18>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S19>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S20>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S21>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S22>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S23>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S24>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S25>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S26>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S27>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S28>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S29>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S30>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S31>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S32>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S33>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S34>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S35>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S36>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S37>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S38>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S39>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S40>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S41>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S42>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S43>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S44>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S45>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S46>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S47>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S48>'  : 'fcsModel/Inner Loop Controller/Angular Rate Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S49>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/Compare To Constant'
//  '<S50>'  : 'fcsModel/Inner Loop Controller/Assemble Angular Rate Ctrl Inputs/EulerRates2BodyRates'
//  '<S51>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant'
//  '<S52>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant1'
//  '<S53>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/Compare To Constant2'
//  '<S54>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem'
//  '<S55>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S56>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S57>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug'
//  '<S58>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S59>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S60>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S61>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S62>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S63>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S64>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S65>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S66>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S67>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S68>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S69>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S70>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S71>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S72>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S73>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S74>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S75>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S76>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S77>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S78>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S79>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S80>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S81>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S82>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S83>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S84>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S85>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S86>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S87>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S88>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S89>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S90>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S91>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S92>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S93>'  : 'fcsModel/Outer Loop Controller/Compare To Constant'
//  '<S94>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl'
//  '<S95>'  : 'fcsModel/Outer Loop Controller/assembleOuterLoopToInnerLoopBus'
//  '<S96>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs'
//  '<S97>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller'
//  '<S98>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller'
//  '<S99>'  : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Assemble Vel Ctrl Inputs/Compare To Constant1'
//  '<S100>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs'
//  '<S101>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control'
//  '<S102>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant'
//  '<S103>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/Compare To Constant1'
//  '<S104>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/baroAltEstimate'
//  '<S105>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter'
//  '<S106>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/baroAltEstimate/Discrete Second Order Filter'
//  '<S107>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/baroAltEstimate/computeAltFromBaro'
//  '<S108>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/baroAltEstimate/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S109>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/baroAltEstimate/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S110>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/Assemble Position Controller Inputs/holdOutputAtCenter/holdOutputAtCenter'
//  '<S111>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block'
//  '<S112>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1'
//  '<S113>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug'
//  '<S114>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S115>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S116>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S117>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S118>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S119>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic'
//  '<S120>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic1'
//  '<S121>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Saturation Dynamic2'
//  '<S122>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S123>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S124>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S125>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S126>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S127>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S128>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S129>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S130>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S131>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S132>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S133>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S134>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic'
//  '<S135>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S136>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S137>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S138>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S139>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S140>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S141>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S142>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S143>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S144>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S145>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic'
//  '<S146>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Saturation Dynamic'
//  '<S147>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S148>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Position Controller/NED Position Control/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S149>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/Assemble Inner Loop Inputs'
//  '<S150>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem'
//  '<S151>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S152>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S153>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug'
//  '<S154>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S155>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S156>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S157>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S158>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S159>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S160>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S161>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S162>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S163>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S164>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S165>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S166>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S167>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S168>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S169>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S170>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S171>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S172>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S173>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S174>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S175>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S176>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S177>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S178>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S179>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S180>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S181>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S182>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S183>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S184>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S185>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S186>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S187>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S188>' : 'fcsModel/Outer Loop Controller/PosAndVelCtrl/Velocity Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S189>' : 'fcsModel/RC Interpreter/Chart'
//  '<S190>' : 'fcsModel/RC Interpreter/Interpret RC In Cmds'


//-
//  Requirements for '<Root>': fcsModel

#endif                                 // RTW_HEADER_fcsModel_h_

//
// File trailer for generated code.
//
// [EOF]
//
