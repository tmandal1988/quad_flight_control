//
// File: fcsModel.h
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
#ifndef RTW_HEADER_fcsModel_h_
#define RTW_HEADER_fcsModel_h_
#include "rtwtypes.h"
#include "fcsModel_types.h"
#include <array>

// Class declaration for model fcsModel
class fcsModel final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for system '<S11>/pidWithDebug'
  struct DW_pidWithDebug_fcsModel_T {
    std::array<real_T, 2> num;
                      // '<S45>/Compute Deriv Filter Numerator And Denominator'
    std::array<real_T, 2> den;
                      // '<S45>/Compute Deriv Filter Numerator And Denominator'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S14>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S46>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S14>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S14>/Unit Delay1'
    real_T DiscreteTransferFcn_states; // '<S45>/Discrete Transfer Fcn'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S14>/Discrete-Time Integrator'
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
  struct DW_CoreSubsys_fcsModel_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S11>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S11>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S11>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S11>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<S10>/For Each Subsystem' 
  struct DW_CoreSubsys_fcsModel_i_T {
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S52>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S52>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S52>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S52>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW_fcsModel_T {
    std::array<DW_CoreSubsys_fcsModel_i_T, 3> CoreSubsys_p;// '<S10>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_T, 3> CoreSubsys;// '<S8>/For Each Subsystem' 
    std::array<real_T, 4> DiscreteTransferFcn_states;// '<S1>/Discrete Transfer Fcn' 
    real_T integratorReset;            // '<S4>/State Machine'
    int32_T durationCounter_1;         // '<S4>/State Machine'
    int32_T durationCounter_1_e;       // '<S4>/State Machine'
    uint16_T temporalCounter_i1;       // '<S4>/State Machine'
    uint8_T is_active_c4_rcInterpreter;// '<S4>/State Machine'
    uint8_T is_c4_rcInterpreter;       // '<S4>/State Machine'
  };

  // Constant parameters (default storage)
  struct ConstP_fcsModel_T {
    // Expression: rcParamsStruct
    //  Referenced by: '<S4>/State Machine'

    struct_ek64hRkWPGOrPAEJD7wTRH StateMachine_rcParamsStruct;

    // Computed Parameter: Constant_Value_h
    //  Referenced by: '<S3>/Constant'

    busOuterLoopToInnerLoop Constant_Value_h;

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
    busInnerLoopCtrlParams ctrlParams; // '<Root>/ctrlParams'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_fcsModel_T {
    std::array<real_T, 4> actuatorsCmds;// '<Root>/actuatorsCmds'
    busFcsDebug fcsDebug;              // '<Root>/fcsDebug'
  };

  // Copy Constructor
  fcsModel(fcsModel const&) = delete;

  // Assignment Operator
  fcsModel& operator= (fcsModel const&) & = delete;

  // Move Constructor
  fcsModel(fcsModel &&) = delete;

  // Move Assignment Operator
  fcsModel& operator= (fcsModel &&) = delete;

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
  static void initialize();

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

  // private member function(s) for subsystem '<S11>/pidWithDebug'
  static void fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, real_T rtu_integratorReset, const busPidParams
    *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
    busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_T *localDW);

  // private member function(s) for subsystem '<S30>/Compute Natural Frequency'
  static void fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps, real_T
    rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps);

  // private member function(s) for subsystem '<S11>/Signal Conditioning Block1'
  static void fcsMod_SignalConditioningBlock1(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW);

  // private member function(s) for subsystem '<Root>'
  boolean_T fcsModel_getStickStates(uint16_T rcInCmds_throttleCmd_nd, uint16_T
    rcInCmds_joystickYCmd_nd, uint16_T rcInCmds_joystickXCmd_nd, uint16_T
    rcInCmds_joystickZCmd_nd, const struct_ek64hRkWPGOrPAEJD7wTRH
    *b_rcParamsStruct);
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
//  Block '<S56>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S56>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S58>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S68>/Data Type Duplicate' : Unused code path elimination
//  Block '<S68>/Data Type Propagation' : Unused code path elimination
//  Block '<S59>/Delay Input2' : Unused code path elimination
//  Block '<S59>/Difference Inputs1' : Unused code path elimination
//  Block '<S59>/Difference Inputs2' : Unused code path elimination
//  Block '<S59>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S69>/Data Type Duplicate' : Unused code path elimination
//  Block '<S69>/Data Type Propagation' : Unused code path elimination
//  Block '<S69>/LowerRelop1' : Unused code path elimination
//  Block '<S69>/Switch' : Unused code path elimination
//  Block '<S69>/Switch2' : Unused code path elimination
//  Block '<S69>/UpperRelop' : Unused code path elimination
//  Block '<S59>/delta fall limit' : Unused code path elimination
//  Block '<S59>/delta rise limit' : Unused code path elimination
//  Block '<S59>/sample time' : Unused code path elimination
//  Block '<S60>/Delay Input2' : Unused code path elimination
//  Block '<S60>/Difference Inputs1' : Unused code path elimination
//  Block '<S60>/Difference Inputs2' : Unused code path elimination
//  Block '<S60>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S70>/Data Type Duplicate' : Unused code path elimination
//  Block '<S70>/Data Type Propagation' : Unused code path elimination
//  Block '<S70>/LowerRelop1' : Unused code path elimination
//  Block '<S70>/Switch' : Unused code path elimination
//  Block '<S70>/Switch2' : Unused code path elimination
//  Block '<S70>/UpperRelop' : Unused code path elimination
//  Block '<S60>/delta fall limit' : Unused code path elimination
//  Block '<S60>/delta rise limit' : Unused code path elimination
//  Block '<S60>/sample time' : Unused code path elimination
//  Block '<S61>/Data Type Duplicate' : Unused code path elimination
//  Block '<S61>/Data Type Propagation' : Unused code path elimination
//  Block '<S62>/Data Type Duplicate' : Unused code path elimination
//  Block '<S62>/Data Type Propagation' : Unused code path elimination
//  Block '<S62>/LowerRelop1' : Unused code path elimination
//  Block '<S62>/Switch' : Unused code path elimination
//  Block '<S62>/Switch2' : Unused code path elimination
//  Block '<S62>/UpperRelop' : Unused code path elimination
//  Block '<S63>/Data Type Duplicate' : Unused code path elimination
//  Block '<S63>/Data Type Propagation' : Unused code path elimination
//  Block '<S63>/LowerRelop1' : Unused code path elimination
//  Block '<S63>/Switch' : Unused code path elimination
//  Block '<S63>/Switch2' : Unused code path elimination
//  Block '<S63>/UpperRelop' : Unused code path elimination
//  Block '<S71>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S71>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S73>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S83>/Data Type Duplicate' : Unused code path elimination
//  Block '<S83>/Data Type Propagation' : Unused code path elimination
//  Block '<S74>/Delay Input2' : Unused code path elimination
//  Block '<S74>/Difference Inputs1' : Unused code path elimination
//  Block '<S74>/Difference Inputs2' : Unused code path elimination
//  Block '<S74>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S84>/Data Type Duplicate' : Unused code path elimination
//  Block '<S84>/Data Type Propagation' : Unused code path elimination
//  Block '<S84>/LowerRelop1' : Unused code path elimination
//  Block '<S84>/Switch' : Unused code path elimination
//  Block '<S84>/Switch2' : Unused code path elimination
//  Block '<S84>/UpperRelop' : Unused code path elimination
//  Block '<S74>/delta fall limit' : Unused code path elimination
//  Block '<S74>/delta rise limit' : Unused code path elimination
//  Block '<S74>/sample time' : Unused code path elimination
//  Block '<S75>/Delay Input2' : Unused code path elimination
//  Block '<S75>/Difference Inputs1' : Unused code path elimination
//  Block '<S75>/Difference Inputs2' : Unused code path elimination
//  Block '<S75>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S85>/Data Type Duplicate' : Unused code path elimination
//  Block '<S85>/Data Type Propagation' : Unused code path elimination
//  Block '<S85>/LowerRelop1' : Unused code path elimination
//  Block '<S85>/Switch' : Unused code path elimination
//  Block '<S85>/Switch2' : Unused code path elimination
//  Block '<S85>/UpperRelop' : Unused code path elimination
//  Block '<S75>/delta fall limit' : Unused code path elimination
//  Block '<S75>/delta rise limit' : Unused code path elimination
//  Block '<S75>/sample time' : Unused code path elimination
//  Block '<S76>/Data Type Duplicate' : Unused code path elimination
//  Block '<S76>/Data Type Propagation' : Unused code path elimination
//  Block '<S77>/Data Type Duplicate' : Unused code path elimination
//  Block '<S77>/Data Type Propagation' : Unused code path elimination
//  Block '<S77>/LowerRelop1' : Unused code path elimination
//  Block '<S77>/Switch' : Unused code path elimination
//  Block '<S77>/Switch2' : Unused code path elimination
//  Block '<S77>/UpperRelop' : Unused code path elimination
//  Block '<S78>/Data Type Duplicate' : Unused code path elimination
//  Block '<S78>/Data Type Propagation' : Unused code path elimination
//  Block '<S78>/LowerRelop1' : Unused code path elimination
//  Block '<S78>/Switch' : Unused code path elimination
//  Block '<S78>/Switch2' : Unused code path elimination
//  Block '<S78>/UpperRelop' : Unused code path elimination
//  Block '<S87>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S90>/Data Type Duplicate' : Unused code path elimination
//  Block '<S90>/Data Type Propagation' : Unused code path elimination
//  Block '<S88>/Data Type Duplicate' : Unused code path elimination
//  Block '<S88>/Data Type Propagation' : Unused code path elimination


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
//  '<S7>'   : 'fcsModel/Allocation/Compare To Constant1'
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
//  '<S52>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem'
//  '<S53>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S54>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S55>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug'
//  '<S56>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S57>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S58>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S59>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S60>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S61>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S62>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S63>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S64>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S65>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S66>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S67>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S68>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S69>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S70>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S71>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S72>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S73>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S74>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S75>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S76>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S77>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S78>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S79>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S80>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S81>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S82>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S83>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S84>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S85>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S86>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S87>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S88>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S89>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S90>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S91>'  : 'fcsModel/Outer Loop Controller/assembleOuterLoopToInnerLoopBus'
//  '<S92>'  : 'fcsModel/RC Interpreter/Interpret RC In Cmds'
//  '<S93>'  : 'fcsModel/RC Interpreter/State Machine'


//-
//  Requirements for '<Root>': fcsModel

#endif                                 // RTW_HEADER_fcsModel_h_

//
// File trailer for generated code.
//
// [EOF]
//
