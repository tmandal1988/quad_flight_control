//
// File: fcsModel.h
//
// Code generated for Simulink model 'fcsModel'.
//
// Model version                  : 1.67
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Fri Oct 27 23:29:18 2023
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
  // Block signals and states (default storage) for system '<S10>/pidWithDebug'
  struct DW_pidWithDebug_fcsModel_T {
    std::array<real_T, 2> num;
                      // '<S44>/Compute Deriv Filter Numerator And Denominator'
    std::array<real_T, 2> den;
                      // '<S44>/Compute Deriv Filter Numerator And Denominator'
    real_T DiscreteTimeIntegrator_DSTATE;// '<S13>/Discrete-Time Integrator'
    real_T DelayInput2_DSTATE;         // '<S45>/Delay Input2'
    real_T UnitDelay_DSTATE;           // '<S13>/Unit Delay'
    real_T UnitDelay1_DSTATE;          // '<S13>/Unit Delay1'
    real_T DiscreteTransferFcn_states; // '<S44>/Discrete Transfer Fcn'
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
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock;// '<S53>/Signal Conditioning Block' 
    DW_SignalConditioningBlock1_f_T SignalConditioningBlock1;// '<S53>/Signal Conditioning Block1' 
    DW_pidWithDebug_fcsModel_T pidWithDebug;// '<S53>/pidWithDebug'
    real_T UnitDelay_DSTATE;           // '<S53>/Unit Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW_fcsModel_T {
    std::array<DW_CoreSubsys_fcsModel_i_T, 3> CoreSubsys_p;// '<S9>/For Each Subsystem' 
    std::array<DW_CoreSubsys_fcsModel_T, 3> CoreSubsys;// '<S7>/For Each Subsystem' 
    std::array<real_T, 4> DiscreteTransferFcn_states;// '<S1>/Discrete Transfer Fcn' 
    int32_T durationCounter_1;         // '<S4>/Chart'
    int32_T durationCounter_1_j;       // '<S4>/Chart'
    uint16_T temporalCounter_i1;       // '<S4>/Chart'
    uint8_T is_active_c1_rcInterpreter;// '<S4>/Chart'
    uint8_T is_c1_rcInterpreter;       // '<S4>/Chart'
    boolean_T rcCheckFlag;             // '<S4>/Chart'
  };

  // Constant parameters (default storage)
  struct ConstP_fcsModel_T {
    // Computed Parameter: Constant_Value
    //  Referenced by: '<S3>/Constant'

    busOuterLoopToInnerLoop Constant_Value;

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

  // private member function(s) for subsystem '<S10>/pidWithDebug'
  static void fcsModel_pidWithDebug(real_T rtu_feedForward, real_T rtu_cmd,
    real_T rtu_meas, boolean_T rtu_integratorReset, const busPidParams
    *rtu_pidParamBus, real_T rtu_trackingCtrlCmd, real_T *rty_ctrlCmd,
    busPidDebug *rty_pidDebug, real_T rtp_sampleTime_s,
    DW_pidWithDebug_fcsModel_T *localDW);

  // private member function(s) for subsystem '<S29>/Compute Natural Frequency'
  static void fcsMode_ComputeNaturalFrequency(real_T rtu_bandwidth_radps, real_T
    rtu_dampingRatio_nd, real_T *rty_naturalFrequency_radps);

  // private member function(s) for subsystem '<S10>/Signal Conditioning Block1'
  static void fcsMod_SignalConditioningBlock1(real_T rtu_input, const
    busSignalConditioningParams *rtu_params, real_T *rty_filteredInput, real_T
    rtp_sampleTime_s, DW_SignalConditioningBlock1_f_T *localDW);

  // private member function(s) for subsystem '<Root>'
  boolean_T fcsModel_checkRcCmds(const busRcInCmds
    *BusConversion_InsertedFor_Chart);
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
//  Block '<S57>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S57>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S59>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S69>/Data Type Duplicate' : Unused code path elimination
//  Block '<S69>/Data Type Propagation' : Unused code path elimination
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
//  Block '<S62>/Data Type Duplicate' : Unused code path elimination
//  Block '<S62>/Data Type Propagation' : Unused code path elimination
//  Block '<S63>/Data Type Duplicate' : Unused code path elimination
//  Block '<S63>/Data Type Propagation' : Unused code path elimination
//  Block '<S63>/LowerRelop1' : Unused code path elimination
//  Block '<S63>/Switch' : Unused code path elimination
//  Block '<S63>/Switch2' : Unused code path elimination
//  Block '<S63>/UpperRelop' : Unused code path elimination
//  Block '<S64>/Data Type Duplicate' : Unused code path elimination
//  Block '<S64>/Data Type Propagation' : Unused code path elimination
//  Block '<S64>/LowerRelop1' : Unused code path elimination
//  Block '<S64>/Switch' : Unused code path elimination
//  Block '<S64>/Switch2' : Unused code path elimination
//  Block '<S64>/UpperRelop' : Unused code path elimination
//  Block '<S72>/Discrete Transfer Fcn' : Unused code path elimination
//  Block '<S72>/Discrete Transfer Fcn1' : Unused code path elimination
//  Block '<S74>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S84>/Data Type Duplicate' : Unused code path elimination
//  Block '<S84>/Data Type Propagation' : Unused code path elimination
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
//  Block '<S77>/Data Type Duplicate' : Unused code path elimination
//  Block '<S77>/Data Type Propagation' : Unused code path elimination
//  Block '<S78>/Data Type Duplicate' : Unused code path elimination
//  Block '<S78>/Data Type Propagation' : Unused code path elimination
//  Block '<S78>/LowerRelop1' : Unused code path elimination
//  Block '<S78>/Switch' : Unused code path elimination
//  Block '<S78>/Switch2' : Unused code path elimination
//  Block '<S78>/UpperRelop' : Unused code path elimination
//  Block '<S79>/Data Type Duplicate' : Unused code path elimination
//  Block '<S79>/Data Type Propagation' : Unused code path elimination
//  Block '<S79>/LowerRelop1' : Unused code path elimination
//  Block '<S79>/Switch' : Unused code path elimination
//  Block '<S79>/Switch2' : Unused code path elimination
//  Block '<S79>/UpperRelop' : Unused code path elimination
//  Block '<S88>/FixPt Data Type Duplicate' : Unused code path elimination
//  Block '<S91>/Data Type Duplicate' : Unused code path elimination
//  Block '<S91>/Data Type Propagation' : Unused code path elimination
//  Block '<S89>/Data Type Duplicate' : Unused code path elimination
//  Block '<S89>/Data Type Propagation' : Unused code path elimination


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
//  '<S53>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem'
//  '<S54>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block'
//  '<S55>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1'
//  '<S56>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug'
//  '<S57>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter'
//  '<S58>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter'
//  '<S59>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic'
//  '<S60>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1'
//  '<S61>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2'
//  '<S62>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic'
//  '<S63>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic1'
//  '<S64>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Saturation Dynamic2'
//  '<S65>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S66>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S67>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S68>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S69>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S70>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S71>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S72>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter'
//  '<S73>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter'
//  '<S74>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic'
//  '<S75>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1'
//  '<S76>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2'
//  '<S77>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic'
//  '<S78>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic1'
//  '<S79>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Saturation Dynamic2'
//  '<S80>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Natural Frequency'
//  '<S81>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Deriv Filter/Compute Numerator And Denominator'
//  '<S82>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Filter Numerator And Denominator'
//  '<S83>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Discrete Second Order Filter/Compute Natural Frequency'
//  '<S84>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S85>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic1/Saturation Dynamic'
//  '<S86>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/Signal Conditioning Block1/Rate Limiter Dynamic2/Saturation Dynamic'
//  '<S87>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter'
//  '<S88>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic'
//  '<S89>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Saturation Dynamic'
//  '<S90>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Discrete First Order Deriv Filter/Compute Deriv Filter Numerator And Denominator'
//  '<S91>'  : 'fcsModel/Inner Loop Controller/Attitude Controller/For Each Subsystem/pidWithDebug/Rate Limiter Dynamic/Saturation Dynamic'
//  '<S92>'  : 'fcsModel/Outer Loop Controller/assembleOuterLoopToInnerLoopBus'
//  '<S93>'  : 'fcsModel/RC Interpreter/Chart'
//  '<S94>'  : 'fcsModel/RC Interpreter/Interpret RC In Cmds'


//-
//  Requirements for '<Root>': fcsModel

#endif                                 // RTW_HEADER_fcsModel_h_

//
// File trailer for generated code.
//
// [EOF]
//
