//#############AUTOGENERATED##########
#include<fcs_out_data.h>

void AssignFcsData(const FcsOutput &fcs_output, DataFields data_to_save[], size_t data_index){

	data_to_save[data_index].actuator_cmds[0] = fcs_output.actuatorsCmds[0];
	data_to_save[data_index].actuator_cmds[1] = fcs_output.actuatorsCmds[1];
	data_to_save[data_index].actuator_cmds[2] = fcs_output.actuatorsCmds[2];
	data_to_save[data_index].actuator_cmds[3] = fcs_output.actuatorsCmds[3];

	data_to_save[data_index].fcsOut[0] = fcs_output.fcsDebug.allocDebug.thrustCmd_N;
	data_to_save[data_index].fcsOut[1] = fcs_output.fcsDebug.allocDebug.xMomCmd_Nm;
	data_to_save[data_index].fcsOut[2] = fcs_output.fcsDebug.allocDebug.yMomCmd_Nm;
	data_to_save[data_index].fcsOut[3] = fcs_output.fcsDebug.allocDebug.zMomCmd_Nm;
	data_to_save[data_index].fcsOut[4] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[0];
	data_to_save[data_index].fcsOut[5] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[1];
	data_to_save[data_index].fcsOut[6] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.cmd[2];
	data_to_save[data_index].fcsOut[7] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[0];
	data_to_save[data_index].fcsOut[8] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[1];
	data_to_save[data_index].fcsOut[9] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.meas[2];
	data_to_save[data_index].fcsOut[10] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[0].output;
	data_to_save[data_index].fcsOut[11] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[0].proportionalOutput;
	data_to_save[data_index].fcsOut[12] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[0].integralOutput;
	data_to_save[data_index].fcsOut[13] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[0].derivativeOutput;
	data_to_save[data_index].fcsOut[14] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[1].output;
	data_to_save[data_index].fcsOut[15] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[1].proportionalOutput;
	data_to_save[data_index].fcsOut[16] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[1].integralOutput;
	data_to_save[data_index].fcsOut[17] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[1].derivativeOutput;
	data_to_save[data_index].fcsOut[18] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[2].output;
	data_to_save[data_index].fcsOut[19] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[2].proportionalOutput;
	data_to_save[data_index].fcsOut[20] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[2].integralOutput;
	data_to_save[data_index].fcsOut[21] = fcs_output.fcsDebug.innerLoopCtrlDebug.angRateCtrlDebug.pidDebug[2].derivativeOutput;
	data_to_save[data_index].fcsOut[22] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[0];
	data_to_save[data_index].fcsOut[23] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[1];
	data_to_save[data_index].fcsOut[24] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.cmd[2];
	data_to_save[data_index].fcsOut[25] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[0];
	data_to_save[data_index].fcsOut[26] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[1];
	data_to_save[data_index].fcsOut[27] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.meas[2];
	data_to_save[data_index].fcsOut[28] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[0].output;
	data_to_save[data_index].fcsOut[29] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[0].proportionalOutput;
	data_to_save[data_index].fcsOut[30] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[0].integralOutput;
	data_to_save[data_index].fcsOut[31] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[0].derivativeOutput;
	data_to_save[data_index].fcsOut[32] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[1].output;
	data_to_save[data_index].fcsOut[33] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[1].proportionalOutput;
	data_to_save[data_index].fcsOut[34] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[1].integralOutput;
	data_to_save[data_index].fcsOut[35] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[1].derivativeOutput;
	data_to_save[data_index].fcsOut[36] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[2].output;
	data_to_save[data_index].fcsOut[37] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[2].proportionalOutput;
	data_to_save[data_index].fcsOut[38] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[2].integralOutput;
	data_to_save[data_index].fcsOut[39] = fcs_output.fcsDebug.innerLoopCtrlDebug.attCtrlDebug.pidDebug[2].derivativeOutput;
	data_to_save[data_index].fcsOut[40] = fcs_output.fcsDebug.outerLoopCtrlDebug.frcCmd_N;
	data_to_save[data_index].fcsOut[41] = static_cast<float>(fcs_output.fcsDebug.state);
	data_to_save[data_index].fcsOut[42] = static_cast<float>(fcs_output.fcsDebug.flightMode);
}