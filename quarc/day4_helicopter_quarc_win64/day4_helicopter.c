/*
 * day4_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day4_helicopter".
 *
 * Model version              : 11.18
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Tue Apr  5 19:25:33 2022
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "day4_helicopter.h"
#include "day4_helicopter_private.h"
#include "day4_helicopter_dt.h"

/* Block signals (default storage) */
B_day4_helicopter_T day4_helicopter_B;

/* Continuous states */
X_day4_helicopter_T day4_helicopter_X;

/* Block states (default storage) */
DW_day4_helicopter_T day4_helicopter_DW;

/* Real-time model */
static RT_MODEL_day4_helicopter_T day4_helicopter_M_;
RT_MODEL_day4_helicopter_T *const day4_helicopter_M = &day4_helicopter_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  day4_helicopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void day4_helicopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum3_n[2];
  real_T rtb_Sum4[6];
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[8];
  real_T rtb_Gain1[6];
  real_T lastTime;
  real_T rtb_Frontgain;
  real_T *lastU;
  int32_T i;
  int32_T i_0;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* set solver stop time */
    if (!(day4_helicopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&day4_helicopter_M->solverInfo,
                            ((day4_helicopter_M->Timing.clockTickH0 + 1) *
        day4_helicopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&day4_helicopter_M->solverInfo,
                            ((day4_helicopter_M->Timing.clockTick0 + 1) *
        day4_helicopter_M->Timing.stepSize0 +
        day4_helicopter_M->Timing.clockTickH0 *
        day4_helicopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(day4_helicopter_M)) {
    day4_helicopter_M->Timing.t[0] = rtsiGetT(&day4_helicopter_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(day4_helicopter_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: day4_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (day4_helicopter_DW.HILReadEncoderTimebase_Task, 1,
         &day4_helicopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          day4_helicopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          day4_helicopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          day4_helicopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    day4_helicopter_B.TravelCounttorad = day4_helicopter_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * day4_helicopter_P.TravelCounttorad_Gain;

    /* Gain: '<S13>/Gain' */
    day4_helicopter_B.Gain = day4_helicopter_P.Gain_Gain *
      day4_helicopter_B.TravelCounttorad;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/lambda_0'
     */
    day4_helicopter_B.travel = day4_helicopter_P.lambda_0_Value +
      day4_helicopter_B.Gain;

    /* Gain: '<S4>/Pitch: Count to rad' */
    day4_helicopter_B.PitchCounttorad = day4_helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    day4_helicopter_B.Gain_i = day4_helicopter_P.Gain_Gain_a *
      day4_helicopter_B.PitchCounttorad;
  }

  /* Gain: '<S14>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  day4_helicopter_B.Gain_d = (day4_helicopter_P.TravelTransferFcn_C *
    day4_helicopter_X.TravelTransferFcn_CSTATE +
    day4_helicopter_P.TravelTransferFcn_D * day4_helicopter_B.TravelCounttorad) *
    day4_helicopter_P.Gain_Gain_l;

  /* Gain: '<S11>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  day4_helicopter_B.Gain_b = (day4_helicopter_P.PitchTransferFcn_C *
    day4_helicopter_X.PitchTransferFcn_CSTATE +
    day4_helicopter_P.PitchTransferFcn_D * day4_helicopter_B.PitchCounttorad) *
    day4_helicopter_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    day4_helicopter_B.ElevationCounttorad = day4_helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * day4_helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S8>/Gain' */
    day4_helicopter_B.Gain_e = day4_helicopter_P.Gain_Gain_lv *
      day4_helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    day4_helicopter_B.elevation = day4_helicopter_B.Gain_e +
      day4_helicopter_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  day4_helicopter_B.Gain_dg = (day4_helicopter_P.ElevationTransferFcn_C *
    day4_helicopter_X.ElevationTransferFcn_CSTATE +
    day4_helicopter_P.ElevationTransferFcn_D *
    day4_helicopter_B.ElevationCounttorad) * day4_helicopter_P.Gain_Gain_n;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      day4_helicopter_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      day4_helicopter_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = day4_helicopter_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = day4_helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[114]) {
      currTimeIndex = 113;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    day4_helicopter_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum3_n[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 115;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum3_n[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 115;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_n[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 115;
          }
        }
      }
    }
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1[0] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.travel;
  rtb_Gain1[1] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_d;
  rtb_Gain1[2] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_i;
  rtb_Gain1[3] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_b;
  rtb_Gain1[4] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.elevation;
  rtb_Gain1[5] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_dg;

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      day4_helicopter_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      day4_helicopter_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = day4_helicopter_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = day4_helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[114]) {
      currTimeIndex = 113;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    day4_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 115;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 115;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum4[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 115;
          }
        }
      }
    }
  }

  for (i = 0; i < 6; i++) {
    /* Sum: '<S5>/Sum4' */
    rtb_Sum4[i] = rtb_Gain1[i] - rtb_Sum4[i];
  }

  for (i = 0; i < 2; i++) {
    /* Sum: '<S5>/Sum3' incorporates:
     *  Gain: '<S5>/K'
     *  Sum: '<S5>/Sum4'
     */
    rtb_Frontgain = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Frontgain += day4_helicopter_P.K[(i_0 << 1) + i] * rtb_Sum4[i_0];
    }

    rtb_Sum3_n[i] -= rtb_Frontgain;

    /* End of Sum: '<S5>/Sum3' */
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  day4_helicopter_B.Sum1 = ((rtb_Sum3_n[0] - rtb_Gain1[2]) *
    day4_helicopter_P.K_pp - day4_helicopter_P.K_pd * rtb_Gain1[3]) +
    day4_helicopter_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (day4_helicopter_X.Integrator_CSTATE >=
      day4_helicopter_P.Integrator_UpperSat) {
    day4_helicopter_X.Integrator_CSTATE = day4_helicopter_P.Integrator_UpperSat;
  } else {
    if (day4_helicopter_X.Integrator_CSTATE <=
        day4_helicopter_P.Integrator_LowerSat) {
      day4_helicopter_X.Integrator_CSTATE =
        day4_helicopter_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' */
  rtb_Frontgain = rtb_Sum3_n[1] - rtb_Gain1[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  day4_helicopter_B.Sum2 = ((day4_helicopter_P.K_ep * rtb_Frontgain +
    day4_helicopter_X.Integrator_CSTATE) - day4_helicopter_P.K_ed * rtb_Gain1[5])
    + day4_helicopter_P.Vs_ff;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* SignalConversion generated from: '<Root>/To File3' */
    rtb_TmpSignalConversionAtToFile[0] = day4_helicopter_B.travel;
    rtb_TmpSignalConversionAtToFile[1] = day4_helicopter_B.Gain_d;
    rtb_TmpSignalConversionAtToFile[2] = day4_helicopter_B.Gain_i;
    rtb_TmpSignalConversionAtToFile[3] = day4_helicopter_B.Gain_b;
    rtb_TmpSignalConversionAtToFile[4] = day4_helicopter_B.elevation;
    rtb_TmpSignalConversionAtToFile[5] = day4_helicopter_B.Gain_dg;
    rtb_TmpSignalConversionAtToFile[6] = day4_helicopter_B.Sum1;
    rtb_TmpSignalConversionAtToFile[7] = day4_helicopter_B.Sum2;

    /* ToFile: '<Root>/To File3' */
    {
      if (!(++day4_helicopter_DW.ToFile3_IWORK.Decimation % 1) &&
          (day4_helicopter_DW.ToFile3_IWORK.Count * (8 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) day4_helicopter_DW.ToFile3_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[8 + 1];
          day4_helicopter_DW.ToFile3_IWORK.Decimation = 0;
          u[0] = day4_helicopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          u[7] = rtb_TmpSignalConversionAtToFile[6];
          u[8] = rtb_TmpSignalConversionAtToFile[7];
          if (fwrite(u, sizeof(real_T), 8 + 1, fp) != 8 + 1) {
            rtmSetErrorStatus(day4_helicopter_M,
                              "Error writing to MAT-file day_4_test6.mat");
            return;
          }

          if (((++day4_helicopter_DW.ToFile3_IWORK.Count) * (8 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file day_4_test6.mat.\n");
          }
        }
      }
    }
  }

  /* If: '<S3>/If' incorporates:
   *  Clock: '<S3>/Clock'
   *  Gain: '<S3>/K_ei'
   *  Inport: '<S7>/In1'
   */
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    rtAction = (int8_T)!(day4_helicopter_M->Timing.t[0] >= 2.0);
    day4_helicopter_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = day4_helicopter_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S7>/Action Port'
     */
    day4_helicopter_B.In1 = day4_helicopter_P.K_ei * rtb_Frontgain;
    if (rtmIsMajorTimeStep(day4_helicopter_M)) {
      srUpdateBC(day4_helicopter_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  rtb_Frontgain = day4_helicopter_M->Timing.t[0];
  if ((day4_helicopter_DW.TimeStampA >= rtb_Frontgain) &&
      (day4_helicopter_DW.TimeStampB >= rtb_Frontgain)) {
    rtb_Frontgain = 0.0;
  } else {
    lastTime = day4_helicopter_DW.TimeStampA;
    lastU = &day4_helicopter_DW.LastUAtTimeA;
    if (day4_helicopter_DW.TimeStampA < day4_helicopter_DW.TimeStampB) {
      if (day4_helicopter_DW.TimeStampB < rtb_Frontgain) {
        lastTime = day4_helicopter_DW.TimeStampB;
        lastU = &day4_helicopter_DW.LastUAtTimeB;
      }
    } else {
      if (day4_helicopter_DW.TimeStampA >= rtb_Frontgain) {
        lastTime = day4_helicopter_DW.TimeStampB;
        lastU = &day4_helicopter_DW.LastUAtTimeB;
      }
    }

    rtb_Frontgain = (day4_helicopter_B.PitchCounttorad - *lastU) /
      (rtb_Frontgain - lastTime);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S12>/Gain' */
  day4_helicopter_B.Gain_l = day4_helicopter_P.Gain_Gain_a1 * rtb_Frontgain;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Frontgain = (day4_helicopter_B.Sum2 - day4_helicopter_B.Sum1) *
    day4_helicopter_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Frontgain > day4_helicopter_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day4_helicopter_B.BackmotorSaturation =
      day4_helicopter_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < day4_helicopter_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day4_helicopter_B.BackmotorSaturation =
      day4_helicopter_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day4_helicopter_B.BackmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Frontgain = (day4_helicopter_B.Sum1 + day4_helicopter_B.Sum2) *
    day4_helicopter_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Frontgain > day4_helicopter_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day4_helicopter_B.FrontmotorSaturation =
      day4_helicopter_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < day4_helicopter_P.FrontmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day4_helicopter_B.FrontmotorSaturation =
      day4_helicopter_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day4_helicopter_B.FrontmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: day4_helicopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      day4_helicopter_DW.HILWriteAnalog_Buffer[0] =
        day4_helicopter_B.FrontmotorSaturation;
      day4_helicopter_DW.HILWriteAnalog_Buffer[1] =
        day4_helicopter_B.BackmotorSaturation;
      result = hil_write_analog(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILWriteAnalog_channels, 2,
        &day4_helicopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void day4_helicopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (day4_helicopter_DW.TimeStampA == (rtInf)) {
    day4_helicopter_DW.TimeStampA = day4_helicopter_M->Timing.t[0];
    lastU = &day4_helicopter_DW.LastUAtTimeA;
  } else if (day4_helicopter_DW.TimeStampB == (rtInf)) {
    day4_helicopter_DW.TimeStampB = day4_helicopter_M->Timing.t[0];
    lastU = &day4_helicopter_DW.LastUAtTimeB;
  } else if (day4_helicopter_DW.TimeStampA < day4_helicopter_DW.TimeStampB) {
    day4_helicopter_DW.TimeStampA = day4_helicopter_M->Timing.t[0];
    lastU = &day4_helicopter_DW.LastUAtTimeA;
  } else {
    day4_helicopter_DW.TimeStampB = day4_helicopter_M->Timing.t[0];
    lastU = &day4_helicopter_DW.LastUAtTimeB;
  }

  *lastU = day4_helicopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    rt_ertODEUpdateContinuousStates(&day4_helicopter_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++day4_helicopter_M->Timing.clockTick0)) {
    ++day4_helicopter_M->Timing.clockTickH0;
  }

  day4_helicopter_M->Timing.t[0] = rtsiGetSolverStopTime
    (&day4_helicopter_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++day4_helicopter_M->Timing.clockTick1)) {
      ++day4_helicopter_M->Timing.clockTickH1;
    }

    day4_helicopter_M->Timing.t[1] = day4_helicopter_M->Timing.clockTick1 *
      day4_helicopter_M->Timing.stepSize1 +
      day4_helicopter_M->Timing.clockTickH1 *
      day4_helicopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void day4_helicopter_derivatives(void)
{
  XDot_day4_helicopter_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_day4_helicopter_T *) day4_helicopter_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += day4_helicopter_P.TravelTransferFcn_A *
    day4_helicopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += day4_helicopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += day4_helicopter_P.PitchTransferFcn_A *
    day4_helicopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += day4_helicopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    day4_helicopter_P.ElevationTransferFcn_A *
    day4_helicopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += day4_helicopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (day4_helicopter_X.Integrator_CSTATE <=
          day4_helicopter_P.Integrator_LowerSat);
  usat = (day4_helicopter_X.Integrator_CSTATE >=
          day4_helicopter_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (day4_helicopter_B.In1 > 0.0)) || (usat &&
       (day4_helicopter_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = day4_helicopter_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void day4_helicopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day4_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &day4_helicopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(day4_helicopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(day4_helicopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
      return;
    }

    if ((day4_helicopter_P.HILInitialize_AIPStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &day4_helicopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (day4_helicopter_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &day4_helicopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = day4_helicopter_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_AIChannels, 8U,
        &day4_helicopter_DW.HILInitialize_AIMinimums[0],
        &day4_helicopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_AOPStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &day4_helicopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (day4_helicopter_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &day4_helicopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = day4_helicopter_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (day4_helicopter_DW.HILInitialize_Card,
         day4_helicopter_P.HILInitialize_AOChannels, 8U,
         &day4_helicopter_DW.HILInitialize_AOMinimums[0],
         &day4_helicopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_AOStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day4_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day4_helicopter_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_AOChannels, 8U,
        &day4_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day4_helicopter_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day4_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day4_helicopter_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (day4_helicopter_DW.HILInitialize_Card,
         day4_helicopter_P.HILInitialize_AOChannels, 8U,
         &day4_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_EIPStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &day4_helicopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = day4_helicopter_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (day4_helicopter_DW.HILInitialize_Card,
         day4_helicopter_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &day4_helicopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_EIStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &day4_helicopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = day4_helicopter_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_EIChannels, 8U,
        &day4_helicopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_POPStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day4_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day4_helicopter_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &day4_helicopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          day4_helicopter_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &day4_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            day4_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            day4_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = day4_helicopter_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            day4_helicopter_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            day4_helicopter_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] = day4_helicopter_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(day4_helicopter_DW.HILInitialize_Card,
          &day4_helicopter_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes, &day4_helicopter_DW.HILInitialize_POSortedFreqs
          [0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(day4_helicopter_DW.HILInitialize_Card,
          &day4_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &day4_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day4_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day4_helicopter_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &day4_helicopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = day4_helicopter_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &day4_helicopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = day4_helicopter_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &day4_helicopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &day4_helicopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &day4_helicopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &day4_helicopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = day4_helicopter_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &day4_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day4_helicopter_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_POChannels, 8U,
        &day4_helicopter_DW.HILInitialize_POSortedFreqs[0],
        &day4_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day4_helicopter_P.HILInitialize_POStart && !is_switching) ||
        (day4_helicopter_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day4_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day4_helicopter_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(day4_helicopter_DW.HILInitialize_Card,
        day4_helicopter_P.HILInitialize_POChannels, 8U,
        &day4_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day4_helicopter_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &day4_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day4_helicopter_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (day4_helicopter_DW.HILInitialize_Card,
         day4_helicopter_P.HILInitialize_POChannels, 8U,
         &day4_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: day4_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (day4_helicopter_DW.HILInitialize_Card,
       day4_helicopter_P.HILReadEncoderTimebase_SamplesI,
       day4_helicopter_P.HILReadEncoderTimebase_Channels, 3,
       &day4_helicopter_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (day4_helicopter_DW.HILReadEncoderTimebase_Task, (t_buffer_overflow_mode)
         (day4_helicopter_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.46784823682248955, 0.29913179319090455,
      0.15048752350343572, 0.023130014521000442, -0.080696640985430959,
      -0.15750989621152992, -0.20277442175028984, -0.23778294351613158,
      -0.26364920992175994, -0.28144129092901066, -0.29217576846920218,
      -0.29680784414061956, -0.29623106297493196, -0.29126896442573424,
      -0.28268182791967195, -0.2711586526834463, -0.25732475660783583,
      -0.24173910635673246, -0.22490098882766779, -0.20724874349686168,
      -0.18916889361050906, -0.170994810623596, -0.15301179617721408,
      -0.13546649404895825, -0.11856384507486635, -0.1024735746095255,
      -0.087334069601914047, -0.073261006322447947, -0.060341217162964644,
      -0.048642678644066954, -0.038214570810246325, -0.029090169213080271,
      -0.021284628505903111, -0.014797346886869721, -0.0096067347972797878,
      -0.0056657570942712007, -0.0028933246688912159, -0.0011619416216426819,
      -0.00028334185863388084, -1.7958405195219585E-7, 3.1344470552099588E-7,
      -3.7243833812917228E-8, -3.7243833812917228E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.2498902065152821, 0.26617674763616739, 0.28201655248111829,
      0.2969498692189223, 0.31039426491307009, 0.32162006169698565,
      0.32971832445114785, 0.3335666451154406, 0.33178575376519837,
      0.32268848877859141, 0.30422072258468086, 0.27389102126336473,
      0.22868701075141265, 0.16497871916144324, 0.078400122850391857,
      -1.6006974542191707E-7, -2.7365997634152278E-7, -7.9763056441242943E-8,
      -2.3030072800441832E-7, -1.2652274110234629E-7, -8.141284561074491E-8,
      -8.0996578586120672E-8, 1.2666496763792769E-7, 2.2004930271401703E-7,
      1.171675829487724E-7, -1.0523643129006975E-7, -1.145498222603564E-7,
      2.6818629061872708E-7, -2.5277749598500126E-7, -1.8994725301261049E-7,
      -3.5698753082089257E-7, -1.1270959650701076E-7, 8.0500523159642024E-8,
      -1.6983277381375507E-7, 2.3229458487777172E-7, -6.5757667149162828E-8,
      6.1028134603606655E-8, -1.7960416578181952E-7, -5.6522012537521482E-8,
      1.0817338845172151E-7, -1.3137139046142607E-7, -1.1500613380366367E-7,
      1.317761043544438E-7, 4.2904834394034227E-8, -1.4361817688976217E-7,
      1.1482561762097056E-7, -6.2971666192638952E-9, 1.7559962904401394E-7,
      1.2961812339089434E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    day4_helicopter_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    day4_helicopter_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    day4_helicopter_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625261, 3.1262155534579983,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595077632935446, 2.7555515879651526,
      2.633904448852213, 2.4960414008875862, 2.344614435015925,
      2.183006792719969, 2.0148815175869066, 1.843794736576035,
      1.6728743462233042, 1.5047360686293036, 1.3414962881941592,
      1.1848192156220627, 1.0359721669402255, 0.89587894537033252,
      0.76516823990293914, 0.64421661468571545, 0.53318660262582651,
      0.43206053933011523, 0.34067073624582073, 0.25872644019912994,
      0.18583793787996089, 0.12153805531887983, 0.06530128435202745,
      0.016560706064074868, -0.0252771355785695, -0.060819225375501849,
      -0.090674798997774933, -0.11544565659610218, -0.13571792807408711,
      -0.15205511643886929, -0.1649923145155392, -0.17503146080695711,
      -0.18263751810920792, -0.18823546673987032, -0.1922080445059027,
      -0.19489417452744306, -0.19658807414864868, -0.19753908723333063,
      -0.19795235008191422, -0.19799048791751137, -0.19777662812147254, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909068423, -0.046506351618112111,
      -0.091625013712135384, -0.14672753935514368, -0.2086939711880792,
      -0.27519058110636668, -0.34459405529608833, -0.41582470131356863,
      -0.48658855645175897, -0.55145219185850791, -0.60570786348664352,
      -0.64643056918382513, -0.67250110053224865, -0.68434712404348819,
      -0.68368156141092251, -0.67255311037600218, -0.65295912174057824,
      -0.62670829028838571, -0.59538819472734894, -0.560372886279572,
      -0.5228428218695732, -0.48380650086889493, -0.4441200482395562,
      -0.40450425318284483, -0.36555921233717786, -0.32777718418676333,
      -0.29155400927667607, -0.25719953024432424, -0.2249470838674095,
      -0.19496231315181034, -0.16735136657057748, -0.14216835918772938,
      -0.11942229448909231, -0.099083430393308983, -0.081089085911939751,
      -0.065348753459128689, -0.051748792306679683, -0.040156585165671638,
      -0.030424229209003246, -0.02239179452264967, -0.015890311064129522,
      -0.010744520086161448, -0.0067755984848224763, -0.0038040523387277789,
      -0.0016530513943343855, -0.00015255134238856177, 0.0008554391841553811,
      0.0015103326430469643, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10602875205865551,
      0.22266037932317656, 0.31888147181640641, 0.38944360631144165,
      0.43795507377677839, 0.46997264230390062, 0.49051724877547076,
      0.5034310014147434, 0.500131901758192, 0.4584314021272094,
      0.38345836556752305, 0.28781253092129711, 0.18425655862970131,
      0.083723169906144257, -0.0047039424931596642, -0.078651641701100908,
      -0.13848282827618658, -0.18553085090234689, -0.22135847355404858,
      -0.24747482694996897, -0.26524816164728987, -0.27589380795586449,
      -0.28048869000687932, -0.27998930926931609, -0.2752486746822535,
      -0.26702894513374908, -0.256011565793345, -0.24280433702783488,
      -0.22794797303425959, -0.21192090753804371, -0.19514362517487038,
      -0.17798387820701325, -0.16076049803997672, -0.14374732354478181,
      -0.12717715433655463, -0.11124665817820928, -0.096119331284722789,
      -0.081929292739710322, -0.068784574889127578, -0.056770180589926736,
      -0.045950002017403686, -0.036368485335282343, -0.028050822054995112,
      -0.021001753258161632, -0.015202453157792247, -0.01060496129136022,
      -0.0071240920666374484, -0.0046285368484352473, -0.0029362892105680755,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462206, 0.46652650905808424,
      0.38488436997291947, 0.28224853798014093, 0.19404586986134692,
      0.12807027410848898, 0.082178425886280437, 0.0516550105570906,
      -0.013196398626205703, -0.16680199852393002, -0.2998921462387456,
      -0.38258333858490384, -0.41422388916638314, -0.4021335548942282,
      -0.35370844959721565, -0.295790796831765, -0.23932474630034273,
      -0.18819209050464122, -0.14331049060680678, -0.10446541358368151,
      -0.0710933387892837, -0.042582585234298261, -0.0183795282040593,
      0.0019975229502527747, 0.01896253834825044, 0.032878918194017774,
      0.0440695173616162, 0.052828915062040412, 0.05942545597430117,
      0.064108261984863607, 0.06710912945269322, 0.068638987871428539,
      0.068893520668146149, 0.068052697980779717, 0.066280676832908766,
      0.063721984633381315, 0.060509307573946014, 0.056760154180049892,
      0.052578871402330925, 0.048057577196803439, 0.04328071429009217,
      0.038326066728485396, 0.0332706531211489, 0.028196275187333912,
      0.023197200401477543, 0.018389967465728106, 0.013923476898891087,
      0.0099822208728088061, 0.0067689905514686857, 0.00442598373487403, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0039045344768012828, 0.010991947016217361,
      0.020653006702096872, 0.032366894000923768, 0.0456795166345904,
      0.06018356435509134, 0.075499706517600754, 0.091258423776363284,
      0.10708193120867888, 0.1225656565485754, 0.13725874416874784,
      0.15064301870754573, 0.16210979127699829, 0.17093386602367985,
      0.17624395851452532, 0.17755568372494962, 0.17578566150504127,
      0.17168383603560991, 0.16586081237407124, 0.15881098271294303,
      0.15093203400167629, 0.14254139959776493, 0.13389011274269536,
      0.12517444167094846, 0.11654565718627714, 0.10811821652734585,
      0.099976608349770832, 0.092181059273760566, 0.084772259011639345,
      0.077775326795970073, 0.071203055509231089, 0.065058610801902381,
      0.059337730786894784, 0.054030527328222271, 0.049122976320275541,
      0.044598085047348482, 0.040436871041213469, 0.036619112651432306,
      0.033123966865921034, 0.029930435581818141, 0.027017723083782981,
      0.024365523857319858, 0.021954224573290035, 0.019765039470385083,
      0.0177801136402147, 0.015982592320012423, 0.014356636955839131,
      0.012887445171453172, 0.011561230906011893, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.015618137907205131, 0.028349650157664313, 0.03864423874351805,
      0.046855549195307593, 0.053250490534666524, 0.058016190882003757,
      0.061264568650037661, 0.063034869035050067, 0.0632940297292624,
      0.061934901359586063, 0.058772350480689677, 0.053537098155191587,
      0.045867090277810239, 0.035296298986726274, 0.021240369963381807,
      0.0052469008416972749, -0.0070800888796333974, -0.016407301877725426,
      -0.02329209464615465, -0.028199318644512927, -0.031515794845066995,
      -0.033562537615645351, -0.0346051474202783, -0.034862684286987612,
      -0.034515137938685235, -0.033709762635725157, -0.032566432710300082,
      -0.031182196304041013, -0.029635201048484935, -0.027987728862677051,
      -0.026289085146955923, -0.024577778829314851, -0.022883520060030384,
      -0.021228813834690052, -0.019630204031786908, -0.018099565091708267,
      -0.016644856024540012, -0.015271033559124649, -0.013980583142045111,
      -0.012774125136411576, -0.01165084999214065, -0.010608796905852485,
      -0.0096451971361192768, -0.008756740411619798, -0.007939703320681531,
      -0.0071900852808091152, -0.00650382145669317, -0.0058767671375438381,
      -0.0053048570617651118, -0.0047841081195396573, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    day4_helicopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    day4_helicopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    day4_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "day_4_test6.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(day4_helicopter_M,
                        "Error creating .mat file day_4_test6.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 8 + 1, 0, "state")) {
      rtmSetErrorStatus(day4_helicopter_M,
                        "Error writing mat file header to file day_4_test6.mat");
      return;
    }

    day4_helicopter_DW.ToFile3_IWORK.Count = 0;
    day4_helicopter_DW.ToFile3_IWORK.Decimation = -1;
    day4_helicopter_DW.ToFile3_PWORK.FilePtr = fp;
  }

  /* Start for If: '<S3>/If' */
  day4_helicopter_DW.If_ActiveSubsystem = -1;

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  day4_helicopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  day4_helicopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  day4_helicopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  day4_helicopter_X.Integrator_CSTATE = day4_helicopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  day4_helicopter_DW.TimeStampA = (rtInf);
  day4_helicopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void day4_helicopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day4_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(day4_helicopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(day4_helicopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((day4_helicopter_P.HILInitialize_AOTerminate && !is_switching) ||
        (day4_helicopter_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day4_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day4_helicopter_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((day4_helicopter_P.HILInitialize_POTerminate && !is_switching) ||
        (day4_helicopter_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day4_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day4_helicopter_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    } else {
      num_final_pwm_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(day4_helicopter_DW.HILInitialize_Card
                         , day4_helicopter_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , day4_helicopter_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &day4_helicopter_DW.HILInitialize_AOVoltages[0]
                         , &day4_helicopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(day4_helicopter_DW.HILInitialize_Card,
            day4_helicopter_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &day4_helicopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(day4_helicopter_DW.HILInitialize_Card,
            day4_helicopter_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &day4_helicopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day4_helicopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(day4_helicopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(day4_helicopter_DW.HILInitialize_Card);
    hil_close(day4_helicopter_DW.HILInitialize_Card);
    day4_helicopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (FILE *) day4_helicopter_DW.ToFile3_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "day_4_test6.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error closing MAT-file day_4_test6.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error reopening MAT-file day_4_test6.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 8 + 1,
           day4_helicopter_DW.ToFile3_IWORK.Count, "state")) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error writing header for state to MAT-file day_4_test6.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error closing MAT-file day_4_test6.mat");
        return;
      }

      day4_helicopter_DW.ToFile3_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  day4_helicopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  day4_helicopter_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  day4_helicopter_initialize();
}

void MdlTerminate(void)
{
  day4_helicopter_terminate();
}

/* Registration function */
RT_MODEL_day4_helicopter_T *day4_helicopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  day4_helicopter_P.Integrator_UpperSat = rtInf;
  day4_helicopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)day4_helicopter_M, 0,
                sizeof(RT_MODEL_day4_helicopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&day4_helicopter_M->solverInfo,
                          &day4_helicopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&day4_helicopter_M->solverInfo, &rtmGetTPtr(day4_helicopter_M));
    rtsiSetStepSizePtr(&day4_helicopter_M->solverInfo,
                       &day4_helicopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&day4_helicopter_M->solverInfo, &day4_helicopter_M->derivs);
    rtsiSetContStatesPtr(&day4_helicopter_M->solverInfo, (real_T **)
                         &day4_helicopter_M->contStates);
    rtsiSetNumContStatesPtr(&day4_helicopter_M->solverInfo,
      &day4_helicopter_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&day4_helicopter_M->solverInfo,
      &day4_helicopter_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&day4_helicopter_M->solverInfo,
      &day4_helicopter_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&day4_helicopter_M->solverInfo,
      &day4_helicopter_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&day4_helicopter_M->solverInfo, (&rtmGetErrorStatus
      (day4_helicopter_M)));
    rtsiSetRTModelPtr(&day4_helicopter_M->solverInfo, day4_helicopter_M);
  }

  rtsiSetSimTimeStep(&day4_helicopter_M->solverInfo, MAJOR_TIME_STEP);
  day4_helicopter_M->intgData.f[0] = day4_helicopter_M->odeF[0];
  day4_helicopter_M->contStates = ((real_T *) &day4_helicopter_X);
  rtsiSetSolverData(&day4_helicopter_M->solverInfo, (void *)
                    &day4_helicopter_M->intgData);
  rtsiSetSolverName(&day4_helicopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = day4_helicopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    day4_helicopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    day4_helicopter_M->Timing.sampleTimes =
      (&day4_helicopter_M->Timing.sampleTimesArray[0]);
    day4_helicopter_M->Timing.offsetTimes =
      (&day4_helicopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    day4_helicopter_M->Timing.sampleTimes[0] = (0.0);
    day4_helicopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    day4_helicopter_M->Timing.offsetTimes[0] = (0.0);
    day4_helicopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(day4_helicopter_M, &day4_helicopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = day4_helicopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    day4_helicopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(day4_helicopter_M, 25.0);
  day4_helicopter_M->Timing.stepSize0 = 0.002;
  day4_helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  day4_helicopter_M->Sizes.checksums[0] = (1114028809U);
  day4_helicopter_M->Sizes.checksums[1] = (78721603U);
  day4_helicopter_M->Sizes.checksums[2] = (806027582U);
  day4_helicopter_M->Sizes.checksums[3] = (3861426760U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    day4_helicopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &day4_helicopter_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(day4_helicopter_M->extModeInfo,
      &day4_helicopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(day4_helicopter_M->extModeInfo,
                        day4_helicopter_M->Sizes.checksums);
    rteiSetTPtr(day4_helicopter_M->extModeInfo, rtmGetTPtr(day4_helicopter_M));
  }

  day4_helicopter_M->solverInfoPtr = (&day4_helicopter_M->solverInfo);
  day4_helicopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&day4_helicopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&day4_helicopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  day4_helicopter_M->blockIO = ((void *) &day4_helicopter_B);

  {
    day4_helicopter_B.TravelCounttorad = 0.0;
    day4_helicopter_B.Gain = 0.0;
    day4_helicopter_B.travel = 0.0;
    day4_helicopter_B.Gain_d = 0.0;
    day4_helicopter_B.PitchCounttorad = 0.0;
    day4_helicopter_B.Gain_i = 0.0;
    day4_helicopter_B.Gain_b = 0.0;
    day4_helicopter_B.ElevationCounttorad = 0.0;
    day4_helicopter_B.Gain_e = 0.0;
    day4_helicopter_B.elevation = 0.0;
    day4_helicopter_B.Gain_dg = 0.0;
    day4_helicopter_B.Sum1 = 0.0;
    day4_helicopter_B.Sum2 = 0.0;
    day4_helicopter_B.Gain_l = 0.0;
    day4_helicopter_B.BackmotorSaturation = 0.0;
    day4_helicopter_B.FrontmotorSaturation = 0.0;
    day4_helicopter_B.In1 = 0.0;
  }

  /* parameters */
  day4_helicopter_M->defaultParam = ((real_T *)&day4_helicopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &day4_helicopter_X;
    day4_helicopter_M->contStates = (x);
    (void) memset((void *)&day4_helicopter_X, 0,
                  sizeof(X_day4_helicopter_T));
  }

  /* states (dwork) */
  day4_helicopter_M->dwork = ((void *) &day4_helicopter_DW);
  (void) memset((void *)&day4_helicopter_DW, 0,
                sizeof(DW_day4_helicopter_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day4_helicopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  day4_helicopter_DW.TimeStampA = 0.0;
  day4_helicopter_DW.LastUAtTimeA = 0.0;
  day4_helicopter_DW.TimeStampB = 0.0;
  day4_helicopter_DW.LastUAtTimeB = 0.0;
  day4_helicopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  day4_helicopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    day4_helicopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  day4_helicopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  day4_helicopter_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  day4_helicopter_M->Sizes.numY = (0); /* Number of model outputs */
  day4_helicopter_M->Sizes.numU = (0); /* Number of model inputs */
  day4_helicopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  day4_helicopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  day4_helicopter_M->Sizes.numBlocks = (65);/* Number of blocks */
  day4_helicopter_M->Sizes.numBlockIO = (17);/* Number of block outputs */
  day4_helicopter_M->Sizes.numBlockPrms = (156);/* Sum of parameter "widths" */
  return day4_helicopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
