/*
 * day4_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day4_helicopter".
 *
 * Model version              : 11.13
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Sun Apr  3 21:05:33 2022
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
  real_T rtb_Sum3_c[2];
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
     *  Constant: '<Root>/Constant'
     */
    day4_helicopter_B.travel = day4_helicopter_P.Constant_Value +
      day4_helicopter_B.Gain;

    /* Gain: '<S4>/Pitch: Count to rad' */
    day4_helicopter_B.PitchCounttorad = day4_helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    day4_helicopter_B.Gain_n = day4_helicopter_P.Gain_Gain_j *
      day4_helicopter_B.PitchCounttorad;
  }

  /* Gain: '<S14>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  day4_helicopter_B.Gain_b = (day4_helicopter_P.TravelTransferFcn_C *
    day4_helicopter_X.TravelTransferFcn_CSTATE +
    day4_helicopter_P.TravelTransferFcn_D * day4_helicopter_B.TravelCounttorad) *
    day4_helicopter_P.Gain_Gain_b;

  /* Gain: '<S11>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  day4_helicopter_B.Gain_k = (day4_helicopter_P.PitchTransferFcn_C *
    day4_helicopter_X.PitchTransferFcn_CSTATE +
    day4_helicopter_P.PitchTransferFcn_D * day4_helicopter_B.PitchCounttorad) *
    day4_helicopter_P.Gain_Gain_f;
  if (rtmIsMajorTimeStep(day4_helicopter_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    day4_helicopter_B.ElevationCounttorad = day4_helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * day4_helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S8>/Gain' */
    day4_helicopter_B.Gain_d = day4_helicopter_P.Gain_Gain_e *
      day4_helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    day4_helicopter_B.elevation = day4_helicopter_B.Gain_d +
      day4_helicopter_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  day4_helicopter_B.Gain_i = (day4_helicopter_P.ElevationTransferFcn_C *
    day4_helicopter_X.ElevationTransferFcn_CSTATE +
    day4_helicopter_P.ElevationTransferFcn_D *
    day4_helicopter_B.ElevationCounttorad) * day4_helicopter_P.Gain_Gain_fa;
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
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
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
              (&rtb_Sum3_c[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum3_c[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
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
            (&rtb_Sum3_c[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 105;
          }
        }
      }
    }
  }

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
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
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
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
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
            pDataValues += 105;
          }
        }
      }
    }
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1[0] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.travel;
  rtb_Gain1[1] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_b;
  rtb_Gain1[2] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_n;
  rtb_Gain1[3] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_k;
  rtb_Gain1[4] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.elevation;
  rtb_Gain1[5] = day4_helicopter_P.Gain1_Gain * day4_helicopter_B.Gain_i;
  for (i = 0; i < 6; i++) {
    /* Sum: '<S5>/Sum4' */
    rtb_Sum4[i] = rtb_Gain1[i] - rtb_Sum4[i];
  }

  for (i = 0; i < 2; i++) {
    /* Sum: '<S5>/Sum3' incorporates:
     *  Gain: '<S5>/Gain'
     *  Sum: '<S5>/Sum4'
     */
    rtb_Frontgain = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Frontgain += day4_helicopter_P.K[(i_0 << 1) + i] * rtb_Sum4[i_0];
    }

    rtb_Sum3_c[i] -= rtb_Frontgain;

    /* End of Sum: '<S5>/Sum3' */
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  day4_helicopter_B.Sum1 = ((rtb_Sum3_c[0] - rtb_Gain1[2]) *
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
  rtb_Frontgain = rtb_Sum3_c[1] - rtb_Gain1[4];

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
    rtb_TmpSignalConversionAtToFile[1] = day4_helicopter_B.Gain_b;
    rtb_TmpSignalConversionAtToFile[2] = day4_helicopter_B.Gain_n;
    rtb_TmpSignalConversionAtToFile[3] = day4_helicopter_B.Gain_k;
    rtb_TmpSignalConversionAtToFile[4] = day4_helicopter_B.elevation;
    rtb_TmpSignalConversionAtToFile[5] = day4_helicopter_B.Gain_i;
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
                              "Error writing to MAT-file day_4_test5.mat");
            return;
          }

          if (((++day4_helicopter_DW.ToFile3_IWORK.Count) * (8 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file day_4_test5.mat.\n");
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
  day4_helicopter_B.Gain_l = day4_helicopter_P.Gain_Gain_ez * rtb_Frontgain;
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
      25.25, 25.5, 25.75, 26.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.46880129875785248, 0.29738098184690115,
      0.14750661331313214, 0.020980874469688594, -0.079395321115377449,
      -0.15004573132484342, -0.18794406340120504, -0.21560745678053195,
      -0.23416176364198316, -0.2446945118962994, -0.24824792926182584,
      -0.24581489470274154, -0.23833256523298543, -0.226685579020395,
      -0.21170471841640587, -0.19416893144397304, -0.17480506165909721,
      -0.15429634144926702, -0.13327672930174753, -0.11234127250945362,
      -0.092045605402452713, -0.072898823188007666, -0.055371538131400277,
      -0.039877974957394896, -0.026765116988997448, -0.016287279790507658,
      -0.00856604445280804, -0.0035356631342305322, -0.00088297334106124316,
      1.0149737605931913E-6, 7.4073186988145789E-7, 1.9308955999310042E-7,
      1.9308955999310042E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24886556374727972,
      0.265150739608283, 0.28100945952432771, 0.29598860962343787,
      0.30951221484434571, 0.32085973388040517, 0.32913203163219512,
      0.33322014640530856, 0.33175957943733053, 0.32308093533655391,
      0.30515223684762416, 0.27550728085282294, 0.23116475973130957,
      0.16853067749761674, 0.083283395727221182, 2.2262445038413904E-7,
      9.4262245654439826E-8, 2.9912121991169006E-7, -2.6861190635851315E-7,
      3.9711782837505577E-7, 2.1263594833066429E-7, 3.3450211427413076E-7,
      -1.8500062962050212E-7, 1.1280385014972778E-7, 2.5385897579406571E-8,
      -3.79020313346709E-8, 1.4998841428945265E-7, 9.9022237782825426E-8,
      -3.5683129943904862E-7, -1.5770476846417069E-7, 4.9935879443163853E-7,
      -2.924359409425171E-7, 7.3436922412886322E-8, -4.343196628558939E-7,
      3.034620577061055E-7, -1.2431436378417895E-7, 6.3297476325745837E-8,
      -3.6476157623878086E-7, 1.3553382676659031E-7, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      25.25, 25.5, 25.75, 26.0 } ;

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
      2.6338976221165358, 2.4960327789823658, 2.3446329698642492,
      2.183098203572758, 2.0150780163733386, 1.8440644608995427,
      1.6730662742130578, 1.5045348726836452, 1.3403895749425978,
      1.1820755235891869, 1.0306272345501266, 0.88672809200187408,
      0.75076328191342823, 0.62286635834902, 0.50296040485741422,
      0.39079480525988181, 0.28597844654127996, 0.18801001282679985,
      0.096305805895242125, 0.010225446110089145, -0.07090429228025609,
      -0.14776739021997443, -0.22103812303694281, -0.29136226898572165,
      -0.35934014261508579, -0.4255117942382, -0.49034499087677158,
      -0.5542269434084397, -0.61746108536080824, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, -0.015002048909068423, -0.046506351618112111, -0.091625013712135384,
      -0.14672753935514368, -0.2086939711880792, -0.27519058110636668,
      -0.34459405529608833, -0.41582470131356863, -0.48661586339446783,
      -0.55145937253668, -0.60559923647246638, -0.64613906516596453,
      -0.67208074879767876, -0.68405422189518361, -0.68399274674594057,
      -0.67412560611764893, -0.6565811909641901, -0.63325620541364369,
      -0.60579315615624019, -0.57559657019301114, -0.54385924035378352,
      -0.51158769425763262, -0.47962381396642317, -0.44866239839012934,
      -0.41926543487440748, -0.3918737348579206, -0.36681682772623092,
      -0.34432143914061192, -0.32451895356138094, -0.3074523917588734,
      -0.2930829312678736, -0.28129658379511546, -0.2719114945174565,
      -0.26468660649245679, -0.25933278655428627, -0.25552781012667264,
      -0.25293656780947388, -0.25123723571073647, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.10602875205865551, 0.22266037932317656, 0.31888147181640641,
      0.38944360631144165, 0.43795507377677839, 0.46997264230390062,
      0.49051724877547076, 0.5034310014147434, 0.500324896800103,
      0.45828915737615084, 0.382639881027644, 0.28652002610423472,
      0.18334591217793197, 0.084623934905913356, -0.00043448287606375082,
      -0.069737181470766629, -0.1239972256853356, -0.1648520895179629,
      -0.19409834337546489, -0.21341793681113277, -0.22430732608084747,
      -0.22808296255519453, -0.22590849815057981, -0.21882346041009598,
      -0.20776651074652239, -0.19359407419059715, -0.17709264979243772,
      -0.1589888150121902, -0.13995640504519349, -0.12061994070409686,
      -0.10155785871994408, -0.083301402457821638, -0.066330226631092712,
      -0.051062749208288731, -0.037838754574905695, -0.0268921201815706,
      -0.018313911042383588, -0.012010230259510321, -0.0076711451969366674, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.42411500823462206, 0.46652650905808424,
      0.38488436997291947, 0.28224853798014093, 0.19404586986134692,
      0.12807027410848898, 0.082178425886280437, 0.051655010557090604,
      -0.012424418458561759, -0.16814295769580836, -0.30259710539402718,
      -0.38447941969363719, -0.41269645570521107, -0.39488790908807447,
      -0.34023367112790842, -0.27721079437881152, -0.2170401768582759,
      -0.16341945533050914, -0.11698501543000801, -0.077278373742671519,
      -0.043557557078858805, -0.015102545897388273, 0.0086978576184588621,
      0.028340150961935283, 0.044227798654294451, 0.056689746223700918,
      0.066005697592637733, 0.072415339120990019, 0.076129639867986862,
      0.077345857364386553, 0.076248327936611154, 0.0730258250484897,
      0.06788470330691572, 0.06106990969121591, 0.052895978533532173,
      0.043786537573340362, 0.034312836556748043, 0.025214723131493068,
      0.01735634025029462, 0.011464076937776417, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0038885244335512457, 0.010947898065094101, 0.020572442899544626,
      0.032244612643481642, 0.045513423888071919, 0.059974643590841557,
      0.075252099113920709, 0.090979651737706127, 0.10678324558559846,
      0.12226252352774972, 0.1369714974728323, 0.15039767726484851,
      0.16193908183164937, 0.17087846338538704, 0.17635400445530858,
      0.17779068774586007, 0.17611267036700712, 0.17207618251060741,
      0.16629705194676209, 0.1592740198771159, 0.15140835771064984,
      0.14302045975181585, 0.13436377780282663, 0.12563657342002277,
      0.11699173650140535, 0.10854503676053517, 0.10038201841561707,
      0.092563740004767608, 0.085131556583397425, 0.078111108115658287,
      0.07151559899571952, 0.0653484765221467, 0.059605704580110878,
      0.054277548891680738, 0.049350097732888551, 0.044806420719949967,
      0.040627568672192714, 0.036793323613225931, 0.033282836176213873, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.015554097734204983, 0.028237494526171424,
      0.0384981793378021, 0.046688678975748063, 0.053075244978361108,
      0.057844878811078554, 0.061109822092316614, 0.062910210495141658,
      0.06321437539156935, 0.061917111768604996, 0.058835895780330355,
      0.053704719168064845, 0.046165618267203466, 0.035757526214950611,
      0.021902164279686196, 0.0057467331622061074, -0.0067120695154118514,
      -0.0161459514255989, -0.023116522255381269, -0.028092128278584639,
      -0.031462648665864339, -0.033551591835335855, -0.034626727795956859,
      -0.034908817531215504, -0.034579347674469692, -0.033786798963480649,
      -0.032652073379672429, -0.03127311364339791, -0.029728733685480714,
      -0.02808179387095654, -0.026382036479755094, -0.024668489894291271,
      -0.022971087768143272, -0.021312622753720551, -0.019709804635168738,
      -0.018174708051754335, -0.016715408191029013, -0.015336980235867149,
      -0.014041949748048232, -0.012831045036862797, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    day4_helicopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    day4_helicopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    day4_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "day_4_test5.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(day4_helicopter_M,
                        "Error creating .mat file day_4_test5.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 8 + 1, 0, "state")) {
      rtmSetErrorStatus(day4_helicopter_M,
                        "Error writing mat file header to file day_4_test5.mat");
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
      char fileName[509] = "day_4_test5.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error closing MAT-file day_4_test5.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error reopening MAT-file day_4_test5.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 8 + 1,
           day4_helicopter_DW.ToFile3_IWORK.Count, "state")) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error writing header for state to MAT-file day_4_test5.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day4_helicopter_M,
                          "Error closing MAT-file day_4_test5.mat");
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

  rtmSetTFinal(day4_helicopter_M, 20.0);
  day4_helicopter_M->Timing.stepSize0 = 0.002;
  day4_helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  day4_helicopter_M->Sizes.checksums[0] = (1881219597U);
  day4_helicopter_M->Sizes.checksums[1] = (3690297916U);
  day4_helicopter_M->Sizes.checksums[2] = (759141438U);
  day4_helicopter_M->Sizes.checksums[3] = (2394120872U);

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
    day4_helicopter_B.Gain_b = 0.0;
    day4_helicopter_B.PitchCounttorad = 0.0;
    day4_helicopter_B.Gain_n = 0.0;
    day4_helicopter_B.Gain_k = 0.0;
    day4_helicopter_B.ElevationCounttorad = 0.0;
    day4_helicopter_B.Gain_d = 0.0;
    day4_helicopter_B.elevation = 0.0;
    day4_helicopter_B.Gain_i = 0.0;
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
