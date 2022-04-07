/*
 * day3_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day3_helicopter".
 *
 * Model version              : 11.14
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Mon Apr  4 15:54:26 2022
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "day3_helicopter.h"
#include "day3_helicopter_private.h"
#include "day3_helicopter_dt.h"

/* Block signals (default storage) */
B_day3_helicopter_T day3_helicopter_B;

/* Continuous states */
X_day3_helicopter_T day3_helicopter_X;

/* Block states (default storage) */
DW_day3_helicopter_T day3_helicopter_DW;

/* Real-time model */
static RT_MODEL_day3_helicopter_T day3_helicopter_M_;
RT_MODEL_day3_helicopter_T *const day3_helicopter_M = &day3_helicopter_M_;

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
  day3_helicopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void day3_helicopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum4[4];
  real_T rtb_Sum;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[6];
  real_T rtb_Gain1_idx_2;
  real_T rtb_Gain1_idx_3;
  real_T *lastU;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    /* set solver stop time */
    if (!(day3_helicopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&day3_helicopter_M->solverInfo,
                            ((day3_helicopter_M->Timing.clockTickH0 + 1) *
        day3_helicopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&day3_helicopter_M->solverInfo,
                            ((day3_helicopter_M->Timing.clockTick0 + 1) *
        day3_helicopter_M->Timing.stepSize0 +
        day3_helicopter_M->Timing.clockTickH0 *
        day3_helicopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(day3_helicopter_M)) {
    day3_helicopter_M->Timing.t[0] = rtsiGetT(&day3_helicopter_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(day3_helicopter_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: day3_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (day3_helicopter_DW.HILReadEncoderTimebase_Task, 1,
         &day3_helicopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          day3_helicopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          day3_helicopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          day3_helicopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    day3_helicopter_B.TravelCounttorad = day3_helicopter_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * day3_helicopter_P.TravelCounttorad_Gain;

    /* Gain: '<S13>/Gain' */
    day3_helicopter_B.Gain = day3_helicopter_P.Gain_Gain *
      day3_helicopter_B.TravelCounttorad;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/Constant'
     */
    day3_helicopter_B.Travel = day3_helicopter_P.Constant_Value +
      day3_helicopter_B.Gain;

    /* Gain: '<S4>/Pitch: Count to rad' */
    day3_helicopter_B.PitchCounttorad = day3_helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S10>/Gain' */
    day3_helicopter_B.Gain_i = day3_helicopter_P.Gain_Gain_a *
      day3_helicopter_B.PitchCounttorad;
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  rtb_Sum = 0.0;
  rtb_Sum += day3_helicopter_P.TravelTransferFcn_C *
    day3_helicopter_X.TravelTransferFcn_CSTATE;
  rtb_Sum += day3_helicopter_P.TravelTransferFcn_D *
    day3_helicopter_B.TravelCounttorad;

  /* Gain: '<S14>/Gain' */
  day3_helicopter_B.Gain_d = day3_helicopter_P.Gain_Gain_l * rtb_Sum;

  /* Sum: '<S3>/Sum' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  rtb_Sum = 0.0;
  rtb_Sum += day3_helicopter_P.PitchTransferFcn_C *
    day3_helicopter_X.PitchTransferFcn_CSTATE;
  rtb_Sum += day3_helicopter_P.PitchTransferFcn_D *
    day3_helicopter_B.PitchCounttorad;

  /* Gain: '<S11>/Gain' */
  day3_helicopter_B.Gain_b = day3_helicopter_P.Gain_Gain_ae * rtb_Sum;
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      day3_helicopter_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      day3_helicopter_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = day3_helicopter_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = day3_helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    day3_helicopter_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Sum = pDataValues[currTimeIndex];
        } else {
          rtb_Sum = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Sum = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    day3_helicopter_B.ElevationCounttorad = day3_helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * day3_helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S8>/Gain' */
    day3_helicopter_B.Gain_e = day3_helicopter_P.Gain_Gain_lv *
      day3_helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    day3_helicopter_B.Sum = day3_helicopter_B.Gain_e +
      day3_helicopter_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S9>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  day3_helicopter_B.Gain_dg = (day3_helicopter_P.ElevationTransferFcn_C *
    day3_helicopter_X.ElevationTransferFcn_CSTATE +
    day3_helicopter_P.ElevationTransferFcn_D *
    day3_helicopter_B.ElevationCounttorad) * day3_helicopter_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_2 = day3_helicopter_P.Gain1_Gain * day3_helicopter_B.Gain_i;
  rtb_Gain1_idx_3 = day3_helicopter_P.Gain1_Gain * day3_helicopter_B.Gain_b;

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      day3_helicopter_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      day3_helicopter_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = day3_helicopter_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = day3_helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    day3_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum4[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum4[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* Sum: '<S5>/Sum4' incorporates:
   *  Gain: '<S2>/Gain1'
   */
  rtb_Sum4[0] = day3_helicopter_P.Gain1_Gain * day3_helicopter_B.Travel -
    rtb_Sum4[0];
  rtb_Sum4[1] = day3_helicopter_P.Gain1_Gain * day3_helicopter_B.Gain_d -
    rtb_Sum4[1];
  rtb_Sum4[2] = rtb_Gain1_idx_2 - rtb_Sum4[2];
  rtb_Sum4[3] = rtb_Gain1_idx_3 - rtb_Sum4[3];

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S5>/K'
   *  Sum: '<S5>/Sum3'
   *  Sum: '<S5>/Sum4'
   *  Sum: '<S6>/Sum2'
   */
  rtb_Sum -= ((day3_helicopter_P.K[0] * rtb_Sum4[0] + day3_helicopter_P.K[1] *
               rtb_Sum4[1]) + day3_helicopter_P.K[2] * rtb_Sum4[2]) +
    day3_helicopter_P.K[3] * rtb_Sum4[3];
  rtb_Sum -= rtb_Gain1_idx_2;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum3'
   */
  day3_helicopter_B.Sum1 = (day3_helicopter_P.K_pp * rtb_Sum -
    day3_helicopter_P.K_pd * rtb_Gain1_idx_3) + day3_helicopter_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (day3_helicopter_X.Integrator_CSTATE >=
      day3_helicopter_P.Integrator_UpperSat) {
    day3_helicopter_X.Integrator_CSTATE = day3_helicopter_P.Integrator_UpperSat;
  } else {
    if (day3_helicopter_X.Integrator_CSTATE <=
        day3_helicopter_P.Integrator_LowerSat) {
      day3_helicopter_X.Integrator_CSTATE =
        day3_helicopter_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   *  Gain: '<S2>/Gain1'
   */
  rtb_Sum = day3_helicopter_P.elevation_ref_Value - day3_helicopter_P.Gain1_Gain
    * day3_helicopter_B.Sum;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Integrator: '<S3>/Integrator'
   *  Sum: '<S3>/Sum1'
   */
  day3_helicopter_B.Sum2 = ((day3_helicopter_P.K_ep * rtb_Sum +
    day3_helicopter_X.Integrator_CSTATE) - day3_helicopter_P.Gain1_Gain *
    day3_helicopter_B.Gain_dg * day3_helicopter_P.K_ed) +
    day3_helicopter_P.Vs_ff;
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    /* SignalConversion generated from: '<Root>/To File3' */
    rtb_TmpSignalConversionAtToFile[0] = day3_helicopter_B.Travel;
    rtb_TmpSignalConversionAtToFile[1] = day3_helicopter_B.Gain_d;
    rtb_TmpSignalConversionAtToFile[2] = day3_helicopter_B.Gain_i;
    rtb_TmpSignalConversionAtToFile[3] = day3_helicopter_B.Gain_b;
    rtb_TmpSignalConversionAtToFile[4] = day3_helicopter_B.Sum1;
    rtb_TmpSignalConversionAtToFile[5] = day3_helicopter_B.Sum2;

    /* ToFile: '<Root>/To File3' */
    {
      if (!(++day3_helicopter_DW.ToFile3_IWORK.Decimation % 1) &&
          (day3_helicopter_DW.ToFile3_IWORK.Count * (6 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) day3_helicopter_DW.ToFile3_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[6 + 1];
          day3_helicopter_DW.ToFile3_IWORK.Decimation = 0;
          u[0] = day3_helicopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          if (fwrite(u, sizeof(real_T), 6 + 1, fp) != 6 + 1) {
            rtmSetErrorStatus(day3_helicopter_M,
                              "Error writing to MAT-file day_3_test11.mat");
            return;
          }

          if (((++day3_helicopter_DW.ToFile3_IWORK.Count) * (6 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file day_3_test11.mat.\n");
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
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    rtAction = (int8_T)!(day3_helicopter_M->Timing.t[0] >= 2.0);
    day3_helicopter_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = day3_helicopter_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S7>/Action Port'
     */
    day3_helicopter_B.In1 = day3_helicopter_P.K_ei * rtb_Sum;
    if (rtmIsMajorTimeStep(day3_helicopter_M)) {
      srUpdateBC(day3_helicopter_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  rtb_Gain1_idx_2 = day3_helicopter_M->Timing.t[0];
  if ((day3_helicopter_DW.TimeStampA >= rtb_Gain1_idx_2) &&
      (day3_helicopter_DW.TimeStampB >= rtb_Gain1_idx_2)) {
    rtb_Gain1_idx_2 = 0.0;
  } else {
    rtb_Gain1_idx_3 = day3_helicopter_DW.TimeStampA;
    lastU = &day3_helicopter_DW.LastUAtTimeA;
    if (day3_helicopter_DW.TimeStampA < day3_helicopter_DW.TimeStampB) {
      if (day3_helicopter_DW.TimeStampB < rtb_Gain1_idx_2) {
        rtb_Gain1_idx_3 = day3_helicopter_DW.TimeStampB;
        lastU = &day3_helicopter_DW.LastUAtTimeB;
      }
    } else {
      if (day3_helicopter_DW.TimeStampA >= rtb_Gain1_idx_2) {
        rtb_Gain1_idx_3 = day3_helicopter_DW.TimeStampB;
        lastU = &day3_helicopter_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_2 = (day3_helicopter_B.PitchCounttorad - *lastU) /
      (rtb_Gain1_idx_2 - rtb_Gain1_idx_3);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S12>/Gain' */
  day3_helicopter_B.Gain_l = day3_helicopter_P.Gain_Gain_a1 * rtb_Gain1_idx_2;
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Gain1_idx_2 = (day3_helicopter_B.Sum2 - day3_helicopter_B.Sum1) *
    day3_helicopter_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Gain1_idx_2 > day3_helicopter_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day3_helicopter_B.BackmotorSaturation =
      day3_helicopter_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_2 < day3_helicopter_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day3_helicopter_B.BackmotorSaturation =
      day3_helicopter_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day3_helicopter_B.BackmotorSaturation = rtb_Gain1_idx_2;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_2 = (day3_helicopter_B.Sum1 + day3_helicopter_B.Sum2) *
    day3_helicopter_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Gain1_idx_2 > day3_helicopter_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day3_helicopter_B.FrontmotorSaturation =
      day3_helicopter_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_2 < day3_helicopter_P.FrontmotorSaturation_LowerSat)
  {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day3_helicopter_B.FrontmotorSaturation =
      day3_helicopter_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day3_helicopter_B.FrontmotorSaturation = rtb_Gain1_idx_2;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: day3_helicopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      day3_helicopter_DW.HILWriteAnalog_Buffer[0] =
        day3_helicopter_B.FrontmotorSaturation;
      day3_helicopter_DW.HILWriteAnalog_Buffer[1] =
        day3_helicopter_B.BackmotorSaturation;
      result = hil_write_analog(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILWriteAnalog_channels, 2,
        &day3_helicopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void day3_helicopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (day3_helicopter_DW.TimeStampA == (rtInf)) {
    day3_helicopter_DW.TimeStampA = day3_helicopter_M->Timing.t[0];
    lastU = &day3_helicopter_DW.LastUAtTimeA;
  } else if (day3_helicopter_DW.TimeStampB == (rtInf)) {
    day3_helicopter_DW.TimeStampB = day3_helicopter_M->Timing.t[0];
    lastU = &day3_helicopter_DW.LastUAtTimeB;
  } else if (day3_helicopter_DW.TimeStampA < day3_helicopter_DW.TimeStampB) {
    day3_helicopter_DW.TimeStampA = day3_helicopter_M->Timing.t[0];
    lastU = &day3_helicopter_DW.LastUAtTimeA;
  } else {
    day3_helicopter_DW.TimeStampB = day3_helicopter_M->Timing.t[0];
    lastU = &day3_helicopter_DW.LastUAtTimeB;
  }

  *lastU = day3_helicopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(day3_helicopter_M)) {
    rt_ertODEUpdateContinuousStates(&day3_helicopter_M->solverInfo);
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
  if (!(++day3_helicopter_M->Timing.clockTick0)) {
    ++day3_helicopter_M->Timing.clockTickH0;
  }

  day3_helicopter_M->Timing.t[0] = rtsiGetSolverStopTime
    (&day3_helicopter_M->solverInfo);

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
    if (!(++day3_helicopter_M->Timing.clockTick1)) {
      ++day3_helicopter_M->Timing.clockTickH1;
    }

    day3_helicopter_M->Timing.t[1] = day3_helicopter_M->Timing.clockTick1 *
      day3_helicopter_M->Timing.stepSize1 +
      day3_helicopter_M->Timing.clockTickH1 *
      day3_helicopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void day3_helicopter_derivatives(void)
{
  XDot_day3_helicopter_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_day3_helicopter_T *) day3_helicopter_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += day3_helicopter_P.TravelTransferFcn_A *
    day3_helicopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += day3_helicopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += day3_helicopter_P.PitchTransferFcn_A *
    day3_helicopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += day3_helicopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    day3_helicopter_P.ElevationTransferFcn_A *
    day3_helicopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += day3_helicopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (day3_helicopter_X.Integrator_CSTATE <=
          day3_helicopter_P.Integrator_LowerSat);
  usat = (day3_helicopter_X.Integrator_CSTATE >=
          day3_helicopter_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (day3_helicopter_B.In1 > 0.0)) || (usat &&
       (day3_helicopter_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = day3_helicopter_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void day3_helicopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day3_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &day3_helicopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(day3_helicopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(day3_helicopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
      return;
    }

    if ((day3_helicopter_P.HILInitialize_AIPStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &day3_helicopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (day3_helicopter_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &day3_helicopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = day3_helicopter_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_AIChannels, 8U,
        &day3_helicopter_DW.HILInitialize_AIMinimums[0],
        &day3_helicopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_AOPStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &day3_helicopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (day3_helicopter_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &day3_helicopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = day3_helicopter_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (day3_helicopter_DW.HILInitialize_Card,
         day3_helicopter_P.HILInitialize_AOChannels, 8U,
         &day3_helicopter_DW.HILInitialize_AOMinimums[0],
         &day3_helicopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_AOStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day3_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day3_helicopter_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_AOChannels, 8U,
        &day3_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day3_helicopter_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day3_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day3_helicopter_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (day3_helicopter_DW.HILInitialize_Card,
         day3_helicopter_P.HILInitialize_AOChannels, 8U,
         &day3_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_EIPStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &day3_helicopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = day3_helicopter_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (day3_helicopter_DW.HILInitialize_Card,
         day3_helicopter_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &day3_helicopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_EIStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &day3_helicopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = day3_helicopter_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_EIChannels, 8U,
        &day3_helicopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_POPStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day3_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day3_helicopter_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &day3_helicopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          day3_helicopter_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &day3_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            day3_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            day3_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = day3_helicopter_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            day3_helicopter_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            day3_helicopter_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] = day3_helicopter_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(day3_helicopter_DW.HILInitialize_Card,
          &day3_helicopter_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes, &day3_helicopter_DW.HILInitialize_POSortedFreqs
          [0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(day3_helicopter_DW.HILInitialize_Card,
          &day3_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &day3_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day3_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day3_helicopter_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &day3_helicopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = day3_helicopter_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &day3_helicopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = day3_helicopter_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &day3_helicopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &day3_helicopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &day3_helicopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &day3_helicopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = day3_helicopter_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &day3_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day3_helicopter_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_POChannels, 8U,
        &day3_helicopter_DW.HILInitialize_POSortedFreqs[0],
        &day3_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day3_helicopter_P.HILInitialize_POStart && !is_switching) ||
        (day3_helicopter_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day3_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day3_helicopter_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(day3_helicopter_DW.HILInitialize_Card,
        day3_helicopter_P.HILInitialize_POChannels, 8U,
        &day3_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day3_helicopter_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &day3_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day3_helicopter_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (day3_helicopter_DW.HILInitialize_Card,
         day3_helicopter_P.HILInitialize_POChannels, 8U,
         &day3_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: day3_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (day3_helicopter_DW.HILInitialize_Card,
       day3_helicopter_P.HILReadEncoderTimebase_SamplesI,
       day3_helicopter_P.HILReadEncoderTimebase_Channels, 3,
       &day3_helicopter_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (day3_helicopter_DW.HILReadEncoderTimebase_Task, (t_buffer_overflow_mode)
         (day3_helicopter_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
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
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.32964146257353211,
      -0.08219588820755408, -0.390160985441269, -0.52359877559825407,
      -0.52359877559829737, -0.52359877559829848, -0.52359877559829859,
      -0.52359877559829859, -0.52359877559829859, -0.5235987755982987,
      -0.52359877559829615, -0.52359877559829837, -0.49033642936304667,
      -0.40501787444629567, -0.32447042262326431, -0.25079671707513257,
      -0.18530807999911131, -0.12865835655490143, -0.080975487106326127,
      -0.041984628103399146, -0.01111884795734297, 0.012384768741018348,
      0.02940439420098151, 0.040867891473357787, 0.047701505304608816,
      0.050791466539744423, 0.050956920530022543, 0.0489324893527906,
      0.045358784036187316, 0.040779265106604679, 0.035641982497751856,
      0.03030488920771679, 0.025043601465379983, 0.020060659666798375,
      0.015495520179333311, 0.011434672106179211, 0.007921421226368297,
      0.0049650132858422014, 0.0025488796997533703, 0.000637880673778457,
      -0.0008155053251506228, -0.0018660394743228448, -0.0025715790555206341,
      -0.002989853942904297, -0.0031760657247065227, -0.0031812020854021794,
      -0.0030509561505469573, -0.0028251429618957324, -0.0025375116354993121,
      -0.0022158607974401834, -0.0018823755361014882, -0.0015541155015651453,
      -0.0012435952523384142, -0.00095940899401347757, -0.00070686211605353044,
      -0.00048858117774264187, -0.00030508209065782665, -0.00015528314289048328,
      -3.695522289504094E-5, 5.2893811869125784E-5, 0.00011770012233114358,
      0.00016108709903484009, 0.0001866658303403268, 0.00019788669932541136,
      0.00019793537119261462, 0.00018966626054817137, 0.00017556672950169716,
      0.000157745671994447, 0.00013794070900463584, 0.00011753888834820092,
      9.760649933043819E-5, 7.8924335932928358E-5, 6.2025442646218565E-5,
      4.723303412768054E-5, 3.469688163670881E-5, 2.442700111449092E-5,
      1.6323961813502486E-5, 1.0205568148236033E-5, 5.8300631057273122E-6,
      2.9163733161752248E-6, 1.1622739847227948E-6, 2.616877995409439E-7,
      -7.7416244748107488E-8, -1.130834783591439E-7, -4.6748792503770176E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    day3_helicopter_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    day3_helicopter_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    day3_helicopter_DW.FromWorkspace_IWORK.PrevIndex = 0;
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
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625142, 3.1262155534579494,
      3.1033093000306571, 3.066627415191995, 3.0144539223957096,
      2.9456562771193773, 2.8595077632960466, 2.7555515879674664,
      2.6335051104928389, 2.4931956060348854, 2.3345185760674338,
      2.1574113214743651, 1.9632257907080559, 1.7564180641970282,
      1.5436868048209913, 1.33200838573368, 1.1270304718834845,
      0.93285514183342755, 0.75228626335724291, 0.58716721224025148,
      0.438677531842949, 0.30755709868998216, 0.19426296508845742,
      0.099074075932674491, 0.021920755730521019, -0.037969151990778535,
      -0.082000899945563588, -0.11201440990535054, -0.13007759875560709,
      -0.13831437791177367, -0.13877705511435465, -0.1333599752794305,
      -0.12374718478502027, -0.11138645055600013, -0.097482802921114886,
      -0.083005938211487149, -0.068706949557068436, -0.055140839500328859,
      -0.042692090304844006, -0.031601243290898258, -0.021990989750957143,
      -0.013890723488906142, -0.0072588654709380843, -0.0020025577978081588,
      0.0020054519065688292, 0.0049127443079569257, 0.0068758641642188917,
      0.0080514552122598615, 0.0085896360162269619, 0.0086293288706718891,
      0.008295244079715788, 0.0076962273758151975, 0.0069246947298248823,
      0.0060569028000332873, 0.005153831824550697, 0.0042624885086320237,
      0.003417467531829468, 0.0026426403048120306, 0.0019528675435259673,
      0.0013556574450566791, 0.00085271337060279439, 0.00044133381876016532,
      0.00011564313272294466, -0.00013235601625935194, -0.00031190754487320503,
      -0.00043279618170196369, -0.00050479782422356329, -0.000537269614873395,
      -0.0005388615962546796, -0.00051733128213879442, -0.00047944287683303121,
      -0.00043093394137479446, -0.00037653382551175717, -0.00032001998052576914,
      -0.00026430019372196672, -0.00021151072665373242, -0.00016312221000391115,
      -0.00012004688717233426, -8.2742365870920245E-5, -5.1308409597692614E-5,
      -2.5574469893841769E-5, -5.1766277601633177E-6, 1.0376611513982647E-5,
      2.1649626028323855E-5, 2.9235600482817835E-5, 3.3722270965622214E-5,
      3.566650464391265E-5, 3.55766278837951E-5, 3.39013890405747E-5,
      3.10244683319768E-5, 2.7263509146696822E-5, 2.2872726563081076E-5,
      1.8048234569211122E-5, 1.2935309345237139E-5, 7.6368591935275476E-6,
      2.2223926847670631E-6, -3.2632356668009043E-6, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048908026831,
      -0.046506351615507022, -0.09162501371031305, -0.14672753935593108,
      -0.2086939711829455, -0.27519058110598277, -0.34459405529323134,
      -0.41582470131328692, -0.48818590989835442, -0.56123801783305816,
      -0.63470811986957276, -0.70842901837219463, -0.77674212306588908,
      -0.82723090604333982, -0.85092503750380555, -0.84671367634979877,
      -0.8199116554015129, -0.77670132020469806, -0.722275513902851,
      -0.6604762044673389, -0.59395872158533325, -0.52448173261770192,
      -0.4531765344093604, -0.38075555661873917, -0.30861328081230216,
      -0.23955963088805216, -0.17612699181998709, -0.12005403983921872,
      -0.072252755405057392, -0.032947116624265152, -0.0018507088076270112,
      0.021668319337823125, 0.03845116197788076, 0.049442936916257681,
      0.05561459054030217, 0.057907458838801912, 0.05719595461695607,
      0.054264440226942129, 0.04979499678208827, 0.044363388055817438,
      0.038441014159713217, 0.032401065048490552, 0.026527432071689054,
      0.021025230692366819, 0.016032038817462409, 0.011629169605786875,
      0.0078524794250010663, 0.004702364192118713, 0.0021527232161369951,
      0.0001587714175167645, -0.0013363391637847089, -0.0023960668154646191,
      -0.0030861305839162624, -0.0034711677191495509, -0.0036122839019729213,
      -0.0035653732636909013, -0.0033800839072083136, -0.0030993089080704268,
      -0.0027590910451350584, -0.0023888403938751783, -0.0020117762978097628,
      -0.0016455182073792612, -0.0013027627441365012, -0.00099199659593792681,
      -0.0007182061144461806, -0.00048355454731884489, -0.00028800657006879652,
      -0.00012988716260278261, -6.3679255287958441E-6, 8.6121256461146967E-5,
      0.00015155362122463279, 0.00019403574182458194, 0.00021760046345421341,
      0.00022605537994435388, 0.00022287914721675579, 0.00021115786827275397,
      0.00019355406659912557, 0.00017230129132586483, 0.00014921808520616459,
      0.00012573582509283499, 0.00010293575881514652, 8.15913685345998E-5,
      6.2212957097588351E-5, 4.5092058057169079E-5, 3.0343897817401947E-5,
      1.7946681931419631E-5, 7.7769347120128763E-6, -3.5950704050567589E-7,
      -6.7009553726487046E-6, -1.1507682834820927E-5, -1.504383674063341E-5,
      -1.7563130334661396E-5, -1.9297967975376671E-5, -2.0451700895911528E-5,
      -2.11938006068429E-5, -2.1657866034981538E-5, -2.194251340628984E-5,
      -2.211524572262389E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875206161289, 0.22266037932533145, 0.31888147181725313,
      0.38944360631191, 0.43795507377723664, 0.46997264230406444,
      0.49051724877546532, 0.50343100141464336, 0.51142138586016028,
      0.51630439857690191, 0.51925862127057065, 0.52103115488679352,
      0.48281093361230909, 0.35683543525828809, 0.16746107181424308,
      -0.029764292226058309, -0.18942644773320871, -0.30539414616636545,
      -0.38466081240189204, -0.43677391648549591, -0.4701201644410804,
      -0.49103682320697251, -0.503957907842915, -0.5118438116659112,
      -0.50987405265616059, -0.48804482446838804, -0.4483176665678606,
      -0.39630220903061592, -0.33784122199671474, -0.27779724316647747,
      -0.21977753396647637, -0.16622350849758255, -0.11861472203408568,
      -0.0776856672462286, -0.043618890731051352, -0.016205117438777505,
      0.0050286400876812243, 0.020718824096048416, 0.031588319286347311,
      0.038388536019832609, 0.041857076805087656, 0.042688053523619596,
      0.041512594598347219, 0.038887457931902625, 0.03528997315590765,
      0.031117798029808519, 0.026692203789390634, 0.022263811366006858,
      0.018019888654750238, 0.014092489772381311, 0.010566870570910769,
      0.00748975030825938, 0.0048771071629601925, 0.0027212954175485082,
      0.0009973552837373445, -0.00033154647512501079, -0.0013095544052306174,
      -0.0019844104593228495, -0.0024045299181540969, -0.0026167901958079609,
      -0.0026649450213495074, -0.0025885723007217365, -0.0024224647078952,
      -0.0021963764467941571, -0.0019350465559391239, -0.0016584276576662171,
      -0.0013820584187647089, -0.0011175275827104425, -0.00087298679265523926,
      -0.00065367821458828246, -0.00046245096391184237, -0.00030024740348943979,
      -0.00016654645248281863, -5.97561206211461E-5, 2.2448399836449973E-5,
      8.28415229217061E-5, 0.00012441694693010596, 0.00015020649870067793,
      0.00016314328484723362, 0.00016596364605425684, 0.0001611421605662855,
      0.00015085399858161264, 0.00013695921098366703, 0.00012100397553638231,
      0.00010423436390327723, 8.761878710528137E-5, 7.1875889265826665E-5,
      5.7505262818891012E-5, 4.4818934865187288E-5, 3.3972113898661505E-5,
      2.4992185271077716E-5, 1.7805404940229863E-5, 1.2261169871741195E-5,
      8.1541436438570258E-6, 5.2448773324442755E-6, 3.2798372090470096E-6,
      2.0117789064633129E-6, 1.2208060413687605E-6, 7.343235263190806E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823508058,
      0.46652650905814247, 0.38488436997386222, 0.28224853797840926,
      0.19404586986032613, 0.12807027410748784, 0.08217842588557929,
      0.05165501055653117, 0.03196153778197907, 0.019532050866961467,
      0.011816890774579094, 0.0070901344647471137, -0.15288088509774941,
      -0.50390199341626907, -0.75749745377558053, -0.7889014561708747,
      -0.63864862202150308, -0.46387079373309209, -0.3170666649445576,
      -0.20845241633388606, -0.13338499182224864, -0.083666635063545991,
      -0.051684338543796442, -0.031543615292220328, 0.0078790360394324164,
      0.087316912751313883, 0.15890863160249885, 0.20806183014870186,
      0.23384394813579815, 0.24017591531969343, 0.23207883679269439,
      0.21421610188103546, 0.1904351458553725, 0.16371621914797746,
      0.13626710606500161, 0.1096550931694316, 0.084935030105616566,
      0.062760736033085573, 0.043477980761461338, 0.027200866934187728,
      0.01387416314071936, 0.0033239068743258168, -0.0047018357015254364,
      -0.010500546666055665, -0.014389939104282511, -0.01668870050512105,
      -0.017702376961568317, -0.017713569694124374, -0.016975690845386818,
      -0.015709595529411355, -0.014102476806453119, -0.012308481050571326,
      -0.010450572580858622, -0.0086232469819272628, -0.0068957605352018836,
      -0.0053156070354806583, -0.0039120317203783107, -0.002699424216354896,
      -0.0016804778353419867, -0.000849041110622517, -0.00019261930211289361,
      0.00030549088246541374, 0.00066443037128795809, 0.00090435304444901879,
      0.0010453195634504979, 0.0011064755930951151, 0.0011054769556172765,
      0.0010581233443033224, 0.00097816316020312983, 0.00087723431227038258,
      0.00076490900273029432, 0.00064881424170254833, 0.00053480380402744384,
      0.00042716132744609694, 0.00032881808183159817, 0.00024157249234016921,
      0.00016630169603238003, 0.00010315820708070762, 5.1747144591242731E-5,
      1.1281444827790785E-5, -1.9285941949462666E-5, -4.1152647936726315E-5,
      -5.5579150395020694E-5, -6.3820941785563227E-5, -6.7078446530984779E-5,
      -6.6462307194093976E-5, -6.2971591355575148E-5, -5.7482505787841469E-5,
      -5.0745311815288416E-5, -4.338728386479286E-5, -3.5919714511550457E-5,
      -2.8747121323016721E-5, -2.2176940274054065E-5, -1.6428104911343013E-5,
      -1.1637065245432845E-5, -7.8601604934477563E-6, -5.0722332106450252E-6,
      -3.1638914604558217E-6, -1.9459300603116726E-6, -1.1834458995398634E-6,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    day3_helicopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    day3_helicopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    day3_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "day_3_test11.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(day3_helicopter_M,
                        "Error creating .mat file day_3_test11.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 6 + 1, 0, "state")) {
      rtmSetErrorStatus(day3_helicopter_M,
                        "Error writing mat file header to file day_3_test11.mat");
      return;
    }

    day3_helicopter_DW.ToFile3_IWORK.Count = 0;
    day3_helicopter_DW.ToFile3_IWORK.Decimation = -1;
    day3_helicopter_DW.ToFile3_PWORK.FilePtr = fp;
  }

  /* Start for If: '<S3>/If' */
  day3_helicopter_DW.If_ActiveSubsystem = -1;

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  day3_helicopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  day3_helicopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  day3_helicopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  day3_helicopter_X.Integrator_CSTATE = day3_helicopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  day3_helicopter_DW.TimeStampA = (rtInf);
  day3_helicopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void day3_helicopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day3_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(day3_helicopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(day3_helicopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((day3_helicopter_P.HILInitialize_AOTerminate && !is_switching) ||
        (day3_helicopter_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day3_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day3_helicopter_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((day3_helicopter_P.HILInitialize_POTerminate && !is_switching) ||
        (day3_helicopter_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day3_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day3_helicopter_P.HILInitialize_POFinal;
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
      result = hil_write(day3_helicopter_DW.HILInitialize_Card
                         , day3_helicopter_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , day3_helicopter_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &day3_helicopter_DW.HILInitialize_AOVoltages[0]
                         , &day3_helicopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(day3_helicopter_DW.HILInitialize_Card,
            day3_helicopter_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &day3_helicopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(day3_helicopter_DW.HILInitialize_Card,
            day3_helicopter_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &day3_helicopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day3_helicopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(day3_helicopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(day3_helicopter_DW.HILInitialize_Card);
    hil_close(day3_helicopter_DW.HILInitialize_Card);
    day3_helicopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (FILE *) day3_helicopter_DW.ToFile3_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "day_3_test11.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error closing MAT-file day_3_test11.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error reopening MAT-file day_3_test11.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 6 + 1,
           day3_helicopter_DW.ToFile3_IWORK.Count, "state")) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error writing header for state to MAT-file day_3_test11.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error closing MAT-file day_3_test11.mat");
        return;
      }

      day3_helicopter_DW.ToFile3_PWORK.FilePtr = (NULL);
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
  day3_helicopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  day3_helicopter_update();
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
  day3_helicopter_initialize();
}

void MdlTerminate(void)
{
  day3_helicopter_terminate();
}

/* Registration function */
RT_MODEL_day3_helicopter_T *day3_helicopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  day3_helicopter_P.Integrator_UpperSat = rtInf;
  day3_helicopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)day3_helicopter_M, 0,
                sizeof(RT_MODEL_day3_helicopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&day3_helicopter_M->solverInfo,
                          &day3_helicopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&day3_helicopter_M->solverInfo, &rtmGetTPtr(day3_helicopter_M));
    rtsiSetStepSizePtr(&day3_helicopter_M->solverInfo,
                       &day3_helicopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&day3_helicopter_M->solverInfo, &day3_helicopter_M->derivs);
    rtsiSetContStatesPtr(&day3_helicopter_M->solverInfo, (real_T **)
                         &day3_helicopter_M->contStates);
    rtsiSetNumContStatesPtr(&day3_helicopter_M->solverInfo,
      &day3_helicopter_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&day3_helicopter_M->solverInfo,
      &day3_helicopter_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&day3_helicopter_M->solverInfo,
      &day3_helicopter_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&day3_helicopter_M->solverInfo,
      &day3_helicopter_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&day3_helicopter_M->solverInfo, (&rtmGetErrorStatus
      (day3_helicopter_M)));
    rtsiSetRTModelPtr(&day3_helicopter_M->solverInfo, day3_helicopter_M);
  }

  rtsiSetSimTimeStep(&day3_helicopter_M->solverInfo, MAJOR_TIME_STEP);
  day3_helicopter_M->intgData.f[0] = day3_helicopter_M->odeF[0];
  day3_helicopter_M->contStates = ((real_T *) &day3_helicopter_X);
  rtsiSetSolverData(&day3_helicopter_M->solverInfo, (void *)
                    &day3_helicopter_M->intgData);
  rtsiSetSolverName(&day3_helicopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = day3_helicopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    day3_helicopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    day3_helicopter_M->Timing.sampleTimes =
      (&day3_helicopter_M->Timing.sampleTimesArray[0]);
    day3_helicopter_M->Timing.offsetTimes =
      (&day3_helicopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    day3_helicopter_M->Timing.sampleTimes[0] = (0.0);
    day3_helicopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    day3_helicopter_M->Timing.offsetTimes[0] = (0.0);
    day3_helicopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(day3_helicopter_M, &day3_helicopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = day3_helicopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    day3_helicopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(day3_helicopter_M, 25.0);
  day3_helicopter_M->Timing.stepSize0 = 0.002;
  day3_helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  day3_helicopter_M->Sizes.checksums[0] = (657371631U);
  day3_helicopter_M->Sizes.checksums[1] = (410321189U);
  day3_helicopter_M->Sizes.checksums[2] = (121825746U);
  day3_helicopter_M->Sizes.checksums[3] = (1249057324U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    day3_helicopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &day3_helicopter_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(day3_helicopter_M->extModeInfo,
      &day3_helicopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(day3_helicopter_M->extModeInfo,
                        day3_helicopter_M->Sizes.checksums);
    rteiSetTPtr(day3_helicopter_M->extModeInfo, rtmGetTPtr(day3_helicopter_M));
  }

  day3_helicopter_M->solverInfoPtr = (&day3_helicopter_M->solverInfo);
  day3_helicopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&day3_helicopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&day3_helicopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  day3_helicopter_M->blockIO = ((void *) &day3_helicopter_B);

  {
    day3_helicopter_B.TravelCounttorad = 0.0;
    day3_helicopter_B.Gain = 0.0;
    day3_helicopter_B.Travel = 0.0;
    day3_helicopter_B.Gain_d = 0.0;
    day3_helicopter_B.PitchCounttorad = 0.0;
    day3_helicopter_B.Gain_i = 0.0;
    day3_helicopter_B.Gain_b = 0.0;
    day3_helicopter_B.ElevationCounttorad = 0.0;
    day3_helicopter_B.Gain_e = 0.0;
    day3_helicopter_B.Sum = 0.0;
    day3_helicopter_B.Gain_dg = 0.0;
    day3_helicopter_B.Sum1 = 0.0;
    day3_helicopter_B.Sum2 = 0.0;
    day3_helicopter_B.Gain_l = 0.0;
    day3_helicopter_B.BackmotorSaturation = 0.0;
    day3_helicopter_B.FrontmotorSaturation = 0.0;
    day3_helicopter_B.In1 = 0.0;
  }

  /* parameters */
  day3_helicopter_M->defaultParam = ((real_T *)&day3_helicopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &day3_helicopter_X;
    day3_helicopter_M->contStates = (x);
    (void) memset((void *)&day3_helicopter_X, 0,
                  sizeof(X_day3_helicopter_T));
  }

  /* states (dwork) */
  day3_helicopter_M->dwork = ((void *) &day3_helicopter_DW);
  (void) memset((void *)&day3_helicopter_DW, 0,
                sizeof(DW_day3_helicopter_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day3_helicopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  day3_helicopter_DW.TimeStampA = 0.0;
  day3_helicopter_DW.LastUAtTimeA = 0.0;
  day3_helicopter_DW.TimeStampB = 0.0;
  day3_helicopter_DW.LastUAtTimeB = 0.0;
  day3_helicopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  day3_helicopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    day3_helicopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  day3_helicopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  day3_helicopter_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  day3_helicopter_M->Sizes.numY = (0); /* Number of model outputs */
  day3_helicopter_M->Sizes.numU = (0); /* Number of model inputs */
  day3_helicopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  day3_helicopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  day3_helicopter_M->Sizes.numBlocks = (66);/* Number of blocks */
  day3_helicopter_M->Sizes.numBlockIO = (17);/* Number of block outputs */
  day3_helicopter_M->Sizes.numBlockPrms = (149);/* Sum of parameter "widths" */
  return day3_helicopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
