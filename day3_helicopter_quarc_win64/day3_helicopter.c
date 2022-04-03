/*
 * day3_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day3_helicopter".
 *
 * Model version              : 11.10
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Sun Apr  3 20:55:07 2022
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
   *  Gain: '<S5>/Gain'
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
                              "Error writing to MAT-file day_3_test10.mat");
            return;
          }

          if (((++day3_helicopter_DW.ToFile3_IWORK.Count) * (6 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file day_3_test10.mat.\n");
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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877560561874,
      0.52359877560243773, 0.52359877559967882, 0.52359877560176293,
      0.52359877560250567, 0.52359877559879708, 0.523598775601886,
      0.48789589211690976, 0.34971187068198772, 0.22862453939368321,
      0.12361824412228262, 0.033621435678587663, -0.042474371213739148,
      -0.10579616716527829, -0.1574723836113574, -0.19861769233059073,
      -0.23031990327360963, -0.25362901646162456, -0.2695483693600485,
      -0.27902775197638452, -0.28295832332992721, -0.282169144438834,
      -0.2774251374368274, -0.26942628311931704, -0.25880787701330971,
      -0.24614167493421862, -0.2319377716052306, -0.21664706939463596,
      -0.20066420798890355, -0.18433083946977491, -0.16793914653320974,
      -0.1517355142929353, -0.13592427812631835, -0.12067148126441607,
      -0.10610858625019387, -0.092336093964110078, -0.07942703263727624,
      -0.0674302871475923, -0.056373745943501641, -0.046267249192563642,
      -0.0371053272442784, -0.028869723270250169, -0.02153169804541244,
      -0.015054118309729803, -0.0093933330500299839, -0.0045008444163011641,
      -0.00032478188550688625, 0.0031888102427672926, 0.0060948608437053631,
      0.0084483235272315715, 0.0103033793675813, 0.011712778199718943,
      0.012727305885465223, 0.013395365138720638, 0.013762657834992087,
      0.013871957186847195, 0.013762958718040452, 0.013472199591795242,
      0.01303303652261667, 0.012475673208138915, 0.011827228942350554,
      0.011111840800780426, 0.010350792510642193, 0.0095626638251978813,
      0.0087634949041146548, 0.0079669608543352366, 0.007184552204344663,
      0.0064257576653433679, 0.00569824607345959, 0.0050080449064697774,
      0.0043597132260854154, 0.00375650731295063, 0.0032005376369405036,
      0.0026929161415080172, 0.0022338931194832989, 0.0018229832210291397,
      0.0014590803648616424, 0.0011405615241664169, 0.00086537953202614482,
      0.000631145201325567, 0.00043519918545331659, 0.0002746741243856432,
      0.00014654773359501494, 4.7687611652746931E-5, -2.51113174962736E-5,
      -7.5095625605481331E-5, -0.00010552709432021246, -0.00011965954800363665,
      -0.00012071641349342688, -0.00011186568160814758, -9.6187572839689928E-5,
      -7.6628292929425967E-5, -5.5930824064076567E-5, -3.653115928026196E-5,
      -2.0407226918428556E-5, -8.87193000354003E-6, -2.3207706419947627E-6,
      2.2423174428354287E-12, 2.2422064205329661E-12, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

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
      3.1415926535897931, 3.1378421413628796, 3.1262155534590508,
      3.1033093000312753, 3.0666274151926292, 3.014453922395651,
      2.9456562771196104, 2.8595077632954782, 2.7558073259615177,
      2.6355434416064383, 2.5007315960809082, 2.3540419646005355,
      2.1984466981367543, 2.0369525052709658, 1.8724212606771795,
      1.70746165134267, 1.5443722376840752, 1.3851193897200611,
      1.2313377619747634, 1.084344677034049, 0.94516260709943922,
      0.81454593599937186, 0.69300953525459641, 0.58085758134157262,
      0.4782116243499428, 0.38503729687455357, 0.30116929848631241,
      0.22633445299366967, 0.16017274341930116, 0.10225630243745035,
      0.0521063860116434, 0.0092083926790348346, -0.026974984917693443,
      -0.056992373034600988, -0.081393980110083461, -0.10072338928314997,
      -0.11551075563777662, -0.12626728283685437, -0.13348085509340527,
      -0.13761270339363377, -0.13909498918597024, -0.13832919408585967,
      -0.13568521025758057, -0.13150103280621078, -0.12608296255207607,
      -0.11970623480413517, -0.11261599705906473, -0.10502856581238532,
      -0.097132899778658077, -0.089092233701274892, -0.081045823524040772,
      -0.073110759948836535, -0.065383813279096856, -0.057943277923115745,
      -0.050850789988729807, -0.044153096035537014, -0.037883755262336727,
      -0.032064761202882382, -0.026708072393367165, -0.021817044476218284,
      -0.017387758834953089, -0.013410245135285638, -0.0098695971009826,
      -0.0067469825029098946, -0.0040205497106447548, -0.0016662342724360706,
      0.00034153012458697707, 0.0020291912825299574, 0.0034235431819580023,
      0.0045512698710344478, 0.0054385625005659976, 0.006110803379750702,
      0.0065923109598145643, 0.0069061397719682423, 0.0070739295367997371,
      0.0071157979105585792, 0.0070502716269881137, 0.0068942511204199966,
      0.0066630040664217216, 0.0063701836417806778, 0.0060278676780718453,
      0.0056466152557060282, 0.0052355376521854248, 0.0048023809139604605,
      0.0043536176608643395, 0.0038945460508059537, 0.0034293941251199531,
      0.0029614280156192506, 0.0024930627148619109, 0.0020259742800656116,
      0.0015612124415638952, 0.0010993125938950771, 0.0006404060254672746,
      0.00018432694173208815, -0.00026928570489996481, -0.00072089440885930365,
      -0.0011709813238870604, -0.0016199838401650167, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.01500204890656333,
      -0.046506351616179109, -0.091625013709456624, -0.14672753935569904,
      -0.20869397118820723, -0.27519058110691275, -0.34459405529640907,
      -0.41480174933452263, -0.48105553742129337, -0.53924738210369716,
      -0.58675852592175626, -0.62238106585451269, -0.6459767714626985,
      -0.65812497837516559, -0.659838437337216, -0.65235765463537054,
      -0.63701139185646583, -0.61512651098007154, -0.58797233976318208,
      -0.55672827973842087, -0.52246668440016464, -0.4861456029790247,
      -0.448607815652071, -0.4105838279663982, -0.37269730990166317,
      -0.33547199355304, -0.29933938197059484, -0.26464683829742924,
      -0.23166576392741245, -0.20059966570318619, -0.17159197333043966,
      -0.14473351038684309, -0.12006955246778667, -0.097606428301749162,
      -0.077317636692324332, -0.059149465418467453, -0.043026108796262758,
      -0.028854289026261422, -0.016527393200825153, -0.0059291431693221566,
      0.0030631804004240485, 0.010575935313022101, 0.016736709805491516,
      0.021672281016505684, 0.0255069109917853, 0.02836095098028896,
      0.030349724986739629, 0.031582664134881734, 0.032162664309489712,
      0.03218564070897284, 0.031740254300854573, 0.030907786678922008,
      0.029762141423938078, 0.028369951737531367, 0.026790775812789543,
      0.025077363092780035, 0.023275976237823973, 0.021426755238051147,
      0.019564111668596267, 0.017717142565056784, 0.015910054798669547,
      0.014162592137211633, 0.012490458392291244, 0.010905731169060752,
      0.00941726175283719, 0.0080310575880885166, 0.0067506446317744436,
      0.0055774075977103616, 0.0045109067563033848, 0.0035491705181318017,
      0.0026889635167396393, 0.0019260303202602554, 0.0012553152486133971,
      0.00067115905932383843, 0.00016747349503473769, -0.0002621051342807392,
      -0.000624082026274541, -0.00092498821599404791, -0.0011712816985660695,
      -0.0013692638548292639, -0.0015250096894597904, -0.0016443104140825373,
      -0.0017326269529012465, -0.0017950530123850362, -0.0018362864402340882,
      -0.0018606077027435496, -0.0018718644380027933, -0.0018734612030289807,
      -0.0018683537391858392, -0.001859047354006321, -0.0018475993906754695,
      -0.0018356262737114278, -0.0018243163349405144, -0.0018144505865283504,
      -0.0018064348158372005, -0.0018003476601109931, -0.0017960100651118315,
      -0.0017930800755657408, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875206157947, 0.22266037932492488, 0.31888147181824933,
      0.3894436063137956, 0.43795507377952358, 0.46997264230631586,
      0.49051724877798242, 0.49620116751189258, 0.46825647032084605,
      0.41127773340780038, 0.3357906189114076, 0.25176650717163274,
      0.16676543549643508, 0.085858886776146726, 0.012110073534532972,
      -0.05287131505877718, -0.1084615242217622, -0.15467397967191054,
      -0.19191531135768641, -0.22082108341782736, -0.24214787054821191,
      -0.25670353161802395, -0.26530274426991723, -0.26873902271412292,
      -0.26776743993909324, -0.263094318736156, -0.25537149877825033,
      -0.24519364879080907, -0.23309763740724576, -0.2195633173813864,
      -0.20501529097683346, -0.18982535821095181, -0.17431543483204842,
      -0.15876078241253483, -0.1433934303325945, -0.1284056957139193,
      -0.11395372671868254, -0.10016100958442298, -0.087121791764035172,
      -0.0749043835025981, -0.0635543086112616, -0.05309728242509032,
      -0.043542001166887168, -0.034882732307646047, -0.027101700128050854,
      -0.020171264612098061, -0.014055895116430461, -0.0087139430087471714,
      -0.0040992197160493182, -0.00016238841325832976, 0.003147821031485476,
      0.0058835632150775474, 0.00809698311637963, 0.0098394649960453773,
      0.011161012314529262, 0.01210974671500864, 0.012731514301768554,
      0.013069587767997826, 0.013164453364339468, 0.013053672225603608,
      0.012771806166277377, 0.012350398696818643, 0.011818002684607198,
      0.011200246771879807, 0.010519933355056965, 0.009797161615069494,
      0.0090494697580514316, 0.0082919912725331058, 0.0075376206278430224,
      0.006797184424246927, 0.0060796145550342517, 0.005392120453069027,
      0.0047403579669591522, 0.004128592845259349, 0.0035598572014762242,
      0.0030360976881895763, 0.0025583144274451408, 0.0021266900281566858,
      0.0017407082719027578, 0.0013992622688442369, 0.0011007520785315528,
      0.0008431719602477461, 0.00062418756856874591, 0.00044120354880905488,
      0.00029142212240529375, 0.0001718933969463432, 7.9558308014426871E-5,
      1.1285327767729392E-5, -3.609760136191742E-5, -6.5773973673688424E-5,
      -8.0909830129094473E-5, -8.4621415423868918E-5, -7.9934325225017311E-5,
      -6.9727339336189686E-5, -5.6652403340295976E-5, -4.3021689350775105E-5,
      -3.0656462330469481E-5, -2.0708047080764125E-5, -1.3505271479918868E-5,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500824011378,
      0.46652650906198218, 0.38488436997254627, 0.28224853798169158,
      0.19404586986348327, 0.1280702741068617, 0.082178425886042836,
      0.022735674934001696, -0.11177878876619858, -0.22791494765357151,
      -0.30194845798715719, -0.33609644696036334, -0.3400042866982827,
      -0.32362619488350647, -0.29499525296694445, -0.25992555437393422,
      -0.22236083665245693, -0.1848498218018628, -0.14896532674260493,
      -0.11562308824098454, -0.085307148521815951, -0.058222644279305805,
      -0.034396850607721977, -0.013745113776794951, 0.003886331100246949,
      0.018692484811857809, 0.030891279831751006, 0.040711399949960911,
      0.048384045534549017, 0.054137280103662984, 0.058192105618514262,
      0.060759731063656029, 0.062039693515435729, 0.062218609678418142,
      0.061469408319385722, 0.059950938475281594, 0.057807875980971431,
      0.055170868536958372, 0.052156871281762981, 0.048869633045788459,
      0.045400299565761894, 0.041828104744963839, 0.0382211250329898,
      0.034637075436858511, 0.031124128718412986, 0.027721742063934213,
      0.024461477982962526, 0.021367808430804026, 0.018458893170828879,
      0.01574732521127949, 0.013240837778919479, 0.010942968734440428,
      0.0088536796052147478, 0.0069699275187056512, 0.0052861892739733555,
      0.0037949376019703216, 0.0024870703470625628, 0.0013522938649174804,
      0.0003794623853893726, -0.00044312455493952579, -0.0011274642373042716,
      -0.0016856298778388182, -0.0021295840488528947, -0.0024710236509044594,
      -0.0027212536672985592, -0.00289108695995347, -0.0029907674280874323,
      -0.0030299139421043645, -0.003017482578752471, -0.0029617448144070744,
      -0.0028702794768879496, -0.0027499764078650129, -0.002607049944466464,
      -0.0024470604868432988, -0.0022749425751742488, -0.0020950380531787756,
      -0.0019111330430095982, -0.0017264975971977706, -0.0015439270250515557,
      -0.0013657840122814676, -0.0011940407612980734, -0.0010303204731883979,
      -0.00087593756676595874, -0.00073193607908741563, -0.00059912570567203236,
      -0.00047811490189565904, -0.00036934035578977954, -0.00027309192104124066,
      -0.00018953171657720415, -0.00011870548930284935, -6.0543425875285191E-5,
      -1.484634122899948E-5, 1.8748360750467846E-5, 4.0827943518984636E-5,
      5.2299743957345837E-5, 5.4522855944650814E-5, 4.9460908085498348E-5,
      3.9793661023502818E-5, 2.8811102451277326E-5, 1.9654628380546675E-5, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    day3_helicopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    day3_helicopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    day3_helicopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "day_3_test10.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(day3_helicopter_M,
                        "Error creating .mat file day_3_test10.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 6 + 1, 0, "state")) {
      rtmSetErrorStatus(day3_helicopter_M,
                        "Error writing mat file header to file day_3_test10.mat");
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
      char fileName[509] = "day_3_test10.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error closing MAT-file day_3_test10.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error reopening MAT-file day_3_test10.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 6 + 1,
           day3_helicopter_DW.ToFile3_IWORK.Count, "state")) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error writing header for state to MAT-file day_3_test10.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day3_helicopter_M,
                          "Error closing MAT-file day_3_test10.mat");
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
  day3_helicopter_M->Sizes.checksums[0] = (3222640785U);
  day3_helicopter_M->Sizes.checksums[1] = (2312860431U);
  day3_helicopter_M->Sizes.checksums[2] = (3541959234U);
  day3_helicopter_M->Sizes.checksums[3] = (685644624U);

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
