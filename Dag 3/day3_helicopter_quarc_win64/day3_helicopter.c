/*
 * day3_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day3_helicopter".
 *
 * Model version              : 11.9
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Tue Mar  8 10:34:10 2022
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
    day3_helicopter_B.Sum3 = day3_helicopter_P.Constant_Value +
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
  rtb_Sum4[0] = day3_helicopter_P.Gain1_Gain * day3_helicopter_B.Sum3 -
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
    rtb_TmpSignalConversionAtToFile[0] = day3_helicopter_B.Sum3;
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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877561042,
      0.52359877560547063, 0.52359877560471013, 0.52359877560581414,
      0.52359877560573531, 0.52359877560107859, 0.52359877560469414,
      0.48122184543729374, 0.34333349193653251, 0.22276859766700974,
      0.11846801006246233, 0.029313828810929032, -0.045848791792632126,
      -0.1081908079852128, -0.15888039888129651, -0.19906656880133505,
      -0.22986542642481195, -0.25234909618885326, -0.26753708427422518,
      -0.27638985687884826, -0.27980436689718269, -0.27861126809970038,
      -0.27357357127134585, -0.26538651725678614, -0.25467846355418589,
      -0.24201260213344922, -0.227889345821299, -0.21274923875198115,
      -0.19697626311798616, -0.18090142991591757, -0.16480655569944946,
      -0.14892814061126292, -0.13346127522635215, -0.11856351503063356,
      -0.1043586716984336, -0.090940479737931734, -0.078376105563411291,
      -0.06670947365317198, -0.055964391194750829, -0.046147458544156639,
      -0.03725075797820615, -0.029254317648096217, -0.022128351400227064,
      -0.015835278270757813, -0.010331528037623605, -0.0055691412814735664,
      -0.00149717401769478, 0.0019370818340090912, 0.0047870500297495511,
      0.0071059809846152344, 0.0089462092477248234, 0.010358544523221824,
      0.011391783449440585, 0.012092329666047119, 0.012503910145538866,
      0.012667376313669654, 0.012620579108762797, 0.01239830781088469,
      0.012032283189604143, 0.011551196257185592, 0.010980784658550791,
      0.010343939468385965, 0.0096608358895212687, 0.00894908204718392,
      0.0082238807445529361, 0.0074981996813354712, 0.0067829462352688719,
      0.0060871434640878963, 0.0054181045011694184, 0.0047816029912023383,
      0.0041820376429974049, 0.003622589365704032, 0.0031053698034743116,
      0.0026315603935511689, 0.0022015413456002442, 0.0018150101776304162,
      0.0014710896476621826, 0.0011684250918045924, 0.00090527131962270246,
      0.00067956932755020016, 0.00048901317197469485, 0.00033110739858843719,
      0.00020321546094015375, 0.00010259959439407051, 2.6452671814225504E-5,
      -2.8077293332318121E-5, -6.3869028184782017E-5, -8.3814709803897713E-5,
      -9.0808915345075292E-5, -8.7735583755188529E-5, -7.74463979563711E-5,
      -6.27210133071765E-5, -4.6196526592434672E-5, -3.0251987160445637E-5,
      -1.6836647033890983E-5, -7.2434850605374024E-6, -1.86109448274685E-6,
      3.7999603463845233E-12, 3.7998493240820608E-12, 0.0, 0.0, 0.0, 0.0, 0.0,
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
      3.1415926535897931, 3.137842141362412, 3.1258780073570924,
      3.1018612272560429, 3.0629751326258323, 3.0073747050680981,
      2.9339699150527871, 2.8421877971597125, 2.7320779029919477,
      2.6049117753631426, 2.4630467016192492, 2.3095358731149775,
      2.1477259592808351, 1.9809395052877874, 1.812263324921493,
      1.6444297146519733, 1.4797663844211129, 1.3201914933962229,
      1.1672350742513766, 1.0220737792065515, 0.88557069855192894,
      0.75831553170414268, 0.64066271192220137, 0.53276647298915214,
      0.43461259458804885, 0.34604691652105862, 0.26680084468582432,
      0.196514095789355, 0.13475490805718077, 0.081037915123034973,
      0.034839854269406029, -0.0043867376294654857, -0.037201697381760185,
      -0.0641671399065157, -0.085838957932654322, -0.10275978229450702,
      -0.11545326963647744, -0.12441958675632359, -0.13013196441952327,
      -0.13303419826525781, -0.13353898023347474, -0.13202695057205469,
      -0.12884636769763078, -0.12431330077350375, -0.1187122576485454,
      -0.11229716860780398, -0.10529265410150442, -0.097895512129479542,
      -0.090276368190068276, -0.08258143759318004, -0.074934356437760752,
      -0.067438043636643477, -0.060176562015317306, -0.053216951702356761,
      -0.04661101376972366, -0.040397026373306617, -0.034601379498019627,
      -0.029240117843601691, -0.0243203844153154, -0.019841760028464786,
      -0.015797496221690757, -0.012175641025815051, -0.0089600586778734777,
      -0.006131345730571839, -0.0036676471108960549, -0.0015453765538244334,
      0.00026015349755180249, 0.0017741869521271387, 0.0030220581963505823,
      0.0040287992089802927, 0.0048188135270125715, 0.0054156109648504413,
      0.0058415971138090972, 0.00611791184362205, 0.0062643112781528175,
      0.0062990880111696507, 0.0062390246543470177, 0.0060993761583404806,
      0.0058938767106339052, 0.00563476738371521, 0.0053328410774984647,
      0.0049975016662815808, 0.0046368346185386294, 0.0042576867047402931,
      0.0038657527418581722, 0.0034656676413418025, 0.00306110232760397,
      0.0026548623723836105, 0.00224898843815372, 0.0018448578255006686,
      0.0014432865479789801, 0.0010446313703061439, 0.0006488910799724505,
      0.00025580583985561779, -0.00013504729025996733, -0.000524165716580862,
      -0.00091204406380456676, -0.001299114257183693, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048905487782,
      -0.04785653602721511, -0.096067120400440659, -0.15554437852606637,
      -0.22240171023163571, -0.29361916006529665, -0.36712847157388667,
      -0.44043957667806394, -0.50866451051092165, -0.56746029497496431,
      -0.6140433140144177, -0.64723965533830075, -0.66714581597468681,
      -0.67470472146731753, -0.6713344410765486, -0.65865332092300033,
      -0.63829956410096222, -0.61182567657878051, -0.58064518017952926,
      -0.54601232261884536, -0.509020667390255, -0.47061127912791756,
      -0.43158495573221356, -0.3926155136045118, -0.35426271226771194,
      -0.31698428734091466, -0.28114699558580497, -0.24703675092864152,
      -0.21486797173655722, -0.18479224341444372, -0.15690636759540774,
      -0.13125983900926869, -0.10786177009894725, -0.08668727210429153,
      -0.0676832974473205, -0.050773949367900188, -0.035865268479658846,
      -0.022849510653014965, -0.011608935382902877, -0.0020191278729877358,
      0.00604811864572314, 0.012722331497504763, 0.018132267696517837,
      0.022404172499752895, 0.025660356163199026, 0.028018058025156761,
      0.029588567887953943, 0.030476575757727848, 0.030779722387406928,
      0.03058832462169284, 0.029985251204409084, 0.029045926485254184,
      0.027838441251817263, 0.0264237517304644, 0.024855949585653212,
      0.023182587501140436, 0.021445046617672704, 0.019678933713127435,
      0.017914497547420458, 0.016177055227091196, 0.01448742078349745,
      0.012862329391763681, 0.011314851789203486, 0.00985479447869437,
      0.0084890822282952, 0.0072221202055149248, 0.0060561338183018395,
      0.0049914849768959869, 0.0040269640505112957, 0.0031600572721271277,
      0.0023871897513570608, 0.0017039445958311941, 0.0011052589192478582,
      0.00058559773813013328, 0.000139106932070144, -0.00024025342729629433,
      -0.000558593984017433, -0.00082199779082783864, -0.0010364373076693432,
      -0.0012077052248703919, -0.0013413576448650264, -0.0014426681909767095,
      -0.001516591655188906, -0.001567735851531095, -0.0016003404020623858,
      -0.001618261254949087, -0.0016249598208807137, -0.0016234957369179086,
      -0.0016165224506117662, -0.0016062851100886442, -0.0015946207106898723,
      -0.0015829611613336788, -0.0015723409604671012, -0.0015634125204625113,
      -0.0015564737052840231, -0.0015515133888949057, -0.0015482807735164715,
      -0.0015463808859864567, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875206426876, 0.23220296701329879, 0.34073399762146195,
      0.42036254467709888, 0.47252208615916713, 0.50333773580237928,
      0.51953573884722715, 0.51813489162690785, 0.48218777565591753,
      0.4155461490943585, 0.3292309875151338, 0.23461906195817528,
      0.14068914079233141, 0.053423457107580452, -0.02381985457635849,
      -0.089625314042699966, -0.14385257964063708, -0.1871073259888914,
      -0.22037183998332943, -0.24477181014742377, -0.26144289117869213,
      -0.27146288680993214, -0.27582314869862345, -0.27542113362297893,
      -0.27106295201686265, -0.26346966988506615, -0.253284291025404,
      -0.24107818173762252, -0.22735664532172545, -0.21256376115878312,
      -0.19708672001931887, -0.1812598690419438, -0.16536861478570375,
      -0.14965326478978785, -0.13431283481529688, -0.11950881414724052,
      -0.1053688625468483, -0.091990405292659283, -0.079444092965533675,
      -0.067777096906870227, -0.057016217320847662, -0.047170787399539338,
      -0.038235362865591371, -0.03019219160556802, -0.023013462513581273,
      -0.016663336295467124, -0.011099763893616887, -0.006276100466039658,
      -0.0021425245953324845, 0.0013527263138474188, 0.0042622925897988218,
      0.0066387883696344963, 0.0085340444668654536, 0.0099984852378720568,
      0.01108062678427868, 0.01182668412095822, 0.0122802753606287,
      0.012482211493367368, 0.012470360943959591, 0.012279578754244458,
      0.011941690940403205, 0.011485525300387156, 0.010936980680980701,
      0.010319127444520015, 0.0096523325924672276, 0.0089544036985371722,
      0.0082407464713558642, 0.0075245314007171826, 0.0068168655383356036,
      0.00612696602155538, 0.0054623324645644011, 0.0048289158159271484,
      0.00423128171381526, 0.0036727667613852821, 0.003155626495861541,
      0.002681174137202591, 0.0022499094772934214, 0.001861637510220282,
      0.0015155766095507062, 0.001210456231431345, 0.00094460426347919846,
      0.00071602425091854549, 0.00052246281601853184, 0.00036146765044153817,
      0.00023043651478205085, 0.00012665774789266671, 4.73429097329614E-5,
      -1.034758220108678E-5, -4.9284520246106744E-5, -7.2353608789432755E-5,
      -8.2439515848964362E-5, -8.240523760505436E-5, -7.5059519414222287E-5,
      -6.3102799742975613E-5, -4.9040891348384186E-5, -3.505761753030967E-5,
      -2.284688663956036E-5, -1.3427679451005758E-5, -7.0115347711174891E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500824421905,
      0.50469685980690948, 0.4341241224332264, 0.31851418822217981,
      0.20863816593004453, 0.12326259857366027, 0.064792012176322411,
      -0.0056033888835972509, -0.14378846388730179, -0.266566506248442,
      -0.3452606463191169, -0.37844770223187763, -0.37571968466485478,
      -0.34906273474136662, -0.30897324673783111, -0.26322183786663383,
      -0.21690906239270646, -0.17301898539357435, -0.1330580559783259,
      -0.097599880656836668, -0.066684324125170669, -0.040079982525095738,
      -0.017441047554808787, 0.0016080603027335784, 0.017432726424670131,
      0.030373128527408408, 0.040741515438858587, 0.048824437151478009,
      0.054886145663975809, 0.059171536652067945, 0.06190816455819844,
      0.063307403909854024, 0.063565017025478412, 0.062861399984335431,
      0.061361719897933, 0.059216082672909252, 0.056559806402129549,
      0.05351382901680185, 0.050185249308554863, 0.046667984235185749,
      0.043043518344636912, 0.039381719685297807, 0.035741698135996577,
      0.032172685040267386, 0.028714916368435244, 0.025400504872669691,
      0.022254289607516339, 0.019294653710392965, 0.016534303483009681,
      0.013981003636875154, 0.011638265103945609, 0.0095059831193863228,
      0.0075810243890300274, 0.0058577630840533851, 0.0043285661857309834,
      0.0029842293467252597, 0.0018143649587217239, 0.00080774453101545524,
      -4.7402197623032052E-5, -0.0007631287588498165, -0.0013515512553609694,
      -0.0018246625600669766, -0.0021941784776332712, -0.0024714129458541572,
      -0.0026671794082339831, -0.0027917155757620593, -0.0028546289087472338,
      -0.0028648602825886251, -0.0028306634495683905, -0.0027595980671569127,
      -0.0026585342280197092, -0.0025336665945982521, -0.0023905364084945515,
      -0.0022340598097825661, -0.002068561062161676, -0.0018978094347056327,
      -0.0017250586396874545, -0.0015530878683715592, -0.0013842436027664241,
      -0.0012204815125553298, -0.0010634078718877373, -0.00091432005033229707,
      -0.000774245739693597, -0.00064398066240278914, -0.00052412454274172651,
      -0.00041511506765874224, -0.00031725935274001717, -0.0002307619678370125,
      -0.00015574775228428261, -9.2276354274134221E-5, -4.0343628335531713E-5,
      1.3711288325885734E-7, 2.9382872679271882E-5, 4.7826878614005335E-5,
      5.6247633526158839E-5, 5.5933095247342343E-5, 4.8842923573182116E-5,
      3.7676828807836912E-5, 2.5664578817597564E-5, 1.5752690330658198E-5, 0.0,
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

  rtmSetTFinal(day3_helicopter_M, 36.0);
  day3_helicopter_M->Timing.stepSize0 = 0.002;
  day3_helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  day3_helicopter_M->Sizes.checksums[0] = (2387751585U);
  day3_helicopter_M->Sizes.checksums[1] = (1885723197U);
  day3_helicopter_M->Sizes.checksums[2] = (1877467796U);
  day3_helicopter_M->Sizes.checksums[3] = (1442077404U);

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
    day3_helicopter_B.Sum3 = 0.0;
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
