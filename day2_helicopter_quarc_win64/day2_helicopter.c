/*
 * day2_helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "day2_helicopter".
 *
 * Model version              : 11.9
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Sun Apr  3 20:40:05 2022
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "day2_helicopter.h"
#include "day2_helicopter_private.h"
#include "day2_helicopter_dt.h"

/* Block signals (default storage) */
B_day2_helicopter_T day2_helicopter_B;

/* Continuous states */
X_day2_helicopter_T day2_helicopter_X;

/* Block states (default storage) */
DW_day2_helicopter_T day2_helicopter_DW;

/* Real-time model */
static RT_MODEL_day2_helicopter_T day2_helicopter_M_;
RT_MODEL_day2_helicopter_T *const day2_helicopter_M = &day2_helicopter_M_;

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
  day2_helicopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void day2_helicopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_Clock;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[4];
  real_T rtb_Backgain;
  real_T rtb_Clock_tmp;
  real_T rtb_Derivative;
  real_T rtb_Frontgain;
  real_T *lastU;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    /* set solver stop time */
    if (!(day2_helicopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&day2_helicopter_M->solverInfo,
                            ((day2_helicopter_M->Timing.clockTickH0 + 1) *
        day2_helicopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&day2_helicopter_M->solverInfo,
                            ((day2_helicopter_M->Timing.clockTick0 + 1) *
        day2_helicopter_M->Timing.stepSize0 +
        day2_helicopter_M->Timing.clockTickH0 *
        day2_helicopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(day2_helicopter_M)) {
    day2_helicopter_M->Timing.t[0] = rtsiGetT(&day2_helicopter_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(day2_helicopter_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: day2_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (day2_helicopter_DW.HILReadEncoderTimebase_Task, 1,
         &day2_helicopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          day2_helicopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          day2_helicopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          day2_helicopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    day2_helicopter_B.TravelCounttorad = day2_helicopter_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * day2_helicopter_P.TravelCounttorad_Gain;

    /* Gain: '<S12>/Gain' */
    day2_helicopter_B.Gain = day2_helicopter_P.Gain_Gain *
      day2_helicopter_B.TravelCounttorad;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/Constant'
     */
    day2_helicopter_B.Travel = day2_helicopter_P.Constant_Value +
      day2_helicopter_B.Gain;

    /* Gain: '<S4>/Pitch: Count to rad' */
    day2_helicopter_B.PitchCounttorad = day2_helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    day2_helicopter_B.Gain_i = day2_helicopter_P.Gain_Gain_a *
      day2_helicopter_B.PitchCounttorad;
  }

  /* Clock: '<S3>/Clock' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  rtb_Clock = 0.0;
  rtb_Clock += day2_helicopter_P.TravelTransferFcn_C *
    day2_helicopter_X.TravelTransferFcn_CSTATE;
  rtb_Clock += day2_helicopter_P.TravelTransferFcn_D *
    day2_helicopter_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  day2_helicopter_B.Gain_d = day2_helicopter_P.Gain_Gain_l * rtb_Clock;

  /* Clock: '<S3>/Clock' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  rtb_Clock = 0.0;
  rtb_Clock += day2_helicopter_P.PitchTransferFcn_C *
    day2_helicopter_X.PitchTransferFcn_CSTATE;
  rtb_Clock += day2_helicopter_P.PitchTransferFcn_D *
    day2_helicopter_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  day2_helicopter_B.Gain_b = day2_helicopter_P.Gain_Gain_ae * rtb_Clock;
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    /* SignalConversion generated from: '<Root>/To File3' */
    rtb_TmpSignalConversionAtToFile[0] = day2_helicopter_B.Travel;
    rtb_TmpSignalConversionAtToFile[1] = day2_helicopter_B.Gain_d;
    rtb_TmpSignalConversionAtToFile[2] = day2_helicopter_B.Gain_i;
    rtb_TmpSignalConversionAtToFile[3] = day2_helicopter_B.Gain_b;

    /* ToFile: '<Root>/To File3' */
    {
      if (!(++day2_helicopter_DW.ToFile3_IWORK.Decimation % 1) &&
          (day2_helicopter_DW.ToFile3_IWORK.Count * (4 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) day2_helicopter_DW.ToFile3_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[4 + 1];
          day2_helicopter_DW.ToFile3_IWORK.Decimation = 0;
          u[0] = day2_helicopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          if (fwrite(u, sizeof(real_T), 4 + 1, fp) != 4 + 1) {
            rtmSetErrorStatus(day2_helicopter_M,
                              "Error writing to MAT-file day_2_r_.mat");
            return;
          }

          if (((++day2_helicopter_DW.ToFile3_IWORK.Count) * (4 + 1))+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file day_2_r_.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      day2_helicopter_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      day2_helicopter_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = day2_helicopter_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = day2_helicopter_M->Timing.t[0];

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

    day2_helicopter_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Clock = pDataValues[currTimeIndex];
        } else {
          rtb_Clock = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Clock = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    day2_helicopter_B.ElevationCounttorad = day2_helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * day2_helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S7>/Gain' */
    day2_helicopter_B.Gain_e = day2_helicopter_P.Gain_Gain_lv *
      day2_helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    day2_helicopter_B.Sum = day2_helicopter_B.Gain_e +
      day2_helicopter_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S8>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  day2_helicopter_B.Gain_dg = (day2_helicopter_P.ElevationTransferFcn_C *
    day2_helicopter_X.ElevationTransferFcn_CSTATE +
    day2_helicopter_P.ElevationTransferFcn_D *
    day2_helicopter_B.ElevationCounttorad) * day2_helicopter_P.Gain_Gain_n;

  /* Clock: '<S3>/Clock' incorporates:
   *  Gain: '<S2>/Gain1'
   *  Sum: '<S5>/Sum2'
   */
  rtb_Clock -= day2_helicopter_P.Gain1_Gain * day2_helicopter_B.Gain_i;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum3'
   */
  rtb_Frontgain = (day2_helicopter_P.K_pp * rtb_Clock -
                   day2_helicopter_P.Gain1_Gain * day2_helicopter_B.Gain_b *
                   day2_helicopter_P.K_pd) + day2_helicopter_P.Vd_ff;

  /* Integrator: '<S3>/Integrator' */
  /* Limited  Integrator  */
  if (day2_helicopter_X.Integrator_CSTATE >=
      day2_helicopter_P.Integrator_UpperSat) {
    day2_helicopter_X.Integrator_CSTATE = day2_helicopter_P.Integrator_UpperSat;
  } else {
    if (day2_helicopter_X.Integrator_CSTATE <=
        day2_helicopter_P.Integrator_LowerSat) {
      day2_helicopter_X.Integrator_CSTATE =
        day2_helicopter_P.Integrator_LowerSat;
    }
  }

  /* Clock: '<S3>/Clock' incorporates:
   *  Integrator: '<S3>/Integrator'
   */
  rtb_Clock = day2_helicopter_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   *  Gain: '<S2>/Gain1'
   */
  rtb_Derivative = day2_helicopter_P.elevation_ref_Value -
    day2_helicopter_P.Gain1_Gain * day2_helicopter_B.Sum;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  rtb_Backgain = ((day2_helicopter_P.K_ep * rtb_Derivative + rtb_Clock) -
                  day2_helicopter_P.Gain1_Gain * day2_helicopter_B.Gain_dg *
                  day2_helicopter_P.K_ed) + day2_helicopter_P.Vs_ff;

  /* Clock: '<S3>/Clock' incorporates:
   *  Derivative: '<S4>/Derivative'
   */
  rtb_Clock_tmp = day2_helicopter_M->Timing.t[0];

  /* Clock: '<S3>/Clock' */
  rtb_Clock = rtb_Clock_tmp;

  /* If: '<S3>/If' incorporates:
   *  Gain: '<S3>/K_ei'
   *  Inport: '<S6>/In1'
   */
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    rtAction = (int8_T)!(rtb_Clock >= 2.0);
    day2_helicopter_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = day2_helicopter_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S3>/If Action Subsystem' incorporates:
     *  ActionPort: '<S6>/Action Port'
     */
    day2_helicopter_B.In1 = day2_helicopter_P.K_ei * rtb_Derivative;
    if (rtmIsMajorTimeStep(day2_helicopter_M)) {
      srUpdateBC(day2_helicopter_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S3>/If Action Subsystem' */
  }

  /* End of If: '<S3>/If' */
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((day2_helicopter_DW.TimeStampA >= rtb_Clock_tmp) &&
      (day2_helicopter_DW.TimeStampB >= rtb_Clock_tmp)) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = day2_helicopter_DW.TimeStampA;
    lastU = &day2_helicopter_DW.LastUAtTimeA;
    if (day2_helicopter_DW.TimeStampA < day2_helicopter_DW.TimeStampB) {
      if (day2_helicopter_DW.TimeStampB < rtb_Clock_tmp) {
        rtb_Derivative = day2_helicopter_DW.TimeStampB;
        lastU = &day2_helicopter_DW.LastUAtTimeB;
      }
    } else {
      if (day2_helicopter_DW.TimeStampA >= rtb_Clock_tmp) {
        rtb_Derivative = day2_helicopter_DW.TimeStampB;
        lastU = &day2_helicopter_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (day2_helicopter_B.PitchCounttorad - *lastU) /
      (rtb_Clock_tmp - rtb_Derivative);
  }

  /* Gain: '<S11>/Gain' */
  day2_helicopter_B.Gain_l = day2_helicopter_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Derivative = (rtb_Backgain - rtb_Frontgain) *
    day2_helicopter_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Derivative > day2_helicopter_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day2_helicopter_B.BackmotorSaturation =
      day2_helicopter_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Derivative < day2_helicopter_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day2_helicopter_B.BackmotorSaturation =
      day2_helicopter_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Back motor: Saturation' */
    day2_helicopter_B.BackmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (rtb_Frontgain + rtb_Backgain) *
    day2_helicopter_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > day2_helicopter_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day2_helicopter_B.FrontmotorSaturation =
      day2_helicopter_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < day2_helicopter_P.FrontmotorSaturation_LowerSat) {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day2_helicopter_B.FrontmotorSaturation =
      day2_helicopter_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S4>/Front motor: Saturation' */
    day2_helicopter_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: day2_helicopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      day2_helicopter_DW.HILWriteAnalog_Buffer[0] =
        day2_helicopter_B.FrontmotorSaturation;
      day2_helicopter_DW.HILWriteAnalog_Buffer[1] =
        day2_helicopter_B.BackmotorSaturation;
      result = hil_write_analog(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILWriteAnalog_channels, 2,
        &day2_helicopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void day2_helicopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (day2_helicopter_DW.TimeStampA == (rtInf)) {
    day2_helicopter_DW.TimeStampA = day2_helicopter_M->Timing.t[0];
    lastU = &day2_helicopter_DW.LastUAtTimeA;
  } else if (day2_helicopter_DW.TimeStampB == (rtInf)) {
    day2_helicopter_DW.TimeStampB = day2_helicopter_M->Timing.t[0];
    lastU = &day2_helicopter_DW.LastUAtTimeB;
  } else if (day2_helicopter_DW.TimeStampA < day2_helicopter_DW.TimeStampB) {
    day2_helicopter_DW.TimeStampA = day2_helicopter_M->Timing.t[0];
    lastU = &day2_helicopter_DW.LastUAtTimeA;
  } else {
    day2_helicopter_DW.TimeStampB = day2_helicopter_M->Timing.t[0];
    lastU = &day2_helicopter_DW.LastUAtTimeB;
  }

  *lastU = day2_helicopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(day2_helicopter_M)) {
    rt_ertODEUpdateContinuousStates(&day2_helicopter_M->solverInfo);
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
  if (!(++day2_helicopter_M->Timing.clockTick0)) {
    ++day2_helicopter_M->Timing.clockTickH0;
  }

  day2_helicopter_M->Timing.t[0] = rtsiGetSolverStopTime
    (&day2_helicopter_M->solverInfo);

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
    if (!(++day2_helicopter_M->Timing.clockTick1)) {
      ++day2_helicopter_M->Timing.clockTickH1;
    }

    day2_helicopter_M->Timing.t[1] = day2_helicopter_M->Timing.clockTick1 *
      day2_helicopter_M->Timing.stepSize1 +
      day2_helicopter_M->Timing.clockTickH1 *
      day2_helicopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void day2_helicopter_derivatives(void)
{
  XDot_day2_helicopter_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_day2_helicopter_T *) day2_helicopter_M->derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += day2_helicopter_P.TravelTransferFcn_A *
    day2_helicopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += day2_helicopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += day2_helicopter_P.PitchTransferFcn_A *
    day2_helicopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += day2_helicopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    day2_helicopter_P.ElevationTransferFcn_A *
    day2_helicopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += day2_helicopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  lsat = (day2_helicopter_X.Integrator_CSTATE <=
          day2_helicopter_P.Integrator_LowerSat);
  usat = (day2_helicopter_X.Integrator_CSTATE >=
          day2_helicopter_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (day2_helicopter_B.In1 > 0.0)) || (usat &&
       (day2_helicopter_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = day2_helicopter_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */
}

/* Model initialize function */
void day2_helicopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day2_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &day2_helicopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(day2_helicopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(day2_helicopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
      return;
    }

    if ((day2_helicopter_P.HILInitialize_AIPStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &day2_helicopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (day2_helicopter_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &day2_helicopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = day2_helicopter_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_AIChannels, 8U,
        &day2_helicopter_DW.HILInitialize_AIMinimums[0],
        &day2_helicopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_AOPStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &day2_helicopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (day2_helicopter_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &day2_helicopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = day2_helicopter_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (day2_helicopter_DW.HILInitialize_Card,
         day2_helicopter_P.HILInitialize_AOChannels, 8U,
         &day2_helicopter_DW.HILInitialize_AOMinimums[0],
         &day2_helicopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_AOStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day2_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day2_helicopter_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_AOChannels, 8U,
        &day2_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day2_helicopter_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day2_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day2_helicopter_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (day2_helicopter_DW.HILInitialize_Card,
         day2_helicopter_P.HILInitialize_AOChannels, 8U,
         &day2_helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_EIPStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &day2_helicopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = day2_helicopter_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (day2_helicopter_DW.HILInitialize_Card,
         day2_helicopter_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &day2_helicopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_EIStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &day2_helicopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = day2_helicopter_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_EIChannels, 8U,
        &day2_helicopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_POPStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day2_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day2_helicopter_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &day2_helicopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          day2_helicopter_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &day2_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            day2_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            day2_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = day2_helicopter_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            day2_helicopter_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            day2_helicopter_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] = day2_helicopter_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(day2_helicopter_DW.HILInitialize_Card,
          &day2_helicopter_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes, &day2_helicopter_DW.HILInitialize_POSortedFreqs
          [0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(day2_helicopter_DW.HILInitialize_Card,
          &day2_helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &day2_helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &day2_helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = day2_helicopter_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &day2_helicopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = day2_helicopter_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &day2_helicopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = day2_helicopter_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &day2_helicopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &day2_helicopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &day2_helicopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &day2_helicopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = day2_helicopter_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &day2_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day2_helicopter_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_POChannels, 8U,
        &day2_helicopter_DW.HILInitialize_POSortedFreqs[0],
        &day2_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((day2_helicopter_P.HILInitialize_POStart && !is_switching) ||
        (day2_helicopter_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day2_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day2_helicopter_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(day2_helicopter_DW.HILInitialize_Card,
        day2_helicopter_P.HILInitialize_POChannels, 8U,
        &day2_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }

    if (day2_helicopter_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &day2_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day2_helicopter_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (day2_helicopter_DW.HILInitialize_Card,
         day2_helicopter_P.HILInitialize_POChannels, 8U,
         &day2_helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: day2_helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (day2_helicopter_DW.HILInitialize_Card,
       day2_helicopter_P.HILReadEncoderTimebase_SamplesI,
       day2_helicopter_P.HILReadEncoderTimebase_Channels, 3,
       &day2_helicopter_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (day2_helicopter_DW.HILReadEncoderTimebase_Task, (t_buffer_overflow_mode)
         (day2_helicopter_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
    }
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "day_2_r_.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(day2_helicopter_M,
                        "Error creating .mat file day_2_r_.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 4 + 1, 0, "state")) {
      rtmSetErrorStatus(day2_helicopter_M,
                        "Error writing mat file header to file day_2_r_.mat");
      return;
    }

    day2_helicopter_DW.ToFile3_IWORK.Count = 0;
    day2_helicopter_DW.ToFile3_IWORK.Decimation = -1;
    day2_helicopter_DW.ToFile3_PWORK.FilePtr = fp;
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

    day2_helicopter_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    day2_helicopter_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    day2_helicopter_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for If: '<S3>/If' */
  day2_helicopter_DW.If_ActiveSubsystem = -1;

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  day2_helicopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  day2_helicopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  day2_helicopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  day2_helicopter_X.Integrator_CSTATE = day2_helicopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  day2_helicopter_DW.TimeStampA = (rtInf);
  day2_helicopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void day2_helicopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: day2_helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(day2_helicopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(day2_helicopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((day2_helicopter_P.HILInitialize_AOTerminate && !is_switching) ||
        (day2_helicopter_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &day2_helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = day2_helicopter_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((day2_helicopter_P.HILInitialize_POTerminate && !is_switching) ||
        (day2_helicopter_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &day2_helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = day2_helicopter_P.HILInitialize_POFinal;
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
      result = hil_write(day2_helicopter_DW.HILInitialize_Card
                         , day2_helicopter_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , day2_helicopter_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &day2_helicopter_DW.HILInitialize_AOVoltages[0]
                         , &day2_helicopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(day2_helicopter_DW.HILInitialize_Card,
            day2_helicopter_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &day2_helicopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(day2_helicopter_DW.HILInitialize_Card,
            day2_helicopter_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &day2_helicopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(day2_helicopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(day2_helicopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(day2_helicopter_DW.HILInitialize_Card);
    hil_close(day2_helicopter_DW.HILInitialize_Card);
    day2_helicopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (FILE *) day2_helicopter_DW.ToFile3_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "day_2_r_.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day2_helicopter_M,
                          "Error closing MAT-file day_2_r_.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(day2_helicopter_M,
                          "Error reopening MAT-file day_2_r_.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 4 + 1,
           day2_helicopter_DW.ToFile3_IWORK.Count, "state")) {
        rtmSetErrorStatus(day2_helicopter_M,
                          "Error writing header for state to MAT-file day_2_r_.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(day2_helicopter_M,
                          "Error closing MAT-file day_2_r_.mat");
        return;
      }

      day2_helicopter_DW.ToFile3_PWORK.FilePtr = (NULL);
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
  day2_helicopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  day2_helicopter_update();
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
  day2_helicopter_initialize();
}

void MdlTerminate(void)
{
  day2_helicopter_terminate();
}

/* Registration function */
RT_MODEL_day2_helicopter_T *day2_helicopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  day2_helicopter_P.Integrator_UpperSat = rtInf;
  day2_helicopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)day2_helicopter_M, 0,
                sizeof(RT_MODEL_day2_helicopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&day2_helicopter_M->solverInfo,
                          &day2_helicopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&day2_helicopter_M->solverInfo, &rtmGetTPtr(day2_helicopter_M));
    rtsiSetStepSizePtr(&day2_helicopter_M->solverInfo,
                       &day2_helicopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&day2_helicopter_M->solverInfo, &day2_helicopter_M->derivs);
    rtsiSetContStatesPtr(&day2_helicopter_M->solverInfo, (real_T **)
                         &day2_helicopter_M->contStates);
    rtsiSetNumContStatesPtr(&day2_helicopter_M->solverInfo,
      &day2_helicopter_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&day2_helicopter_M->solverInfo,
      &day2_helicopter_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&day2_helicopter_M->solverInfo,
      &day2_helicopter_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&day2_helicopter_M->solverInfo,
      &day2_helicopter_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&day2_helicopter_M->solverInfo, (&rtmGetErrorStatus
      (day2_helicopter_M)));
    rtsiSetRTModelPtr(&day2_helicopter_M->solverInfo, day2_helicopter_M);
  }

  rtsiSetSimTimeStep(&day2_helicopter_M->solverInfo, MAJOR_TIME_STEP);
  day2_helicopter_M->intgData.f[0] = day2_helicopter_M->odeF[0];
  day2_helicopter_M->contStates = ((real_T *) &day2_helicopter_X);
  rtsiSetSolverData(&day2_helicopter_M->solverInfo, (void *)
                    &day2_helicopter_M->intgData);
  rtsiSetSolverName(&day2_helicopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = day2_helicopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    day2_helicopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    day2_helicopter_M->Timing.sampleTimes =
      (&day2_helicopter_M->Timing.sampleTimesArray[0]);
    day2_helicopter_M->Timing.offsetTimes =
      (&day2_helicopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    day2_helicopter_M->Timing.sampleTimes[0] = (0.0);
    day2_helicopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    day2_helicopter_M->Timing.offsetTimes[0] = (0.0);
    day2_helicopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(day2_helicopter_M, &day2_helicopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = day2_helicopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    day2_helicopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(day2_helicopter_M, 25.0);
  day2_helicopter_M->Timing.stepSize0 = 0.002;
  day2_helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  day2_helicopter_M->Sizes.checksums[0] = (2750373368U);
  day2_helicopter_M->Sizes.checksums[1] = (1394997924U);
  day2_helicopter_M->Sizes.checksums[2] = (1010672269U);
  day2_helicopter_M->Sizes.checksums[3] = (4127653228U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    day2_helicopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)
      &day2_helicopter_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(day2_helicopter_M->extModeInfo,
      &day2_helicopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(day2_helicopter_M->extModeInfo,
                        day2_helicopter_M->Sizes.checksums);
    rteiSetTPtr(day2_helicopter_M->extModeInfo, rtmGetTPtr(day2_helicopter_M));
  }

  day2_helicopter_M->solverInfoPtr = (&day2_helicopter_M->solverInfo);
  day2_helicopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&day2_helicopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&day2_helicopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  day2_helicopter_M->blockIO = ((void *) &day2_helicopter_B);

  {
    day2_helicopter_B.TravelCounttorad = 0.0;
    day2_helicopter_B.Gain = 0.0;
    day2_helicopter_B.Travel = 0.0;
    day2_helicopter_B.Gain_d = 0.0;
    day2_helicopter_B.PitchCounttorad = 0.0;
    day2_helicopter_B.Gain_i = 0.0;
    day2_helicopter_B.Gain_b = 0.0;
    day2_helicopter_B.ElevationCounttorad = 0.0;
    day2_helicopter_B.Gain_e = 0.0;
    day2_helicopter_B.Sum = 0.0;
    day2_helicopter_B.Gain_dg = 0.0;
    day2_helicopter_B.Gain_l = 0.0;
    day2_helicopter_B.BackmotorSaturation = 0.0;
    day2_helicopter_B.FrontmotorSaturation = 0.0;
    day2_helicopter_B.In1 = 0.0;
  }

  /* parameters */
  day2_helicopter_M->defaultParam = ((real_T *)&day2_helicopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &day2_helicopter_X;
    day2_helicopter_M->contStates = (x);
    (void) memset((void *)&day2_helicopter_X, 0,
                  sizeof(X_day2_helicopter_T));
  }

  /* states (dwork) */
  day2_helicopter_M->dwork = ((void *) &day2_helicopter_DW);
  (void) memset((void *)&day2_helicopter_DW, 0,
                sizeof(DW_day2_helicopter_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      day2_helicopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  day2_helicopter_DW.TimeStampA = 0.0;
  day2_helicopter_DW.LastUAtTimeA = 0.0;
  day2_helicopter_DW.TimeStampB = 0.0;
  day2_helicopter_DW.LastUAtTimeB = 0.0;
  day2_helicopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  day2_helicopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    day2_helicopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  day2_helicopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  day2_helicopter_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  day2_helicopter_M->Sizes.numY = (0); /* Number of model outputs */
  day2_helicopter_M->Sizes.numU = (0); /* Number of model inputs */
  day2_helicopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  day2_helicopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  day2_helicopter_M->Sizes.numBlocks = (62);/* Number of blocks */
  day2_helicopter_M->Sizes.numBlockIO = (15);/* Number of block outputs */
  day2_helicopter_M->Sizes.numBlockPrms = (145);/* Sum of parameter "widths" */
  return day2_helicopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
