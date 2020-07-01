/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findNearestLocation_initialize.c
 *
 * Code generation for function 'findNearestLocation_initialize'
 *
 */

/* Include files */
#include "findNearestLocation_initialize.h"
#include "_coder_findNearestLocation_mex.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;

/* Function Definitions */
void findNearestLocation_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtLicenseCheckR2012b(&st, "Statistics_Toolbox", 2);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (findNearestLocation_initialize.c) */
