/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findNearestLocation_terminate.c
 *
 * Code generation for function 'findNearestLocation_terminate'
 *
 */

/* Include files */
#include "findNearestLocation_terminate.h"
#include "_coder_findNearestLocation_mex.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void findNearestLocation_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void findNearestLocation_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (findNearestLocation_terminate.c) */
