/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_findNearestLocation_mex.c
 *
 * Code generation for function '_coder_findNearestLocation_mex'
 *
 */

/* Include files */
#include "_coder_findNearestLocation_mex.h"
#include "_coder_findNearestLocation_api.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "findNearestLocation_initialize.h"
#include "findNearestLocation_mexutil.h"
#include "findNearestLocation_terminate.h"

/* Variable Definitions */
static jmp_buf emlrtJBEnviron;

/* Function Declarations */
MEXFUNCTION_LINKAGE void findNearestLocation_mexFunction(int32_T nlhs, mxArray
  *plhs[2], int32_T nrhs, const mxArray *prhs[2]);

/* Function Definitions */
void findNearestLocation_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T
  nrhs, const mxArray *prhs[2])
{
  const mxArray *outputs[2];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        19, "findNearestLocation");
  }

  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 19,
                        "findNearestLocation");
  }

  /* Call the function. */
  findNearestLocation_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexAtExit(&findNearestLocation_atexit);

  /* Initialize the memory manager. */
  omp_init_lock(&emlrtLockGlobal);
  omp_init_nest_lock(&emlrtNestLockGlobal);

  /* Module initialization. */
  findNearestLocation_initialize();
  st.tls = emlrtRootTLSGlobal;
  emlrtSetJmpBuf(&st, &emlrtJBEnviron);
  if (setjmp(emlrtJBEnviron) == 0) {
    /* Dispatch the entry-point. */
    findNearestLocation_mexFunction(nlhs, plhs, nrhs, prhs);

    /* Module termination. */
    findNearestLocation_terminate();
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&emlrtNestLockGlobal);
  } else {
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&emlrtNestLockGlobal);
    emlrtReportParallelRunTimeError(&st);
  }
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal,
                     &emlrtLockerFunction, omp_get_num_procs());
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_findNearestLocation_mex.c) */
