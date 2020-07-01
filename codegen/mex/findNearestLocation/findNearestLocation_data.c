/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findNearestLocation_data.c
 *
 * Code generation for function 'findNearestLocation_data'
 *
 */

/* Include files */
#include "findNearestLocation_data.h"
#include "findNearestLocation.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
omp_lock_t emlrtLockGlobal;
omp_nest_lock_t emlrtNestLockGlobal;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131594U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "findNearestLocation",               /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2433290357U, 2237796540U, 4066813863U, 833189415U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo r_emlrtRSI = { 21,         /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo s_emlrtRSI = { 587,        /* lineNo */
  "merge_pow2_block",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

emlrtRSInfo t_emlrtRSI = { 589,        /* lineNo */
  "merge_pow2_block",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

emlrtRSInfo u_emlrtRSI = { 617,        /* lineNo */
  "merge_pow2_block",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

emlrtRSInfo w_emlrtRSI = { 506,        /* lineNo */
  "merge_block",                       /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

/* End of code generation (findNearestLocation_data.c) */
