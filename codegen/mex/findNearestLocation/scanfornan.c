/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * scanfornan.c
 *
 * Code generation for function 'scanfornan'
 *
 */

/* Include files */
#include "scanfornan.h"
#include "findNearestLocation.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_scanfornan(const emlrtStack *sp, const real_T X[4000], boolean_T nanobs
                  [2000])
{
  int32_T i;
  int32_T qq;
  boolean_T nanflag;
  int32_T jj;
  boolean_T exitg1;
  jmp_buf * volatile emlrtJBStack;
  for (i = 0; i < 2000; i++) {
    nanobs[i] = true;
  }

  emlrtEnterParallelRegion(sp, omp_in_parallel());
  emlrtPushJmpBuf(sp, &emlrtJBStack);

#pragma omp parallel for \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(nanflag,jj,exitg1)

  for (qq = 0; qq < 2000; qq++) {
    nanflag = false;
    jj = 0;
    exitg1 = false;
    while ((!exitg1) && (jj < 2)) {
      if (muDoubleScalarIsNaN(X[jj + (qq << 1)])) {
        nanflag = true;
        exitg1 = true;
      } else {
        jj++;
      }
    }

    if (nanflag) {
      nanobs[qq] = false;
    }
  }

  emlrtPopJmpBuf(sp, &emlrtJBStack);
  emlrtExitParallelRegion(sp, omp_in_parallel());
}

void scanfornan(const emlrtStack *sp, const real_T X[2086], boolean_T nanobs
                [1043])
{
  int32_T i;
  int32_T qq;
  boolean_T nanflag;
  int32_T jj;
  boolean_T exitg1;
  jmp_buf * volatile emlrtJBStack;
  for (i = 0; i < 1043; i++) {
    nanobs[i] = true;
  }

  emlrtEnterParallelRegion(sp, omp_in_parallel());
  emlrtPushJmpBuf(sp, &emlrtJBStack);

#pragma omp parallel for \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(nanflag,jj,exitg1)

  for (qq = 0; qq < 1043; qq++) {
    nanflag = false;
    jj = 0;
    exitg1 = false;
    while ((!exitg1) && (jj < 2)) {
      if (muDoubleScalarIsNaN(X[jj + (qq << 1)])) {
        nanflag = true;
        exitg1 = true;
      } else {
        jj++;
      }
    }

    if (nanflag) {
      nanobs[qq] = false;
    }
  }

  emlrtPopJmpBuf(sp, &emlrtJBStack);
  emlrtExitParallelRegion(sp, omp_in_parallel());
}

/* End of code generation (scanfornan.c) */
