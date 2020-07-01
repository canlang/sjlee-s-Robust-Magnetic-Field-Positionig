/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findNearestLocation.c
 *
 * Code generation for function 'findNearestLocation'
 *
 */

/* Include files */
#include "findNearestLocation.h"
#include "pdist2.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 2,     /* lineNo */
  "findNearestLocation",               /* fcnName */
  "/Users/sjlee/git/magnetic_field_pf_exp/findNearestLocation.m"/* pathName */
};

/* Function Definitions */
void findNearestLocation(const emlrtStack *sp, const real_T locs_1[2086], const
  real_T locs_2[4000], real_T distance[2000], real_T b_index[2000])
{
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &emlrtRSI;
  pdist2(&st, locs_1, locs_2, distance, b_index);
}

/* End of code generation (findNearestLocation.c) */
