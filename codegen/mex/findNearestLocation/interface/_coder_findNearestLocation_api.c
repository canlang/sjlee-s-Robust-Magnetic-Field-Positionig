/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_findNearestLocation_api.c
 *
 * Code generation for function '_coder_findNearestLocation_api'
 *
 */

/* Include files */
#include "_coder_findNearestLocation_api.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2086];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *locs_2,
  const char_T *identifier))[4000];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4000];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2086];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *locs_1,
  const char_T *identifier))[2086];
static const mxArray *emlrt_marshallOut(const real_T u[2000]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4000];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2086]
{
  real_T (*y)[2086];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *locs_2,
  const char_T *identifier))[4000]
{
  real_T (*y)[4000];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(locs_2), &thisId);
  emlrtDestroyArray(&locs_2);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4000]
{
  real_T (*y)[4000];
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2086]
{
  real_T (*ret)[2086];
  static const int32_T dims[2] = { 1043, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[2086])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *locs_1,
  const char_T *identifier))[2086]
{
  real_T (*y)[2086];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(locs_1), &thisId);
  emlrtDestroyArray(&locs_1);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u[2000])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 1, 2000 };

  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4000]
{
  real_T (*ret)[4000];
  static const int32_T dims[2] = { 2000, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[4000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  void findNearestLocation_api(const mxArray * const prhs[2], int32_T nlhs,
  const mxArray *plhs[2])
{
  real_T (*distance)[2000];
  real_T (*b_index)[2000];
  real_T (*locs_1)[2086];
  real_T (*locs_2)[4000];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  distance = (real_T (*)[2000])mxMalloc(sizeof(real_T [2000]));
  b_index = (real_T (*)[2000])mxMalloc(sizeof(real_T [2000]));

  /* Marshall function inputs */
  locs_1 = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "locs_1");
  locs_2 = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "locs_2");

  /* Invoke the target function */
  findNearestLocation(&st, *locs_1, *locs_2, *distance, *b_index);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*distance);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(*b_index);
  }
}

/* End of code generation (_coder_findNearestLocation_api.c) */
