/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sortIdx.c
 *
 * Code generation for function 'sortIdx'
 *
 */

/* Include files */
#include "sortIdx.h"
#include "eml_int_forloop_overflow_check.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo v_emlrtRSI = { 499, /* lineNo */
  "merge_block",                       /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 507, /* lineNo */
  "merge_block",                       /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 514, /* lineNo */
  "merge_block",                       /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 530,/* lineNo */
  "merge",                             /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 561,/* lineNo */
  "merge",                             /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { 0,     /* iFirst */
  31,                                  /* iLast */
  490,                                 /* lineNo */
  24,                                  /* colNo */
  "",                                  /* aName */
  "merge_block",                       /* fName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m",/* pName */
  1                                    /* checkKind */
};

/* Function Declarations */
static void merge(const emlrtStack *sp, int32_T idx[1043], real_T x[1043],
                  int32_T offset, int32_T np, int32_T nq, int32_T iwork[1043],
                  real_T xwork[1043]);

/* Function Definitions */
static void merge(const emlrtStack *sp, int32_T idx[1043], real_T x[1043],
                  int32_T offset, int32_T np, int32_T nq, int32_T iwork[1043],
                  real_T xwork[1043])
{
  int32_T n_tmp;
  int32_T j;
  int32_T p;
  int32_T iout;
  int32_T q;
  int32_T exitg1;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if ((np != 0) && (nq != 0)) {
    n_tmp = np + nq;
    st.site = &ab_emlrtRSI;
    if ((1 <= n_tmp) && (n_tmp > 2147483646)) {
      b_st.site = &r_emlrtRSI;
      check_forloop_overflow_error(&b_st, true);
    }

    for (j = 0; j < n_tmp; j++) {
      iout = offset + j;
      iwork[j] = idx[iout];
      xwork[j] = x[iout];
    }

    p = 0;
    q = np;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[q]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[q];
        x[iout] = xwork[q];
        if (q + 1 < n_tmp) {
          q++;
        } else {
          q = iout - p;
          st.site = &bb_emlrtRSI;
          if ((p + 1 <= np) && (np > 2147483646)) {
            b_st.site = &r_emlrtRSI;
            check_forloop_overflow_error(&b_st, true);
          }

          for (j = p + 1; j <= np; j++) {
            iout = q + j;
            idx[iout] = iwork[j - 1];
            x[iout] = xwork[j - 1];
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

void merge_block(const emlrtStack *sp, int32_T idx[1043], real_T x[1043],
                 int32_T offset, int32_T n, int32_T preSortLevel, int32_T iwork
                 [1043], real_T xwork[1043])
{
  int32_T nPairs;
  int32_T bLen;
  int32_T tailOffset;
  int32_T nTail;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  if ((preSortLevel < 0) || (preSortLevel > 31)) {
    emlrtDynamicBoundsCheckR2012b(preSortLevel, 0, 31, &emlrtBCI, sp);
  }

  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        st.site = &v_emlrtRSI;
        merge(&st, idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      st.site = &x_emlrtRSI;
      merge(&st, idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    st.site = &y_emlrtRSI;
    merge(&st, idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

void merge_pow2_block(int32_T idx[1043], real_T x[1043], int32_T offset)
{
  int32_T b;
  int32_T bLen;
  int32_T bLen2;
  int32_T nPairs;
  int32_T k;
  int32_T blockOffset;
  int32_T j;
  int32_T p;
  int32_T iout;
  int32_T q;
  int32_T iwork[256];
  real_T xwork[256];
  int32_T exitg1;
  for (b = 0; b < 6; b++) {
    bLen = 1 << (b + 2);
    bLen2 = bLen << 1;
    nPairs = 256 >> (b + 3);
    for (k = 0; k < nPairs; k++) {
      blockOffset = offset + k * bLen2;
      for (j = 0; j < bLen2; j++) {
        iout = blockOffset + j;
        iwork[j] = idx[iout];
        xwork[j] = x[iout];
      }

      p = 0;
      q = bLen;
      iout = blockOffset - 1;
      do {
        exitg1 = 0;
        iout++;
        if (xwork[p] <= xwork[q]) {
          idx[iout] = iwork[p];
          x[iout] = xwork[p];
          if (p + 1 < bLen) {
            p++;
          } else {
            exitg1 = 1;
          }
        } else {
          idx[iout] = iwork[q];
          x[iout] = xwork[q];
          if (q + 1 < bLen2) {
            q++;
          } else {
            iout -= p;
            for (j = p + 1; j <= bLen; j++) {
              q = iout + j;
              idx[q] = iwork[j - 1];
              x[q] = xwork[j - 1];
            }

            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }
}

/* End of code generation (sortIdx.c) */
