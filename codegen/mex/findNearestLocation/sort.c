/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "eml_int_forloop_overflow_check.h"
#include "findNearestLocation.h"
#include "findNearestLocation_data.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo i_emlrtRSI = { 72,  /* lineNo */
  "sort",                              /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sort.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 105, /* lineNo */
  "sortIdx",                           /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 308, /* lineNo */
  "block_merge_sort",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 316, /* lineNo */
  "block_merge_sort",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 317, /* lineNo */
  "block_merge_sort",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 325, /* lineNo */
  "block_merge_sort",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 333, /* lineNo */
  "block_merge_sort",                  /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 420, /* lineNo */
  "initialize_vector_sort",            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 427, /* lineNo */
  "initialize_vector_sort",            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/eml/+coder/+internal/sortIdx.m"/* pathName */
};

/* Function Definitions */
void sort(const emlrtStack *sp, real_T x[1043], int32_T idx[1043])
{
  real_T x4[4];
  int16_T idx4[4];
  real_T xwork[1043];
  int32_T nNaNs;
  int32_T ib;
  int32_T k;
  int8_T perm[4];
  int32_T i2;
  int32_T quartetOffset;
  int32_T i3;
  int32_T iwork[1043];
  int32_T i4;
  real_T d;
  real_T d1;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &i_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  memset(&idx[0], 0, 1043U * sizeof(int32_T));
  b_st.site = &j_emlrtRSI;
  c_st.site = &k_emlrtRSI;
  x4[0] = 0.0;
  idx4[0] = 0;
  x4[1] = 0.0;
  idx4[1] = 0;
  x4[2] = 0.0;
  idx4[2] = 0;
  x4[3] = 0.0;
  idx4[3] = 0;
  memset(&xwork[0], 0, 1043U * sizeof(real_T));
  nNaNs = 0;
  ib = 0;
  for (k = 0; k < 1043; k++) {
    if (muDoubleScalarIsNaN(x[k])) {
      idx[1042 - nNaNs] = k + 1;
      xwork[1042 - nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (int16_T)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        quartetOffset = k - nNaNs;
        if (x4[0] <= x4[1]) {
          ib = 1;
          i2 = 2;
        } else {
          ib = 2;
          i2 = 1;
        }

        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }

        d = x4[ib - 1];
        d1 = x4[i3 - 1];
        if (d <= d1) {
          d = x4[i2 - 1];
          if (d <= d1) {
            perm[0] = (int8_T)ib;
            perm[1] = (int8_T)i2;
            perm[2] = (int8_T)i3;
            perm[3] = (int8_T)i4;
          } else if (d <= x4[i4 - 1]) {
            perm[0] = (int8_T)ib;
            perm[1] = (int8_T)i3;
            perm[2] = (int8_T)i2;
            perm[3] = (int8_T)i4;
          } else {
            perm[0] = (int8_T)ib;
            perm[1] = (int8_T)i3;
            perm[2] = (int8_T)i4;
            perm[3] = (int8_T)i2;
          }
        } else {
          d1 = x4[i4 - 1];
          if (d <= d1) {
            if (x4[i2 - 1] <= d1) {
              perm[0] = (int8_T)i3;
              perm[1] = (int8_T)ib;
              perm[2] = (int8_T)i2;
              perm[3] = (int8_T)i4;
            } else {
              perm[0] = (int8_T)i3;
              perm[1] = (int8_T)ib;
              perm[2] = (int8_T)i4;
              perm[3] = (int8_T)i2;
            }
          } else {
            perm[0] = (int8_T)i3;
            perm[1] = (int8_T)i4;
            perm[2] = (int8_T)ib;
            perm[3] = (int8_T)i2;
          }
        }

        i3 = perm[0] - 1;
        idx[quartetOffset - 3] = idx4[i3];
        i4 = perm[1] - 1;
        idx[quartetOffset - 2] = idx4[i4];
        ib = perm[2] - 1;
        idx[quartetOffset - 1] = idx4[ib];
        i2 = perm[3] - 1;
        idx[quartetOffset] = idx4[i2];
        x[quartetOffset - 3] = x4[i3];
        x[quartetOffset - 2] = x4[i4];
        x[quartetOffset - 1] = x4[ib];
        x[quartetOffset] = x4[i2];
        ib = 0;
      }
    }
  }

  if (ib > 0) {
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    d_st.site = &p_emlrtRSI;
    if (ib > 2147483646) {
      e_st.site = &r_emlrtRSI;
      check_forloop_overflow_error(&e_st, true);
    }

    for (k = 0; k < ib; k++) {
      i3 = perm[k] - 1;
      i4 = ((k - nNaNs) - ib) + 1043;
      idx[i4] = idx4[i3];
      x[i4] = x4[i3];
    }
  }

  i2 = (nNaNs >> 1) + 1043;
  d_st.site = &q_emlrtRSI;
  for (k = 0; k <= i2 - 1044; k++) {
    ib = (k - nNaNs) + 1043;
    i3 = idx[ib];
    idx[ib] = idx[1042 - k];
    idx[1042 - k] = i3;
    x[ib] = xwork[1042 - k];
    x[1042 - k] = xwork[ib];
  }

  if ((nNaNs & 1) != 0) {
    ib = i2 - nNaNs;
    x[ib] = xwork[ib];
  }

  memset(&iwork[0], 0, 1043U * sizeof(int32_T));
  ib = 2;
  if (1043 - nNaNs > 1) {
    i2 = (1043 - nNaNs) >> 8;
    if (i2 > 0) {
      c_st.site = &l_emlrtRSI;
      for (ib = 0; ib < i2; ib++) {
        c_st.site = &m_emlrtRSI;
        merge_pow2_block(idx, x, ib << 8);
      }

      ib = i2 << 8;
      i2 = 1043 - (nNaNs + ib);
      if (i2 > 0) {
        c_st.site = &n_emlrtRSI;
        merge_block(&c_st, idx, x, ib, i2, 2, iwork, xwork);
      }

      ib = 8;
    }

    c_st.site = &o_emlrtRSI;
    merge_block(&c_st, idx, x, 0, 1043 - nNaNs, ib, iwork, xwork);
  }
}

/* End of code generation (sort.c) */
