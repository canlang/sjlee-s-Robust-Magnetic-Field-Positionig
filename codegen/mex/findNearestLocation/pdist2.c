/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pdist2.c
 *
 * Code generation for function 'pdist2'
 *
 */

/* Include files */
#include "pdist2.h"
#include "findNearestLocation.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "scanfornan.h"
#include "sort.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo b_emlrtRSI = { 361, /* lineNo */
  "pdist2",                            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/eml/pdist2.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 362, /* lineNo */
  "pdist2",                            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/eml/pdist2.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 448, /* lineNo */
  "pdist2",                            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/eml/pdist2.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 459, /* lineNo */
  "pdist2",                            /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/eml/pdist2.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 27,  /* lineNo */
  "@(temp) sqrt(temp)",                /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/stats/+stats/+coder/+distutils/disthandles.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 743, /* lineNo */
  "partialSort",                       /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/stats/eml/pdist2.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 27,  /* lineNo */
  "sort",                              /* fcnName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/lib/matlab/datafun/sort.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 13,  /* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/Applications/MATLAB_R2020a.app/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

/* Function Declarations */
static void partialSort(const emlrtStack *sp, const real_T Din[1043], real_T
  *Dout, real_T *Iout);

/* Function Definitions */
static void partialSort(const emlrtStack *sp, const real_T Din[1043], real_T
  *Dout, real_T *Iout)
{
  real_T D[1043];
  int32_T iidx[1043];
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &g_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  memcpy(&D[0], &Din[0], 1043U * sizeof(real_T));
  b_st.site = &h_emlrtRSI;
  sort(&b_st, D, iidx);
  *Dout = D[0];
  *Iout = iidx[0];
}

void pdist2(const emlrtStack *sp, const real_T Xin[2086], const real_T Yin[4000],
            real_T D[2000], real_T b_I[2000])
{
  int32_T i;
  int32_T X_tmp;
  real_T X[2086];
  boolean_T logIndX[1043];
  real_T Y[4000];
  boolean_T logIndY[2000];
  int32_T ii;
  real_T tempSum;
  real_T tempSum_tmp;
  real_T d[1043];
  int32_T b_i;
  int32_T qq;
  int32_T i1;
  jmp_buf * volatile emlrtJBStack;
  emlrtStack st;
  emlrtStack b_st;
  jmp_buf b_emlrtJBEnviron;
  emlrtStack c_st;
  emlrtStack d_st;
  boolean_T emlrtHadParallelError = false;
  st.prev = sp;
  st.tls = sp->tls;
  for (i = 0; i < 1043; i++) {
    X_tmp = i << 1;
    X[X_tmp] = Xin[i];
    X[X_tmp + 1] = Xin[i + 1043];
  }

  for (i = 0; i < 2000; i++) {
    X_tmp = i << 1;
    Y[X_tmp] = Yin[i];
    Y[X_tmp + 1] = Yin[i + 2000];
  }

  st.site = &b_emlrtRSI;
  scanfornan(&st, X, logIndX);
  st.site = &c_emlrtRSI;
  b_scanfornan(&st, Y, logIndY);
  emlrtEnterParallelRegion(sp, omp_in_parallel());
  emlrtPushJmpBuf(sp, &emlrtJBStack);

#pragma omp parallel \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(tempSum,tempSum_tmp,d,b_emlrtJBEnviron,d_st,b_i,qq,i1) \
 firstprivate(b_st,c_st,emlrtHadParallelError)

  {
    if (setjmp(b_emlrtJBEnviron) == 0) {
      b_st.prev = sp;
      b_st.tls = emlrtAllocTLS(sp, omp_get_thread_num());
      b_st.site = NULL;
      emlrtSetJmpBuf(&b_st, &b_emlrtJBEnviron);
      c_st.prev = &b_st;
      c_st.tls = b_st.tls;
      d_st.prev = &c_st;
      d_st.tls = c_st.tls;
    } else {
      emlrtHadParallelError = true;
    }

#pragma omp for nowait

    for (ii = 0; ii < 2000; ii++) {
      if (emlrtHadParallelError)
        continue;
      if (setjmp(b_emlrtJBEnviron) == 0) {
        for (b_i = 0; b_i < 1043; b_i++) {
          d[b_i] = rtNaN;
        }

        if (logIndY[ii]) {
          for (qq = 0; qq < 1043; qq++) {
            if (logIndX[qq]) {
              b_i = qq << 1;
              i1 = ii << 1;
              tempSum = X[b_i] - Y[i1];
              tempSum_tmp = X[b_i + 1] - Y[i1 + 1];
              tempSum = tempSum * tempSum + tempSum_tmp * tempSum_tmp;
              c_st.site = &d_emlrtRSI;
              d_st.site = &f_emlrtRSI;
              if (tempSum < 0.0) {
                emlrtErrorWithMessageIdR2018a(&d_st, &emlrtRTEI,
                  "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
              }

              d[qq] = muDoubleScalarSqrt(tempSum);
            }
          }
        }

        c_st.site = &e_emlrtRSI;
        partialSort(&c_st, d, &tempSum_tmp, &tempSum);
        D[ii] = tempSum_tmp;
        b_I[ii] = tempSum;
      } else {
        emlrtHadParallelError = true;
      }
    }
  }

  emlrtPopJmpBuf(sp, &emlrtJBStack);
  emlrtExitParallelRegion(sp, omp_in_parallel());
}

/* End of code generation (pdist2.c) */
