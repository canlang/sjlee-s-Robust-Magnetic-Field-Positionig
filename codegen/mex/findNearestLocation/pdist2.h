/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * pdist2.h
 *
 * Code generation for function 'pdist2'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "omp.h"
#include "findNearestLocation_types.h"

/* Function Declarations */
void pdist2(const emlrtStack *sp, const real_T Xin[2086], const real_T Yin[4000],
            real_T D[2000], real_T b_I[2000]);

/* End of code generation (pdist2.h) */
