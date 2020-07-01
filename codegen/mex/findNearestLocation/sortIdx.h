/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sortIdx.h
 *
 * Code generation for function 'sortIdx'
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
void merge_block(const emlrtStack *sp, int32_T idx[1043], real_T x[1043],
                 int32_T offset, int32_T n, int32_T preSortLevel, int32_T iwork
                 [1043], real_T xwork[1043]);
void merge_pow2_block(int32_T idx[1043], real_T x[1043], int32_T offset);

/* End of code generation (sortIdx.h) */
