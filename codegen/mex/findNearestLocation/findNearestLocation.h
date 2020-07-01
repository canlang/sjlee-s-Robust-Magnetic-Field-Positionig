/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * findNearestLocation.h
 *
 * Code generation for function 'findNearestLocation'
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
void findNearestLocation(const emlrtStack *sp, const real_T locs_1[2086], const
  real_T locs_2[4000], real_T distance[2000], real_T b_index[2000]);

/* End of code generation (findNearestLocation.h) */
