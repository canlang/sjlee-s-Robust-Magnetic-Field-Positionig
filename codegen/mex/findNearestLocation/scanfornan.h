/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * scanfornan.h
 *
 * Code generation for function 'scanfornan'
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
void b_scanfornan(const emlrtStack *sp, const real_T X[4000], boolean_T nanobs
                  [2000]);
void scanfornan(const emlrtStack *sp, const real_T X[2086], boolean_T nanobs
                [1043]);

/* End of code generation (scanfornan.h) */
