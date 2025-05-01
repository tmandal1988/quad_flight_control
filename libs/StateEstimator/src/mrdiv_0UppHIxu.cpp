//
// File: mrdiv_0UppHIxu.cpp
//
// Code generated for Simulink model 'stateEstimator'.
//
// Model version                  : 1.375
// Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
// C/C++ source code generated on : Tue Apr 29 15:53:54 2025
//
#include "rtwtypes.h"
#include "mrdiv_0UppHIxu.h"
#include <cstring>
#include <cmath>

// Function for MATLAB Function: '<S1>/EKF'
void mrdiv_0UppHIxu(const real32_T A[138], const real32_T B[36], real32_T Y[138])
{
  int32_T c_ix;
  int32_T ijA;
  int32_T ix;
  int32_T iy;
  int32_T jj;
  int32_T kBcol;
  real32_T b_A[36];
  real32_T smax;
  int8_T ipiv[6];
  std::memcpy(&b_A[0], &B[0], 36U * sizeof(real32_T));
  std::memcpy(&Y[0], &A[0], 138U * sizeof(real32_T));
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  ipiv[4] = 5;
  ipiv[5] = 6;
  for (int32_T d_j{0}; d_j < 5; d_j++) {
    jj = d_j * 7;
    iy = 0;
    ix = jj;
    smax = std::abs(b_A[jj]);
    for (kBcol = 2; kBcol <= 6 - d_j; kBcol++) {
      real32_T s;
      ix++;
      s = std::abs(b_A[ix]);
      if (s > smax) {
        iy = kBcol - 1;
        smax = s;
      }
    }

    if (b_A[jj + iy] != 0.0F) {
      if (iy != 0) {
        iy += d_j;
        ipiv[d_j] = static_cast<int8_T>(iy + 1);
        smax = b_A[d_j];
        b_A[d_j] = b_A[iy];
        b_A[iy] = smax;
        smax = b_A[d_j + 6];
        b_A[d_j + 6] = b_A[iy + 6];
        b_A[iy + 6] = smax;
        smax = b_A[d_j + 12];
        b_A[d_j + 12] = b_A[iy + 12];
        b_A[iy + 12] = smax;
        smax = b_A[d_j + 18];
        b_A[d_j + 18] = b_A[iy + 18];
        b_A[iy + 18] = smax;
        smax = b_A[d_j + 24];
        b_A[d_j + 24] = b_A[iy + 24];
        b_A[iy + 24] = smax;
        smax = b_A[d_j + 30];
        b_A[d_j + 30] = b_A[iy + 30];
        b_A[iy + 30] = smax;
      }

      iy = (jj - d_j) + 6;
      for (ix = jj + 1; ix < iy; ix++) {
        b_A[ix] /= b_A[jj];
      }
    }

    iy = jj;
    ix = jj + 6;
    for (kBcol = 0; kBcol <= 4 - d_j; kBcol++) {
      if (b_A[ix] != 0.0F) {
        int32_T c;
        smax = -b_A[ix];
        c_ix = jj + 1;
        ijA = iy + 7;
        c = (iy - d_j) + 12;
        while (ijA + 1 <= c) {
          b_A[ijA] += b_A[c_ix] * smax;
          c_ix++;
          ijA++;
        }
      }

      ix += 6;
      iy += 6;
    }
  }

  for (int32_T d_j{0}; d_j < 6; d_j++) {
    jj = 23 * d_j;
    iy = 6 * d_j;
    for (ix = 0; ix < d_j; ix++) {
      kBcol = 23 * ix;
      if (b_A[ix + iy] != 0.0F) {
        for (c_ix = 0; c_ix < 23; c_ix++) {
          ijA = c_ix + jj;
          Y[ijA] -= b_A[ix + iy] * Y[c_ix + kBcol];
        }
      }
    }

    smax = 1.0F / b_A[d_j + iy];
    for (iy = 0; iy < 23; iy++) {
      ijA = iy + jj;
      Y[ijA] *= smax;
    }
  }

  for (int32_T d_j{5}; d_j >= 0; d_j--) {
    jj = 23 * d_j;
    iy = 6 * d_j - 1;
    for (ix = d_j + 2; ix < 7; ix++) {
      kBcol = (ix - 1) * 23;
      if (b_A[ix + iy] != 0.0F) {
        for (c_ix = 0; c_ix < 23; c_ix++) {
          ijA = c_ix + jj;
          Y[ijA] -= b_A[ix + iy] * Y[c_ix + kBcol];
        }
      }
    }
  }

  for (int32_T d_j{4}; d_j >= 0; d_j--) {
    int8_T ipiv_0;
    ipiv_0 = ipiv[d_j];
    if (d_j + 1 != ipiv_0) {
      for (iy = 0; iy < 23; iy++) {
        smax = Y[23 * d_j + iy];
        ijA = (ipiv_0 - 1) * 23 + iy;
        Y[iy + 23 * d_j] = Y[ijA];
        Y[ijA] = smax;
      }
    }
  }
}

//
// File trailer for generated code.
//
// [EOF]
//
