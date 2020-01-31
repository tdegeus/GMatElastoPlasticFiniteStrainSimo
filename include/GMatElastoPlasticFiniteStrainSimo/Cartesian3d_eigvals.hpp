// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrices
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_EIGVALS_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_EIGVALS_HPP

#include <stdio.h>
#include <math.h>
#include "Cartesian3d.h"

// Macros
#define SQR(x) ((x)*(x)) // x^2

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

// ----------------------------------------------------------------------------
inline int dsyevj3(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double sd, so;                  // Sums of diagonal resp. off-diagonal elements
  double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
  double g, h, z, theta;          // More temporary storage
  double thresh;

  // Initialize Q to the identity matrix
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }

  // Initialize w to diag(A)
  for (int i=0; i < n; i++)
    w[i] = A[i][i];

  // Calculate SQR(tr(A))
  sd = 0.0;
  for (int i=0; i < n; i++)
    sd += fabs(w[i]);
  sd = SQR(sd);

  // Main iteration loop
  for (int nIter=0; nIter < 50; nIter++)
  {
    // Test for convergence
    so = 0.0;
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
        so += fabs(A[p][q]);
    if (so == 0.0)
      return 0;

    if (nIter < 4)
      thresh = 0.2 * so / SQR(n);
    else
      thresh = 0.0;

    // Do sweep
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
      {
        g = 100.0 * fabs(A[p][q]);
        if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
                       &&  fabs(w[q]) + g == fabs(w[q]))
        {
          A[p][q] = 0.0;
        }
        else if (fabs(A[p][q]) > thresh)
        {
          // Calculate Jacobi transformation
          h = w[q] - w[p];
          if (fabs(h) + g == fabs(h))
          {
            t = A[p][q] / h;
          }
          else
          {
            theta = 0.5 * h / A[p][q];
            if (theta < 0.0)
              t = -1.0 / (sqrt(1.0 + SQR(theta)) - theta);
            else
              t = 1.0 / (sqrt(1.0 + SQR(theta)) + theta);
          }
          c = 1.0/sqrt(1.0 + SQR(t));
          s = t * c;
          z = t * A[p][q];

          // Apply Jacobi transformation
          A[p][q] = 0.0;
          w[p] -= z;
          w[q] += z;
          for (int r=0; r < p; r++)
          {
            t = A[r][p];
            A[r][p] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=p+1; r < q; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=q+1; r < n; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[q][r];
            A[q][r] = s*t + c*A[q][r];
          }

          // Update eigenvectors
          for (int r=0; r < n; r++)
          {
            t = Q[r][p];
            Q[r][p] = c*t - s*Q[r][q];
            Q[r][q] = s*t + c*Q[r][q];
          }
        }
      }
  }

  return -1;
}
// ----------------------------------------------------------------------------

template <class U, class V, class W>
void eig(const U& A, V& vec, W& val)
{
    double a[3][3];
    double Q[3][3];
    double w[3];

    std::copy(A.begin(), A.end(), &a[0][0]);

    // use the 'Jacobi' algorithm, which is accurate but not very fast
    // (in practice the faster 'hybrid' "dsyevh3" is too inaccurate for finite elements)
    int succes = dsyevj3(a, Q, w);

    (void)(succes);

    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(succes == 0);

    std::copy(&Q[0][0], &Q[0][0] + 3 * 3, vec.begin());
    std::copy(&w[0], &w[0] + 3, val.begin());
}


template <class U, class V, class W>
void inv_eig(const U& vec, const V& val, W& A)
{
    A(0,0) = val(0)*vec(0,0)*vec(0,0) + val(1)*vec(0,1)*vec(0,1) + val(2)*vec(0,2)*vec(0,2);
    A(0,1) = val(0)*vec(0,0)*vec(1,0) + val(1)*vec(0,1)*vec(1,1) + val(2)*vec(0,2)*vec(1,2);
    A(0,2) = val(0)*vec(0,0)*vec(2,0) + val(1)*vec(0,1)*vec(2,1) + val(2)*vec(0,2)*vec(2,2);
    A(1,1) = val(0)*vec(1,0)*vec(1,0) + val(1)*vec(1,1)*vec(1,1) + val(2)*vec(1,2)*vec(1,2);
    A(1,2) = val(0)*vec(1,0)*vec(2,0) + val(1)*vec(1,1)*vec(2,1) + val(2)*vec(1,2)*vec(2,2);
    A(2,2) = val(0)*vec(2,0)*vec(2,0) + val(1)*vec(2,1)*vec(2,1) + val(2)*vec(2,2)*vec(2,2);
    A(1,0) = A(0,1);
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);
}


}} // namespace ...

#endif

