/*
 SphGLLTools

 Author: Caio Ciardelli, University of São Paulo, October 2020

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

----------------------------------------------------------------------------------------------- */

#include <math.h>
#include "constants.h"

long double powerl (long double v, unsigned n)
{
  /* Computes v raised to the power of n */
  long double r = 1.L;

  for (unsigned i = 0; i < n; i++)

    r *= v;

  return r;
}

double roundV (double v)
{
  /* Rounds a value to the numerical precision set
     on the constant ROUND */
  double pw = 1;

  for (unsigned i = 0; i < ROUND; i++) pw *= 10;

  return round (v * pw) / pw;
}

long double factorial (unsigned n)
{
  /* Computes the factorial of n */
  long double r = 1;

  for (unsigned i = 2; i <= n; i++)

    r *= i;

  return r;
}

long double factorial2 (unsigned n)
{
  /* Computes the double factorial of n */
  long double r = 1;

  for (int i = n; i >= 0; i -= 2)
  {
    if (i == 0 || i == 1) return r;

    else r *= i;
  }

  return r;
}

struct Pair
{
  /* Used to store the permutations carried out
     during the matrix inversion */
  unsigned j1;
  unsigned j2;
};

static void swapRows (unsigned n,
                      double M[n][n], double T[n][n],
                      unsigned i1, unsigned i2)
{
  /* Swaps two rows of a square matrix */
  double tp;

  for (unsigned j = 0; j < n; j++)
  {
    tp = M[i1][j];
    M[i1][j] = M[i2][j];
    M[i2][j] = tp;

    tp = T[i1][j];
    T[i1][j] = T[i2][j];
    T[i2][j] = tp;
  }
}

static void swapRowsl (unsigned n,
                       long double M[n][n], long double T[n][n],
                       unsigned i1, unsigned i2)
{
  /* Swaps two rows of a square matrix */
  long double tp;

  for (unsigned j = 0; j < n; j++)
  {
    tp = M[i1][j];
    M[i1][j] = M[i2][j];
    M[i2][j] = tp;

    tp = T[i1][j];
    T[i1][j] = T[i2][j];
    T[i2][j] = tp;
  }
}

static void swapColumns (unsigned n,
                         double M[n][n], double T[n][n],
                         unsigned j1, unsigned j2)
{
  /* Swaps two columns of a square matrix */
  double tp;

  for (unsigned i = 0; i < n; i++)
  {
    tp = M[i][j1];
    M[i][j1] = M[i][j2];
    M[i][j2] = tp;

    tp = T[i][j1];
    T[i][j1] = T[i][j2];
    T[i][j2] = tp;
  }
}

static void swapColumnsl (unsigned n,
                          long double M[n][n], long double T[n][n],
                          unsigned j1, unsigned j2)
{
  /* Swaps two columns of a square matrix */
  long double tp;

  for (unsigned i = 0; i < n; i++)
  {
    tp = M[i][j1];
    M[i][j1] = M[i][j2];
    M[i][j2] = tp;

    tp = T[i][j1];
    T[i][j1] = T[i][j2];
    T[i][j2] = tp;
  }
}

static void findPivot (unsigned n, double M[n][n], unsigned k,
                       unsigned *imax, unsigned *jmax)
{
  /* Finds the pivot of a square matrix */
  double max = fabs (M[k][k]);

  for (unsigned i = k; i < n; i++)

    for (unsigned j = k; j < n; j++)
    {
      double v = fabs (M[i][j]);

      if (v > max)
      {
          max = v;
        *imax = i;
        *jmax = j;
      }
    }
}

static void findPivotl (unsigned n, long double M[n][n], unsigned k,
                        unsigned *imax, unsigned *jmax)
{
  /* Finds the pivot of a square matrix */
  long double max = fabsl (M[k][k]);

  for (unsigned i = k; i < n; i++)

    for (unsigned j = k; j < n; j++)
    {
      long double v = fabsl (M[i][j]);

      if (v > max)
      {
          max = v;
        *imax = i;
        *jmax = j;
      }
    }
}

static void pivotMatrix (unsigned n,
                         double M[n][n], double T[n][n],
                         struct Pair S[n], unsigned k)
{
  /* Carries out complete pivoting */
  unsigned imax = k;
  unsigned jmax = k;

  findPivot (n, M, k, &imax, &jmax);

  if (imax != k) swapRows (n, M, T, k, imax);
  if (jmax != k) swapColumns (n, M, T, k, jmax);

  S[k].j1 = k;
  S[k].j2 = jmax;
}

static void pivotMatrixl (unsigned n,
                          long double M[n][n], long double T[n][n],
                          struct Pair S[n], unsigned k)
{
  /* Carries out complete pivoting */
  unsigned imax = k;
  unsigned jmax = k;

  findPivotl (n, M, k, &imax, &jmax);

  if (imax != k) swapRowsl (n, M, T, k, imax);
  if (jmax != k) swapColumnsl (n, M, T, k, jmax);

  S[k].j1 = k;
  S[k].j2 = jmax;
}

static void identityMatrix (unsigned n, double T[n][n])
{
  /* Creates an identity matrix */
  for (unsigned i = 0; i < n; i++)

    for (unsigned j = 0; j < n; j++)

      T[i][j] = (i == j) ? 1 : 0;
}

static void identityMatrixl (unsigned n, long double T[n][n])
{
  /* Creates an identity matrix */
  for (unsigned i = 0; i < n; i++)

    for (unsigned j = 0; j < n; j++)

      T[i][j] = (i == j) ? 1.L : 0.L;
}

static void unscramble (unsigned int n, double M[n][n], struct Pair S[n])
{
  /* Unscrambles inverse matrix */
  double tp;

  for (int k = n - 1; k >= 0; k--)
  {
    if (S[k].j1 == S[k].j2) continue;

    for (unsigned i = 0; i < n; i++)
    {
      tp = M[i][S[k].j1];
      M[i][S[k].j1] = M[i][S[k].j2];
      M[i][S[k].j2] = tp;
    }

    for (unsigned j = 0; j < n; j++)
    {
      tp = M[S[k].j1][j];
      M[S[k].j1][j] = M[S[k].j2][j];
      M[S[k].j2][j] = tp;
    }
  }
}

static void unscramblel (unsigned int n, long double M[n][n], struct Pair S[n])
{
  /* Unscrambles inverse matrix */
  long double tp;

  for (int k = n - 1; k >= 0; k--)
  {
    if (S[k].j1 == S[k].j2) continue;

    for (unsigned i = 0; i < n; i++)
    {
      tp = M[i][S[k].j1];
      M[i][S[k].j1] = M[i][S[k].j2];
      M[i][S[k].j2] = tp;
    }

    for (unsigned j = 0; j < n; j++)
    {
      tp = M[S[k].j1][j];
      M[S[k].j1][j] = M[S[k].j2][j];
      M[S[k].j2][j] = tp;
    }
  }
}

unsigned gaussJordan (unsigned n, double M[n][n])
{
  /* Computes the inverse of a matrix using Gauss-Jordan's algorithm */
  double T[n][n];
  struct Pair S[n];

  identityMatrix (n, T);

  for (unsigned k = 0; k < n; k++)
  {
    pivotMatrix (n, M, T, S, k);

    double pivot = M[k][k];

    if (fabs (pivot) < TOLERANCE) return 1;

    for (unsigned j = 0; j < n; j++)
    {
      M[k][j] /= pivot;
      T[k][j] /= pivot;
    }

    for (unsigned i = 0; i < n; i++)
    {
      double f = M[i][k];

      if (fabs (f) > TOLERANCE)

        for (unsigned j = 0; j < n; j++)

          if (i != k)
          {
            M[i][j] -= f * M[k][j];
            T[i][j] -= f * T[k][j];
          }

      M[i][k] = T[i][k];
    }
  }

  unscramble (n, M, S);

  return 0;
}

unsigned gaussJordanl (unsigned n, long double M[n][n])
{
  /* Computes the inverse of a matrix using Gauss-Jordan's algorithm */
  long double T[n][n];
  struct Pair S[n];

  identityMatrixl (n, T);

  for (unsigned k = 0; k < n; k++)
  {
    pivotMatrixl (n, M, T, S, k);

    long double pivot = M[k][k];

    if (fabsl (pivot) < TOLERANCE) return 1;

    for (unsigned j = 0; j < n; j++)
    {
      M[k][j] /= pivot;
      T[k][j] /= pivot;
    }

    for (unsigned i = 0; i < n; i++)
    {
      long double f = M[i][k];

      if (fabsl (f) > TOLERANCE)

        for (unsigned j = 0; j < n; j++)

          if (i != k)
          {
            M[i][j] -= f * M[k][j];
            T[i][j] -= f * T[k][j];
          }

      M[i][k] = T[i][k];
    }
  }

  unscramblel (n, M, S);

  return 0;
}

double simpson (double x1, double x2, unsigned nx, double f[nx])
{
  /* 1D integration using 1/3 Simpson's rule (nx must be odd!) */
  double dx = (x2 - x1) / (nx - 1);

  double sum1 = f[1], sum2 = 0;

  for (unsigned i = 2; i < nx - 2; i += 2)
  {
    sum1 += f[i + 1]; sum2 += f[i];
  }

  return (dx / 3) * (f[0] + f[nx - 1] + 4 * sum1 + 2 * sum2);
}

double bSplines (double x, unsigned nnt, unsigned dg,
                 double T[nnt], unsigned i)
{
  /* Function to compute the b-splines */
  if (dg == 0) return (T[i] <= x && x < T[i + 1]) ? 1.0 : 0.0;

  else
  {
    double a = T[i + dg] - T[i], b = T[i + dg + 1] - T[i + 1];

    double c1 = (a > 0) ? (x - T[i]) / a : 0.0;
    double c2 = (b > 0) ? (T[i + dg + 1] - x) / b : 0.0;

    return c1 * bSplines (x, nnt, dg - 1, T, i) +
           c2 * bSplines (x, nnt, dg - 1, T, i + 1);
  }
}

