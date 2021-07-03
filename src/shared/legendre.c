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
#include "exmath.h"
#include "legendre.h"
#include "constants.h"

double lgPn (unsigned n, double x)
{
  /* Computes the Lengendre polynomials Pn (x)*/
  switch (n)
  {
    case 0:

      return 1.0;

    case 1:

      return x;

    default:
    {
      double fPn = 1.0, sPn = x, nPn = 0.0;

      for (unsigned i = 2; i <= n; i++)
      {
        nPn = ((2 * i - 1) * x * sPn - (i - 1) * fPn) / i;

        fPn = sPn; sPn = nPn;
      }

      return nPn;
    }
  }
}

long double lgPnl (unsigned n, long double x)
{
  /* Computes the Lengendre polynomials Pn (x)*/
  switch (n)
  {
    case 0:

      return 1.L;

    case 1:

      return x;

    default:
    {
      long double fPn = 1.L, sPn = x, nPn = 0.L;

      for (unsigned i = 2; i <= n; i++)
      {
        nPn = ((2 * i - 1) * x * sPn - (i - 1) * fPn) / i;

        fPn = sPn; sPn = nPn;
      }

      return nPn;
    }
  }
}

static void gLLNodesAndWeights (unsigned n, double X[n], double W[n])
{
  /* Computes nodes and weights for GLL integration */
  double k1 = 1 - (3.0 / 8.0) * (n - 2) / cube (n - 1);
  double k2 = 4 * (n - 1) + 1;

  X[0] = -1.0; X[n - 1] = 1.0;
  W[0] =  2.0 / ((n * (n - 1))); W[n - 1] = W[0];

  unsigned n_2 = n / 2;

  for (unsigned i = 2; i <= n_2; i++)
  {
    double x = k1 * cos ((4 * i - 3) * PI / k2);

    for (unsigned itr = 0; itr < MAX_HITR; itr++)
    {
      double k3 = 1.0 / (1.0 - square (x));

      double y = lgPn (n - 1, x);

      double y2 = dlgPn (n - 1, x, y, k3);
      double y3 = d2lgPn (n - 1, x, y, y2, k3);
      double y4 = d3lgPn (n - 1, x, y2, y3, k3);

      double dx = -2.0 * (y2 * y3) / (2.0 * square (y3) - y2 * y4);

      x += dx;

      if (fabs (dx) < TOLERANCE) break;
    }

      X[i - 1] = -x;
      X[n - i] = x;

      W[i - 1] = 2.0 / ((n * (n - 1)) * square (lgPn (n - 1, X[i - 1])));
      W[n - i] = W[i - 1];
    }

    if (n % 2)
    {
      X[n_2] = 0.0;
      W[n_2] = 2.0 / ((n * (n - 1)) * square (lgPn (n - 1, X[n_2])));
    }
}

void gLLNodesAndWeights3D (double X[NX], double Y[NX], double Z[NX],
                           double W[NX][NY][NZ])
{
  /* Computes nodes and weights for 3D GLL integration */
  double Wx[NX], Wy[NY], Wz[NZ];

  gLLNodesAndWeights (NX, X, Wx);
  gLLNodesAndWeights (NY, Y, Wy);
  gLLNodesAndWeights (NZ, Z, Wz);

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)

        W[i][j][k] = Wx[i] * Wy[j] * Wz[k];
}

static void gLLNodesAndWeightsl (unsigned n, long double X[n], long double W[n])
{
  /* Computes nodes and weights for GLL integration */
  long double k1 = 1 - (3.L / 8.L) * (n - 2) / cubel (n - 1);
  long double k2 = 4 * (n - 1) + 1;

  X[0] = -1.L; X[n - 1] = 1.L;
  W[0] =  2.L / ((n * (n - 1))); W[n - 1] = W[0];

  unsigned n_2 = n / 2;

  for (unsigned i = 2; i <= n_2; i++)
  {
    long double x = k1 * cosl ((4 * i - 3) * PI / k2);

    for (unsigned itr = 0; itr < MAX_HITR; itr++)
    {
      long double k3 = 1.L / (1.L - squarel (x));

      long double y = lgPnl (n - 1, x);

      long double y2 = dlgPnl (n - 1, x, y, k3);
      long double y3 = d2lgPnl (n - 1, x, y, y2, k3);
      long double y4 = d3lgPnl (n - 1, x, y2, y3, k3);

      long double dx = -2.L * (y2 * y3) / (2.L * squarel (y3) - y2 * y4);

      x += dx;

      if (fabsl (dx) < TOLERANCE) break;
    }

      X[i - 1] = -x;
      X[n - i] = x;

      W[i - 1] = 2.L / ((n * (n - 1)) * squarel (lgPnl (n - 1, X[i - 1])));
      W[n - i] = W[i - 1];
    }

    if (n % 2)
    {
      X[n_2] = 0.L;
      W[n_2] = 2.L / ((n * (n - 1)) * squarel (lgPnl (n - 1, X[n_2])));
    }
}

void gLLNodesAndWeights3Dl (long double X[NX],
                            long double Y[NX],
                            long double Z[NX], long double W[NX][NY][NZ])
{
  /* Computes nodes and weights for 3D GLL integration */
  long double Wx[NX], Wy[NY], Wz[NZ];

  gLLNodesAndWeightsl (NX, X, Wx);
  gLLNodesAndWeightsl (NY, Y, Wy);
  gLLNodesAndWeightsl (NZ, Z, Wz);

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)

        W[i][j][k] = Wx[i] * Wy[j] * Wz[k];
}

static inline long double sgn (unsigned n)
{
  /* Condon–Shortley phase */
  return (n % 2) ? -1.L : 1.L;
}

static long double factor (unsigned m, unsigned n)
{
  /* Auxiliary function to compute the normalization
     factor */
  if (m == 0) return K1;

  long double r = n - m + 1;

  for (unsigned i = n - m + 2; i <= n + m; i++)

    r *= i;

  return K1 / r;
}

static long double nml (unsigned m, unsigned n)
{
  /* Normalization factor */
  long double k1 = (m == 0) ? 1.L : 2.L;
  long double k2 = (2 * n + 1) / (4.L * PI);
  long double k3 = factor (m, n);

  return sqrtl (k1 * k2 * k3);
}

void nmlFactors (unsigned nmax, unsigned nlg, long double nF[nlg])
{
  /* Precomputes normalization factors */
  for (unsigned n = 0; n <= nmax; n++)
  {
    for (unsigned m = 0; m <= n; m++)

      nF[mN2I (m, n)] = nml (m , n);
  }
}

void lgPmn (long double x, unsigned n, unsigned nmax, long double P[nmax + 2])
{
  /* Computes the Associated Lengendre polynomials */
  long double sint = (fabsl (x) == 1.L) ? 0.L : sqrtl (1.L - square (x));
  long double cott = (fabsl (x) == 1.L) ? 0.L : x / sint;

  if (n == 0) P[1] = 1.L;

  else if (n == 1)
  {
    P[0] = P[1];
    P[1] = x;
    P[2] = -sint;
  }

  else
  {
    long double P0 = P[0];

    P[0] = P[1];
    P[1] = ((2 * n - 1) * x * P[1] - (n - 1) * P0) / n;
    P[n + 1] = sgn (n) * factorial2 (2 * n - 1) * powl (sint, n);
    P[n] *= x * (2 * (n - 1) + 1);

    for (unsigned m = n - 1; m > 1; m--)

      P[m] = (-2.0 * m * cott * P[m + 1] - P[m + 2]) / ((n + m) * (n - m + 1));
  }
}

void normalize (unsigned n, unsigned nlg, unsigned nmax,
                long double nF[nlg],
                long double P[nmax + 2], long double nP[nmax + 2])
{
  /* Normalizes the Associated Legendre polynomials */
  for (unsigned m = 0; m <= n; m++)

    nP[m + 1] = K2 * nF[mN2I (m, n)] * P[m + 1];
}

