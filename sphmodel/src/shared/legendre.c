/*
 SphGLLTools/SphModel

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
  /* Computes the associated Legendre polynomials */
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
  /* Normalizes the associated Legendre polynomials */
  for (unsigned m = 0; m <= n; m++)

    nP[m + 1] = K2 * nF[mN2I (m, n)] * P[m + 1];
}

