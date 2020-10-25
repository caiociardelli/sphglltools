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
#include "constants.h"

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

