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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "exmath.h"
#include "structs_smooth.h"

static unsigned getIndex (double r, unsigned n, struct Profile p[n])
{
  /* Converts radius to index */
  if (n == 1 || r > p[1].r) return 0;

  for (unsigned i = 1; i < n - 1; i++)

    if (r <= p[i].r && r > p[i + 1].r) return i;

  return n - 1;
}

static void slidingCorrelation (unsigned m, double weights[m],
                                struct Smooth sm[NP])
{
  /* Applies smoothing kernel to the smoothing profile */
  int n = (int) NP;
  unsigned hw = m / 2;

  struct Smooth c[NP];

  memcpy (c, sm, NP * sizeof (struct Smooth));

  for (unsigned i = 0; i < NP; i++)
  {
    sm[i].sigma_h = 0;
    sm[i].sigma_v = 0;

    for (unsigned j = 0; j < m; j++)
    {
      int k = i + j - hw;

      if (k >= n) k = n - 1; else if (k < 0) k = 0;

      sm[i].sigma_h += c[k].sigma_h * weights[j];
      sm[i].sigma_v += c[k].sigma_v * weights[j];
    }

    sm[i].search_h = TRUNCATE * sm[i].sigma_h;
    sm[i].search_v = TRUNCATE * sm[i].sigma_v;
    sm[i].sg_h     = 0.5 / square (sm[i].sigma_h);
    sm[i].sg_v     = 0.5 / square (sm[i].sigma_v);
    sm[i].sc_h     = 1.0 / square (sm[i].search_h);
    sm[i].sc_v     = 1.0 / square (sm[i].search_v);
  }
}

void smoothPf (double sigma, struct Smooth sm[NP])
{
  /* Smooths smoothing profile */
  sigma = NP * sigma / EARTH_R;

  unsigned lw = (unsigned) (TRUNCATE * sigma + 0.5);
  unsigned m  = 2 * lw + 1;

  double sum = 1.0, sigma2 = square (sigma), w, weights[m];

  weights[lw] = 1.0;

  for (unsigned i = 1; i < lw + 1; i++)
  {
    w = exp (-0.5 * square (i) / sigma2);

    weights[lw + i] = w;
    weights[lw - i] = w;

    sum += 2.0 * w;
  }

  for (unsigned i = 0; i < m; i++)

    weights[i] /= sum;

  slidingCorrelation (m, weights, sm);
}

void interpolate (unsigned n, struct Profile p[n],
                  struct Smooth sm[NP])
{
  /* Interpolates smoothing profile */
  const double dr = (MAX_SURFACE_R - CMB_R) / (NP - 1);

  for (unsigned i = 0; i < NP; i++)
  {
    unsigned j = getIndex (CMB_R + i * dr, n, p);

    sm[i].sigma_h = p[j].sigma_h;
    sm[i].sigma_v = p[j].sigma_v;
  }
}

