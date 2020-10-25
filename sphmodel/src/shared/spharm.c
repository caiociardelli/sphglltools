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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"

void radialBasis (double r1, double r2, unsigned nr,
                  unsigned dg, unsigned nnt, unsigned ns,
                  double T[nnt], double R[nr][ns])
{
  /* Computes radial basis */
  double dr = (r2 - r1) / (nr - 1);

  for (unsigned s = 0; s < ns; s++)
  {
    double r = r1;

    for (unsigned i = 0; i < nr; i++)
    {
      R[i][s] = bSplines (r, nnt, dg, T, s);

      r += dr;
    }
  }

  R[nr - 1][ns - 1] = 1.0;
}

void basisMatrix (double r1, double r2, unsigned Nr,
                  unsigned dg, unsigned ns,
                  double R[Nr][ns], double Bm[ns][ns])
{
  /* Computes the radial basis matrix */
  double dr = (r2 - r1) / (Nr - 1);

  double f[Nr];

  for (unsigned i = 0; i < ns; i++)

    for (unsigned j = 0; j < ns; j++)
    {
      Bm[i][j] = 0;

      if (abs (i - j) > dg) continue;

      double r = r1;

      for (unsigned k = 0; k < Nr; k++)
      {
        f[k] = R[k][i] * R[k][j] * square (r);

        r += dr;
      }

      Bm[i][j] = simpson (r1, r2, Nr, f);
    }
}

void polarBasis (unsigned nt, double Theta[nt],
                 unsigned n,
                 unsigned nlg, unsigned nmax,
                 long double nF[nlg],
                 long double P[nt][nmax + 2],
                 long double nP[nt][nmax + 2])
{
  /* Computes polar basis */
  for (unsigned i = 0; i < nt; i++)
  {
    double t = Theta[i];

    lgPmn (cosl (t), n, nmax, P[i]);
    normalize (n, nlg, nmax, nF, P[i], nP[i]);
  }
}

void polarBasisPf (double t, unsigned n,
                   unsigned nlg, unsigned nmax,
                   long double nF[nlg],
                   long double P[nmax + 2],
                   long double nP[nmax + 2])
{
  /* Computes polar basis for 1D
     profile */
  lgPmn (cos (t), n, nmax, P);
  normalize (n, nlg, nmax, nF, P, nP);
}

void azimuthalBasis (unsigned np, double Phi[np],
                     unsigned nmax,
                     double C[np][nmax + 1],
                     double S[np][nmax + 1])
{
  /* Computes azimuthal basis */
  for (unsigned i = 0; i < np; i++)
  {
    double p = Phi[i];

    for (unsigned m = 0; m <= nmax; m++)
    {
      C[i][m] = cos (m * p);
      S[i][m] = sin (m * p);
    }
  }
}

void azimuthalBasisPf (double p, unsigned nmax,
                       double C[nmax + 1],
                       double S[nmax + 1])
{
  /* Computes azimuthal basis for
     1D profile */
  for (unsigned m = 0; m <= nmax; m++)
  {
    C[m] = cos (m * p);
    S[m] = sin (m * p);
  }
}

void initializeArray (unsigned n, double M[n])
{
  /* Initializes array to recreate model */
  for (unsigned i = 0; i < n; i++)

    M[i] = 0.0;
}

void initialize2DArray (unsigned n1, unsigned n2, double M[n1][n2])
{
  /* Initializes 2D array to recreate model */
  for (unsigned i = 0; i < n1; i++)

    for (unsigned j = 0; j < n2; j++)

      M[i][j] = 0.0;
}

void initializeMeanArray (unsigned n, struct MeanModel M[n])
{
  /* Initializes mean model array */
  for (unsigned i = 0; i < n; i++)
  {
      M[i].vam = 0.0;
  }
}

