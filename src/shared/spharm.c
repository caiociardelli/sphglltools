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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"

void radialBasis (double r1, double r2, unsigned Nr,
                  unsigned dg, unsigned nnt, unsigned ns,
                  double T[nnt], double R[Nr][ns])
{
  /* Computes the radial basis */
  double dr = (r2 - r1) / (Nr - 1);

  for (unsigned s = 0; s < ns; s++)
  {
    double r = r1;

    for (unsigned i = 0; i < Nr; i++)
    {
      R[i][s] = bSplines (r, nnt, dg, T, s);

      r += dr;
    }
  }

  R[Nr - 1][ns - 1] = 1.0;
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

void polarBasis (double t1, double t2, unsigned Nt, unsigned n,
                 unsigned nlg, unsigned nmax, long double nF[nlg],
                 long double P[Nt][nmax + 2], long double nP[Nt][nmax + 2])
{
  /* Computes the polar basis */
  double dt = (t2 - t1) / (Nt - 1);
  double t  = t1;

  for (unsigned i = 0; i < Nt; i++)
  {
    lgPmn (cosl (t), n, nmax, P[i]);
    normalize (n, nlg, nmax, nF, P[i], nP[i]);

    t += dt;
  }
}

void azimuthalBasis (double p1, double p2, unsigned Np,
                     unsigned nmax,
                     double C[Np][nmax + 1],
                     double S[Np][nmax + 1])
{
  /* Computes the azimuthal basis */
  double dp = (p2 - p1) / (Np - 1);
  double p  = p1;

  for (unsigned k = 0; k < Np; k++)
  {
    for (unsigned m = 0; m <= nmax; m++)
    {
      C[k][m] = cos (m * p);
      S[k][m] = sin (m * p);
    }

    p += dp;
  }
}

void precomputeIntegral (double W[NX][NY][NZ],
                         double M[NEL][NX][NY][NZ],
                         double J[NEL][NX][NY][NZ],
                         double K[NEL][NX][NY][NZ])
{
  /* Precomputes part of the integral to save numerical operations */
  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)

          K[el][i][j][k] = W[i][j][k] * M[el][i][j][k] * J[el][i][j][k];
}

void preComputeRadialIndices (double r1, double r2,
                              struct SphericalPoint Tm[NEL][NX][NY][NZ],
                              unsigned Ri[NEL][NX][NY][NZ])
{
  /* Precomputes indices for the radial basis functions */
  double dr_i = (NPT - 1) / (r2 - r1);

  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          struct SphericalPoint p = Tm[el][i][j][k];

          Ri[el][i][j][k] = v2Index (p.r, r1, dr_i);
        }
}

void preComputeIndices (double r1, double r2,
                        double t1, double t2,
                        double p1, double p2,
                        struct SphericalPoint Tm[NEL][NX][NY][NZ],
                        struct BasisIndices Bi[NEL][NX][NY][NZ])
{
  /* Precomputes indices for the basis functions */
  double dr_i = (NPT - 1) / (r2 - r1);
  double dt_i = (NPT - 1) / (t2 - t1);
  double dp_i = (NPT - 1) / (p2 - p1);

  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          struct SphericalPoint p = Tm[el][i][j][k];

          Bi[el][i][j][k].ri = v2Index (p.r, r1, dr_i);
          Bi[el][i][j][k].ti = v2Index (p.theta, t1, dt_i);
          Bi[el][i][j][k].pi = v2Index (p.phi, p1, dp_i);
        }
}

void initializeArray (double Em[NEL][NX][NY][NZ])
{
  /* Initializes array to recreate the model */
  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)

          Em[el][i][j][k] = 0.0;
}

double computeRMSD (bool In[NEL],
                    double M[NEL][NX][NY][NZ],
                    double Em[NEL][NX][NY][NZ],
                    double W[NX][NY][NZ], double J[NEL][NX][NY][NZ])
{
  /* Computes the rmsd between the original and the recreated model */
  double rs, rs_l = 0.0;
  double v,  v_l  = 0.0;

  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          if (!In[el]) continue;

          rs_l += W[i][j][k] * square (M[el][i][j][k] - Em[el][i][j][k])
                             * J[el][i][j][k];
          v_l  += W[i][j][k] * J[el][i][j][k];
        }

  MPI_Allreduce (&rs_l, &rs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&v_l, &v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sqrt (rs / v);
}

