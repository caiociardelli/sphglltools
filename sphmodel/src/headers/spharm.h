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

#include "structs.h"

#ifndef SPHARM_H
#define SPHARM_H
void radialBasis (double r1, double r2, unsigned nr,
                  unsigned dg, unsigned nnt, unsigned ns,
                  double T[nnt], double R[nr][ns]);
void basisMatrix (double r1, double r2, unsigned Nr,
                  unsigned dg, unsigned ns,
                  double R[Nr][ns], double Bm[ns][ns]);

void polarBasis (unsigned nt, double Theta[nt],
                 unsigned n,
                 unsigned nlg, unsigned nmax,
                 long double nF[nlg],
                 long double P[nt][nmax + 2],
                 long double nP[nt][nmax + 2]);
void polarBasisPf (double t, unsigned n,
                   unsigned nlg, unsigned nmax,
                   long double nF[nlg],
                   long double P[nmax + 2],
                   long double nP[nmax + 2]);

void azimuthalBasis (unsigned np, double Phi[np],
                     unsigned nmax,
                     double C[np][nmax + 1],
                     double S[np][nmax + 1]);
void azimuthalBasisPf (double p, unsigned nmax,
                       double C[nmax + 1],
                       double S[nmax + 1]);

void initializeArray (unsigned n, double M[n]);
void initialize2DArray (unsigned n1, unsigned n2, double M[n1][n2]);
void initializeMeanArray (unsigned n, struct MeanModel M[n]);
#endif
