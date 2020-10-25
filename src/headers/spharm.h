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

#include "structs.h"

#ifndef SPHARM_H
#define SPHARM_H
void radialBasis (double r1, double r2, unsigned Nr,
                  unsigned dg, unsigned nnt, unsigned ns,
                  double T[nnt], double R[Nr][ns]);
void basisMatrix (double r1, double r2, unsigned Nr,
                  unsigned dg, unsigned ns,
                  double R[Nr][ns], double Bm[ns][ns]);

void polarBasis (double t1, double t2, unsigned Nt,
                 unsigned n, unsigned nlg, unsigned nmax,
                 long double nF[nlg],
                 long double P[Nt][nmax + 2],
                 long double nP[Nt][nmax + 2]);

void azimuthalBasis (double p1, double p2, unsigned Np,
                     unsigned nmax,
                     double C[Np][nmax + 1],
                     double S[Np][nmax + 1]);

void precomputeIntegral (double W[NX][NY][NZ],
                         double M[NEL][NX][NY][NZ],
                         double J[NEL][NX][NY][NZ],
                         double K[NEL][NX][NY][NZ]);

void preComputeRadialIndices (double r1, double r2,
                              struct SphericalPoint Tm[NEL][NX][NY][NZ],
                              unsigned Ri[NEL][NX][NY][NZ]);
void preComputeIndices (double r1, double r2,
                        double t1, double t2,
                        double p1, double p2,
                        struct SphericalPoint Tm[NEL][NX][NY][NZ],
                        struct BasisIndices Bi[NEL][NX][NY][NZ]);

void initializeArray (double Em[NEL][NX][NY][NZ]);

double computeRMSD (bool In[NEL],
                    double M[NEL][NX][NY][NZ],
                    double Em[NEL][NX][NY][NZ],
                    double W[NX][NY][NZ], double J[NEL][NX][NY][NZ]);
#endif
