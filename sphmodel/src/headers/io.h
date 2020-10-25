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
#include "constants.h"

#ifndef IO_H
#define IO_H
unsigned readMeanModelHeader (char *prm, unsigned *nl);
unsigned readMeanModel (char *prm, unsigned nl,
                        struct MeanModel Mm[nl]);

unsigned readBlockModelHeader (char *prm,
                               unsigned *np_b,
                               unsigned *nt_b,
                               unsigned *nr_b);
unsigned readBlockModel (char *prm,
                         unsigned np_b,
                         unsigned nt_b,
                         unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b]);
unsigned readBSplinesHeader (double *rmin, double *rmax,
                             unsigned *nn, unsigned *dg,
                             unsigned zone);
unsigned readKnots (double r1, double r2, unsigned dg,
                    unsigned zone, unsigned nnt, double T[nnt]);

unsigned readCoefficientsHeader (unsigned zone, char *prm,
                                 unsigned *N, unsigned *ns);
unsigned readCoefficients (unsigned zone, char *prm, unsigned N,
                           unsigned ns, unsigned nlg,
                           double A[ns][nlg], double B[ns][nlg]);

void createProfilePf (double r1, double r2,
                      unsigned nr, double R[nr]);
void createOutputGridLL (unsigned nt, unsigned np,
                         double Theta[nt], double Phi[np]);
void createOutputGridDD (double r1, double r2,
                         double t1, double t2,
                         double p1, double p2,
                         unsigned nr, unsigned nd,
                         double R[nr], double Delta[nd],
                         double Theta[nd], double Phi[nd]);

unsigned writeExpansionLL (char *prm, bool dvv,
                           struct SphericalBoundaries *sb,
                           double r,
                           unsigned np, unsigned nt,
                           double M[np][nt]);
unsigned writeExpansionDD (char *argv[], char *prm, bool dvv,
                           unsigned nr, unsigned nd,
                           double R[nr],
                           double Delta[nd],
                           double M[nr][nd]);
unsigned writeExpansionPf (char *argv[], bool dvv,
                           unsigned nr,
                           double R[nr], double M[nr]);
unsigned writeModel1D (char *prm, unsigned nr,
                       double R[nr], struct MeanModel Mo[nr]);

unsigned writePowSpec1D (char *prm, unsigned nmax,
                         double Pws1D2[nmax + 1],
                         double Pws1D3[nmax + 1],
                         double Pws1D4[nmax + 1]);
unsigned writePowSpec2D (char *prm, unsigned nmax,
                         unsigned Nr,
                         double Pws2D2[Nr][nmax + 1],
                         double Pws2D3[Nr][nmax + 1],
                         double Pws2D4[Nr][nmax + 1]);

unsigned checkBlockModelHeaderIO (unsigned rvalue);
unsigned checkBlockModelIO (unsigned rvalue);
unsigned checkBSplinesHeaderIO (unsigned rvalue);
unsigned checkKnotsIO (unsigned rvalue);
unsigned checkMeanModelHeaderIO (unsigned rvalue);
unsigned checkMeanModelIO (unsigned rvalue);
unsigned checkCoefficientsHeaderIO (unsigned rvalue);
unsigned checkCoefficientsIO (unsigned rvalue);
unsigned checkExpansionIO (unsigned rvalue);
unsigned checkPowSpecIO (unsigned rvalue);

void helpMenuLL (void);
void helpMenuDD (void);
void helpMenuPf (void);
#endif
