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

#include "constants.h"
#include "config.h"

#ifndef IO_H
#define IO_H
void initializeElNode (struct ElNode **lle);
void initializePmNode (struct PmNode **llm);

unsigned scanMeshAndModel (int ic, char *prm, unsigned nelm, struct Boundaries *gb,
                           struct ElNode *lle, struct PmNode *llm, unsigned *nel);

void toArrayElAndPmNodes (unsigned nel, struct ElNode *lle, struct PmNode *llm,
                          struct Point Ti[nel][NX][NY][NZ],
                          struct Parameters M[nel][NX][NY][NZ]);

unsigned readOutputMesh (int ic, char *prm, struct Boundaries *gb,
                         struct Point To[NEL2][NX][NY][NZ]);

unsigned readMeshModelAndJacobian (int ic, char *input, char *prm,
                                   struct SphericalPoint Tm[NEL][NX][NY][NZ],
                                   double M[NEL][NX][NY][NZ],
                                   double J[NEL][NX][NY][NZ]);

unsigned readBSplinesHeader (double *rmin, double *rmax,
                             unsigned *nn, unsigned *dg,
                             unsigned zone);
unsigned readKnots (double r1, double r2, unsigned dg,
                    unsigned zone, unsigned nnt, double T[nnt]);

unsigned readMeanModelHeader (char *prm, unsigned *nl);
unsigned readMeanModel (char *prm, unsigned nl,
                        struct MeanModel Mm[nl]);
double aMeanModel (double r, unsigned nl, struct MeanModel Mm[nl]);

unsigned computeGatheringArrays (int ic, int nc,
                                 unsigned n1, unsigned n2, unsigned n3,
                                 int k[nc], int kk[nc]);

void createOutputGridBk (int ic, int nc, int k[nc], int kk[nc],
                         long double r1, long double r2,
                         unsigned np, unsigned np_l, unsigned nt, unsigned nr,
                         struct Boundaries *gb, struct Point To[np_l][nt][nr]);
void createOutputGridDD (int ic, int nc, int k[nc], int kk[nc],
                         unsigned nd, unsigned nd_l, unsigned nr,
                         long double r1, long double r2,
                         long double t1, long double t2,
                         long double p1, long double p2,
                         long double R[nr], long double Delta_l[nd_l],
                         struct Boundaries *gb, struct Point To[nd_l][nr]);
void createOutputGridLL (int ic, int nc, int k[nc], int kk[nc],
                         long double r,
                         unsigned np, unsigned np_l, unsigned nt,
                         struct Boundaries *gb, struct Point To[np_l][nt]);
void createOutputGridPf (long double r1, long double r2,
                         long double t, long double p,
                         unsigned nr,
                         long double R[nr],
                         struct Boundaries *gb, struct Point To[nr]);
void createOutputGrid1D (int ic, int nc, int k[nc], int kk[nc],
                         long double r1, long double r2,
                         unsigned np, unsigned np_l, unsigned nt, unsigned nr,
                         struct Boundaries *gb,
                         long double R[nr],
                         long double Theta[nt],
                         long double Phi[np_l]);

unsigned writeModelBk (char *prm,
                       unsigned np, unsigned nt, unsigned nr,
                       struct Parameters Mo[np][nt][nr]);
unsigned writeModelDD (char *argv[], unsigned nd, unsigned nr,
                       long double R[nr], long double Delta[nd],
                       struct Parameters Mo[nd][nr]);
unsigned writeModelGLL (int ic, char *prm,
                        struct Parameters M[NEL2][NX][NY][NZ]);
unsigned writeModelLL (char *prm, struct Boundaries *gb,
                       long double depth, unsigned np, unsigned nt,
                       struct Parameters Mo[np][nt]);
unsigned writeModelPf (char *argv[], unsigned nr,
                       long double R[nr], struct Parameters Mo[nr]);
unsigned writeModel1D (char *prm, unsigned nr,
                       long double R[nr], struct Means Mo[nr]);

unsigned writeCoefficients (char *output, char *prm, unsigned zone,
                            unsigned nmax, unsigned ns, unsigned nlg,
                            double A[ns][nlg], double B[ns][nlg]);

unsigned writeKnotsAndCoeffs (char *output, char *prm, unsigned zone,
                              unsigned nnt, double T[nnt],
                              unsigned ns, double C[ns]);

unsigned checkTopoIO (unsigned rvalue);
unsigned checkTopoAndModelIO (unsigned rvalue);
unsigned checkModelIO (unsigned rvalue);
unsigned checkBSplinesHeaderIO (unsigned rvalue);
unsigned checkKnotsIO (unsigned rvalue);
unsigned checkMeanModelHeaderIO (unsigned rvalue);
unsigned checkMeanModelIO (unsigned rvalue);
unsigned checkCoefficientsIO (unsigned rvalue);
#endif
