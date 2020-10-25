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
#include "config_smooth.h"

#ifndef IO_SMOOTH_H
#define IO_SMOOTH_H
void initializeElNode (struct ElNode **lle);
void initializePmNode (struct PmNode **llm);

unsigned readInputMeshAndModel (int ic, char *prm1, char *prm2,
                                unsigned n, struct Profile p[n],
                                struct Boundaries *gb, struct ElNode *lle,
                                struct PmNode *llm, unsigned *nel);

unsigned readOutputMeshAndModel (int ic, char *prm1, char *prm2,
                                 struct Boundaries *gb,
                                 struct Point To[NEL][NX][NY][NZ],
                                 struct Parameters Mo[NEL][NX][NY][NZ]);

unsigned scanMeshAndModel (int ic, char *prm1, char *prm2,
                           unsigned n, struct Profile p[n],
                           struct Boundaries *gb,
                           struct ElNode *lle, struct PmNode *llm,
                           unsigned *nel);

void toArrayElAndPmNodes (unsigned nel, struct ElNode *lle,
                          struct PmNode *llm,
                          struct Point Ti[nel][NX][NY][NZ],
                          struct Parameters M[nel][NX][NY][NZ]);

unsigned readNumberOfPoints (unsigned *n);
unsigned readSmoothingProfile (unsigned n, double *sigma,
                               struct Profile p[n]);

unsigned writeModel (int ic, char *prm,
                     struct Parameters M[NEL][NX][NY][NZ]);

unsigned checkProfileIO (unsigned rvalue);
unsigned checkTopoAndModelIO (unsigned rvalue);
unsigned checkModelIO (unsigned rvalue);
#endif
