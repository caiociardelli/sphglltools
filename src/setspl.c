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

-----------------------------------------------------------------------------------------------

 SETSPL

 USAGE
   mpiexec -n 24 bin/setspl PARAMETER ZONE INPUT_DIRECTORY OUTPUT_DIRECTORY

 EXAMPLE
   mpiexec -n 24 bin/setspl vsv 4 data/OUTPUT/ .

 COMMAND-LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vp, vs, rho, eta, vsv, etc.)
   ZONE                   - zone to be expanded (2, 3 or 4)
   INPUT_DIRECTORY        - directory containing the input files
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads the model parameter, the zone (2: upper mantle, 3: transition zone or 4: lower mantle),
   and the input and output directory names from the command line and computes B-splines knots
   and coefficients. The routine writes the output to a file called PARAMETER_ZZONE_KC.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"
#include "boundaries.h"
#include "io.h"
#include "coordinates.h"
#include "spharm.h"
#include "constants.h"
#include "config.h"

static void gllIntegrate (unsigned Ri[NEL][NX][NY][NZ], bool In[NEL],
                          unsigned ns, unsigned s,
                          double K[NEL][NX][NY][NZ], double Rb[NPT][ns],
                          double *intg)
{
  /* Computes volumetric integral through all spectral elements */
  for (unsigned el = 0; el < NEL; el++)
  {
    if (!In[el]) continue;

      for (unsigned i = 0; i < NX; i++)

        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            unsigned ri = Ri[el][i][j][k];

            if (!Rb[ri][s]) continue;

            double v = K[el][i][j][k] * Rb[ri][s];

            *intg += v;
          }
  }
}

static void computeCoefficients (unsigned ns,
                                 unsigned Ri[NEL][NX][NY][NZ],
                                 bool In[NEL], double Bm[ns][ns],
                                 double K[NEL][NX][NY][NZ],
                                 double Rb[NPT][ns], double C[ns])
{
  /* Computes the coefficients of the expansion */
  double Cl[ns];

  double nf = 1.0 / (4.0 * PI);

  for (unsigned s1 = 0; s1 < ns; s1++)
  {
    Cl[s1] = 0.0;

    for (unsigned s2 = 0; s2 < ns; s2++)
    {
      double intg = 0;

      gllIntegrate (Ri, In, ns, s2, K, Rb, &intg);

      Cl[s1] += Bm[s1][s2] * intg;
    }

    Cl[s1] *= nf;
  }

  MPI_Allreduce (Cl, C, ns, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

static unsigned meanModel (double r1, double r2,
                           unsigned nl,
                           struct MeanModel Mm[nl],
                           double Rmm[NPT])
{
  /* Returns the arithmetic mean model value */
  double r = r1, dr = (r2 - r1) / (NPT - 1);

  for (unsigned i = 0; i < NPT; i++)

    Rmm[i] = aMeanModel (r += dr, nl, Mm);

  return 0;
}

static void computeExpansion (double r1, double r2,
                              unsigned ns, unsigned nnt, unsigned dg,
                              double T[nnt], double C[ns], double Em[NPT])
{
  /* Recreates 1D model from the coefficients */
  double Rb[NPT][ns];

  radialBasis (r1, r2, NPT, dg, nnt, ns, T, Rb);

  for (unsigned r = 0; r < NPT; r++)
  {
    Em[r] = 0.0;

    for (unsigned s = 0; s < ns; s++)

      if (Rb[r][s]) Em[r] += Rb[r][s] * C[s];
  }
}

static double misfit (double r1, double r2,
                      double Rmm[NPT], double Em[NPT])
{
  /* Auxiliary funtion to compute the misfit */
  double Qui[NPT];

  for (unsigned r = 0; r < NPT; r++)

    Qui[r] = square (Em[r] - Rmm[r]);

  return simpson (r1, r2, NPT, Qui);
}

static unsigned computeMisfit (double r1, double r2,
                               unsigned ns, unsigned nnt, unsigned dg,
                               double T[nnt], bool In[NEL],
                               double Rmm[NPT], double Rb[NPT][ns],
                               double Bm[ns][ns], unsigned Ri[NEL][NX][NY][NZ],
                               double K[NEL][NX][NY][NZ],
                               double C[ns], double *q)
{
  /* Computes the misfit between the mean model and the recreated
     1D model from the B-spline */
  double Em[NPT];

  radialBasis (r1, r2, NPT, dg, nnt, ns, T, Rb);
  basisMatrix (r1, r2, NPT, dg, ns, Rb, Bm);

  if (gaussJordan (ns, Bm)) return 1;

  computeCoefficients (ns, Ri, In, Bm, K, Rb, C);
  computeExpansion (r1, r2, ns, nnt, dg, T, C, Em);

  *q = misfit (r1, r2, Rmm, Em);

  return 0;
}

static void helpMenu (void)
{
  char *help_menu = "\n SETSPL"

                    "\n\n USAGE"
                    "\n    mpiexec -n 24 bin/setspl PARAMETER ZONE INPUT_DIRECTORY OUTPUT_DIRECTORY"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 24 bin/setspl vsv 4 data/OUTPUT/ ."

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vp, vs, rho, eta, vsv, etc.)"
                    "\n    ZONE                   - zone to be expanded (2, 3 or 4)"
                    "\n    INPUT_DIRECTORY        - directory containing the input files"
                    "\n    OUTPUT_DIRECTORY       - directory where the routine will write the output files"

                    "\n\n DESCRIPTION"
                    "\n    Reads the model parameter, the zone (2: upper mantle, 3: transition zone or 4: lower mantle),"
                    "\n    and the input and output directory names from the command line and computes B-splines knots"
                    "\n    and coefficients. The routine writes the output to a file called PARAMETER_ZZONE_KC.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  int ic;
  int nc;

  MPI_Init (NULL, NULL);

  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  if (argc != 5 && ic == 0)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command line...\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  char *prm = argv[1];
  unsigned zone = atoi (argv[2]);

  double X[NX], Y[NY], Z[NZ];
  double W[NX][NY][NZ];

  gLLNodesAndWeights3D (X, Y, Z, W);

  struct SphericalPoint Tm[NEL][NX][NY][NZ];

  double M[NEL][NX][NY][NZ];
  double J[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nReading topological and model"
                                " files and computing jacobian...\n");

  if (checkTopoAndModelIO (readMeshModelAndJacobian (ic, argv[3], prm,
                                                     Tm, M, J)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  unsigned nn;
  unsigned dg;

  double rmin;
  double rmax;

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "Reading B-splines information from"
                     " 'knots_Zone%u.dat' file...\n", zone);

    if (checkBSplinesHeaderIO (readBSplinesHeader (&rmin, &rmax,
                                                   &nn, &dg,
                                                   zone)))

      MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Bcast (&nn, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast (&dg, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast (&rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  bool In[NEL];

  double r1, r2;

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nDetecting boundaries of the domain...\n");

  detectBoundaries (rmin, rmax, Tm, In, &r1, &r2);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Expanding from depth: %.1lf km to %.1lf km\n",
                        r2Depth (r1), r2Depth (r2));

  unsigned ns  = nnAndDg2Ns (nn, dg);
  unsigned nnt = nnAndDg2Nnt (nn, dg);

  double T[nnt];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\nCreating knots array...\n");

    if (checkKnotsIO (readKnots (r1, r2, dg,
                                 zone, nnt, T))) MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Bcast (T, nnt, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double K[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing part of the integrals...\n");

  precomputeIntegral (W, M, J, K);

  unsigned Ri[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing indices...\n");

  preComputeRadialIndices (r1, r2, Tm, Ri);

  unsigned nl;

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\nReading mean model...\n");

    if (checkMeanModelHeaderIO (readMeanModelHeader (prm, &nl)))

      MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Bcast (&nl, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  struct MeanModel Mm[nl];

  if (ic == 0 && checkMeanModelIO (readMeanModel (prm, nl, Mm)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  MPI_Datatype PARAMETERS, types[] = {MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_DOUBLE};
  int blocks[] = {1, 1, 1};

  MPI_Aint displ[] = {offsetof (struct MeanModel, r),
                      offsetof (struct MeanModel, vam),
                      offsetof (struct MeanModel, vgm)};

  MPI_Type_create_struct (3, blocks, displ, types, &PARAMETERS);
  MPI_Type_commit (&PARAMETERS);

  MPI_Bcast (Mm, nl, PARAMETERS, 0, MPI_COMM_WORLD);

  MPI_Type_free (&PARAMETERS);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing mean model...\n");

  double Rmm[NPT];

  if (meanModel (r1, r2, nl, Mm, Rmm))
  {
    fprintf (stderr, "Error: could not read mean model file!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  double Rb[NPT][ns];
  double Bm[ns][ns];

  double q = INFINITY, C[ns];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nComputing misfit...\n");

  if (computeMisfit (r1, r2, ns, nnt, dg, T, In, Rmm, Rb, Bm, Ri, K, C, &q))
  {
    fprintf (stderr, "Error: radial basis matrix is singular!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  if (ic == 0) fprintf (stderr, "\nMisfit = %E\n", q);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\nWriting B-splines knots and coefficients...\n");

    if (checkCoefficientsIO (writeKnotsAndCoeffs (argv[4], prm, zone,
                                                  nnt, T, ns, C)))

      MPI_Abort (MPI_COMM_WORLD, 1);

    fprintf (stderr, "Done!\n");
  }

  MPI_Finalize ();

  return 0;
}

