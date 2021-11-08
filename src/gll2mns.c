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

 GLL2MNS

 USAGE
   mpiexec -n 24 bin/gll2mns PARAMETER NMAX ZONE INPUT_DIRECTORY OUTPUT_DIRECTORY

 EXAMPLE
   mpiexec -n 24 bin/gll2mns vsv 20 4 data/OUTPUT/ .

 COMMAND-LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vp, vs, rho, eta, vsv, etc.)
   NMAX                   - maximum degree of the spherical harmonics expansion
   ZONE                   - zone to be expanded (2, 3 or 4)
   INPUT_DIRECTORY        - directory containing the input files
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads the model parameter, the maximum degree of the expansion, the zone (2: upper mantle,
   3: transition zone or 4: lower mantle), and the input and output directory names from the
   command line and computes the B-splines + spherical harmonics coefficients. The routine
   writes the output to a file called mns_ZZONE_PARAMER.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "spharm.h"
#include "structs.h"
#include "boundaries.h"
#include "io.h"
#include "progress.h"
#include "constants.h"

static void gllIntegrate (struct BasisIndices Bi[NEL][NX][NY][NZ],
                          bool In[NEL], unsigned nmax, unsigned ns,
                          unsigned m, unsigned s,
                          double K[NEL][NX][NY][NZ], double Rb[NPT][ns],
                          long double nP[NPT][nmax + 2],
                          double C[NPT][nmax + 1], double S[NPT][nmax + 1],
                          double *intgA, double *intgB)
{
  /* Computes volumetric integral through all spectral elements */
  *intgA = 0.0, *intgB = 0.0;

  for (unsigned el = 0; el < NEL; el++)
  {
    if (!In[el]) continue;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          struct BasisIndices b = Bi[el][i][j][k];

          if (!Rb[b.ri][s]) continue;

          double v = nP[b.ti][m + 1] * Rb[b.ri][s] * K[el][i][j][k];

          *intgA += v * C[b.pi][m];
          *intgB += v * S[b.pi][m];
        }
  }
}

static void computeCoefficients (int ic,
                                 double t1, double t2,
                                 struct BasisIndices Bi[NEL][NX][NY][NZ],
                                 bool In[NEL], unsigned ns, unsigned nmax,
                                 unsigned nlg, double Bm[ns][ns],
                                 double K[NEL][NX][NY][NZ], double Rb[NPT][ns],
                                 long double nF[nlg],
                                 double C[NPT][nmax + 1], double S[NPT][nmax + 1],
                                 double A[ns][nlg], double B[ns][nlg])
{
  /* Computes the coefficients of the expansion */
  double Al[ns][nlg];
  double Bl[ns][nlg];

  long double P[NPT][nmax + 2];
  long double nP[NPT][nmax + 2];

  unsigned i = 0; clock_t starttime = clock ();

  for (unsigned n = 0; n <= nmax; n++)
  {
    polarBasis (t1, t2, NPT, n, nlg, nmax, nF, P, nP);

    for (unsigned m = 0; m <= n; m++)
    {
      unsigned mni = mN2I (m, n);

      for (unsigned s1 = 0; s1 < ns; s1++)
      {
        Al[s1][mni] = 0.0;
        Bl[s1][mni] = 0.0;

        for (unsigned s2 = 0; s2 < ns; s2++)
        {
          double intgA, intgB;

          gllIntegrate (Bi, In, nmax, ns, m, s2, K,
                        Rb, nP, C, S, &intgA, &intgB);

          Al[s1][mni] += Bm[s1][s2] * intgA;
          Bl[s1][mni] += Bm[s1][s2] * intgB;
        }

        MPI_Barrier (MPI_COMM_WORLD);

        if (ic == 0) progressBar (i++, 10, ns * nlg, starttime);
      }
    }
  }

  MPI_Allreduce (Al, A, ns * nlg, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (Bl, B, ns * nlg, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

static void expandOnGLLMesh (double t1, double t2,
                             unsigned nmax, unsigned ns, unsigned nlg,
                             struct BasisIndices Bi[NEL][NX][NY][NZ],
                             bool In[NEL], double Rb[NPT][ns],
                             long double nF[nlg],
                             double C[NPT][nmax + 1], double S[NPT][nmax + 1],
                             double A[ns][nlg], double B[ns][nlg],
                             double Em[NEL][NX][NY][NZ])
{
  /* Recreates model on the mesh from the coefficients */
  long double P[NPT][nmax + 2];
  long double nP[NPT][nmax + 2];

  for (unsigned n = 0; n <= nmax; n++)
  {
    polarBasis (t1, t2, NPT, n, nlg, nmax, nF, P, nP);

    for (unsigned m = 0; m <= n; m++)

      for (unsigned el = 0; el < NEL; el++)
      {
        if (!In[el]) continue;

        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)

            for (unsigned k = 0; k < NZ; k++)
            {
              struct BasisIndices b = Bi[el][i][j][k];

              double Cmn = nP[b.ti][m + 1] * C[b.pi][m];
              double Smn = nP[b.ti][m + 1] * S[b.pi][m];

              unsigned mni = mN2I (m, n);

              for (unsigned s = 0; s < ns; s++)

                if (Rb[b.ri][s]) Em[el][i][j][k] += Rb[b.ri][s] *
                                                   (Cmn * A[s][mni] +
                                                    Smn * B[s][mni]);
            }
      }
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n GLL2MNS"

                    "\n\n USAGE"
                    "\n    mpiexec -n 24 bin/gll2mns PARAMETER NMAX ZONE INPUT_DIRECTORY OUTPUT_DIRECTORY"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 24 bin/gll2mns vsv 20 4 data/OUTPUT/ ."

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vp, vs, rho, eta, vsv, etc.)"
                    "\n    NMAX                   - maximum degree of the spherical harmonics expansion"
                    "\n    ZONE                   - zone to be expanded (2, 3 or 4)"
                    "\n    INPUT_DIRECTORY        - directory containing the input files"
                    "\n    OUTPUT_DIRECTORY       - directory where the routine will write the output files"

                    "\n\n DESCRIPTION"
                    "\n    Reads the model parameter, the maximum degree of the expansion, the zone (2: upper mantle,"
                    "\n    3: transition zone or 4: lower mantle), and the input and output directory names from the"
                    "\n    command line and computes the B-splines + spherical harmonics coefficients. The routine"
                    "\n    writes the output to a file called mns_ZZONE_PARAMER.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  int ic;
  int nc;

  MPI_Init (NULL, NULL);

  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  if (argc != 6 && ic == 0)
  {
    fprintf (stderr, "Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  unsigned nmax = atoi (argv[2]);
  unsigned zone = atoi (argv[3]);

  double t1 =  0;
  double t2 =  PI;
  double p1 = -PI;
  double p2 =  PI;

  unsigned nlg = N2Nlg (nmax);

  double X[NX], Y[NY], Z[NZ];
  double W[NX][NY][NZ];

  gLLNodesAndWeights3D (X, Y, Z, W);

  struct SphericalPoint Tm[NEL][NX][NY][NZ];

  double M[NEL][NX][NY][NZ];
  double J[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nReading topological and model"
                                " files and computing jacobian...\n");

  if (checkTopoAndModelIO (readMeshModelAndJacobian (ic, argv[4], argv[1],
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

  double Rb[NPT][ns];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Creating radial basis...\n");

  radialBasis (r1, r2, NPT, dg, nnt, ns, T, Rb);

  double Bm[ns][ns];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing radial basis matrix...\n");

  basisMatrix (r1, r2, NPT, dg, ns, Rb, Bm);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Inverting radial basis matrix...\n");

  if (gaussJordan (ns, Bm))
  {
    fprintf (stderr, "Error: radial basis matrix is singular!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  long double nF[nlg];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing normalization factors...\n");

  nmlFactors (nmax, nlg, nF);

  double C[NPT][nmax + 1];
  double S[NPT][nmax + 1];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing azimuthal basis...\n");

  azimuthalBasis (p1, p2, NPT, nmax, C, S);

  double K[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing part of the integrals...\n");

  precomputeIntegral (W, M, J, K);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Precomputing indices...\n");

  struct BasisIndices Bi[NEL][NX][NY][NZ];

  preComputeIndices (r1, r2, t1, t2, p1, p2, Tm, Bi);

  double A[ns][nlg];
  double B[ns][nlg];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nComputing coefficients...\n\n");

  computeCoefficients (ic, t1, t2, Bi, In, ns, nmax,
                       nlg, Bm, K, Rb, nF, C, S, A, B);

  MPI_Barrier (MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\n\nWriting coefficients...\n");

    if (checkCoefficientsIO (writeCoefficients (argv[5], argv[1], zone,
                                                nmax, ns, nlg,
                                                A, B))) MPI_Abort (MPI_COMM_WORLD, 1);
  }

  double Em[NEL][NX][NY][NZ];

  initializeArray (Em);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nExpanding model onto the same"
                                " GLL mesh to evaluate the error...\n");

  expandOnGLLMesh (t1, t2, nmax, ns, nlg, Bi, In, Rb, nF, C, S, A, B, Em);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Calculating RMSD...\n");

  double rmsd = computeRMSD (In, M, Em, W, J);

  if (ic == 0)
  {
    fprintf (stderr, "RMSD: %.2E\n", rmsd);
    fprintf (stderr, "Done!\n");
  }

  MPI_Finalize ();

  return 0;
}

