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

-----------------------------------------------------------------------------------------------

 BKMNS2LL

 USAGE
   ./bin/bkmns2ll PARAMETER DEPTH RESOLUTION NMAX

 EXAMPLE
   ./bin/bkmns2ll vsv 100 0.5 10

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)
   DEPTH                  - depth in which the depth slice will be created
   RESOLUTION             - spatial distance (in degrees) between both latitude and longitude
                            grid points of the depth slice
   NMAX (optional)        - maximum degree of the spherical harmonics expansion

 DESCRIPTION
   Reads the desired model parameter, the depth, and the horizontal resolution from the command
   line, and creates a vertical cross-section of the model up to the requested spherical harmonics
   degree. In case you don't provide NMAX, all the coefficients, the routine will expand all the
   coefficients. If you want the perturbations instead of the absolute values, just add a 'd' at
   beginning of the parameter code (e.g., dvs, drho, etc). The routine writes the output to a file
   called PARAMETER_DEPTH_DS.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "structs.h"
#include "spharm.h"
#include "io.h"
#include "constants.h"

static double gMeanModel (double r, unsigned nl, struct MeanModel Mm[nl])
{
  /* Computes the geometric mean model value using
     linear interpolation */
  int n = (int) nl;
  int i = (int) (nl * (r - Mm[0].r) / (Mm[n - 1].r - Mm[0].r) + 0.5);

  if (i <= 0) return Mm[0].vgm; else if (i > n - 1) return Mm[n - 1].vgm;

  double m = (Mm[i].vgm - Mm[i - 1].vgm) / (Mm[i].r - Mm[i - 1].r);

  return Mm[i].vgm + m * (r - Mm[i].r);
}

static double block2Crust (struct SphericalBoundaries *sb,
                           double rmin, double rmax,
                           double radius, double theta, double phi,
                           unsigned np_b, unsigned nt_b, unsigned nr_b,
                           double Bm[np_b][nt_b][nr_b])
{
  /* Computes parameter value from block model */
  int np = (int) np_b;
  int nt = (int) nt_b;
  int nr = (int) nr_b;

  if (phi    < sb->pmin) phi    = sb->pmin;
  if (phi    > sb->pmax) phi    = sb->pmax;
  if (theta  < sb->tmin) theta  = sb->tmin;
  if (theta  > sb->tmax) theta  = sb->tmax;
  if (radius < rmin)     radius = rmin;
  if (radius > rmax)     radius = rmax;

  double p = phi + PI;
  double t = theta;
  double r = radius - rmin;

  double dp = (sb->pmax - sb->pmin) / (np - 1);
  double dt = (sb->tmax - sb->tmin) / (nt - 1);
  double dr = (rmax - rmin) / (nr - 1);

  int i0 = (int) (p / dp);
  int j0 = (int) (t / dt);
  int k0 = (int) (r / dr);

  double p0 = i0 * dp;
  double t0 = j0 * dt;
  double r0 = k0 * dr;

  int i1 = (i0 < np - 1) ? i0 + 1 : i0;
  int j1 = (j0 < nt - 1) ? j0 + 1 : j0;
  int k1 = (k0 < nr - 1) ? k0 + 1 : k0;

  double p1 = p0 + dp;
  double t1 = t0 + dt;
  double r1 = r0 + dr;

  double f000 = Bm[i0][j0][k0];
  double f100 = Bm[i1][j0][k0];
  double f010 = Bm[i0][j1][k0];
  double f110 = Bm[i1][j1][k0];
  double f001 = Bm[i0][j0][k1];
  double f101 = Bm[i1][j0][k1];
  double f011 = Bm[i0][j1][k1];
  double f111 = Bm[i1][j1][k1];

  double sum = 0; unsigned n = 0;

  if (f000 >= WATER_LEVEL) {sum += f000; n++;}
  if (f100 >= WATER_LEVEL) {sum += f100; n++;}
  if (f010 >= WATER_LEVEL) {sum += f010; n++;}
  if (f110 >= WATER_LEVEL) {sum += f110; n++;}
  if (f001 >= WATER_LEVEL) {sum += f001; n++;}
  if (f101 >= WATER_LEVEL) {sum += f101; n++;}
  if (f011 >= WATER_LEVEL) {sum += f011; n++;}
  if (f111 >= WATER_LEVEL) {sum += f111; n++;}

  if (n == 0) return 0;

  if (n < 8)
  {
    double mean = sum / n;

    if (f000 < WATER_LEVEL) f000 = mean;
    if (f100 < WATER_LEVEL) f100 = mean;
    if (f010 < WATER_LEVEL) f010 = mean;
    if (f110 < WATER_LEVEL) f110 = mean;
    if (f001 < WATER_LEVEL) f001 = mean;
    if (f101 < WATER_LEVEL) f101 = mean;
    if (f011 < WATER_LEVEL) f011 = mean;
    if (f111 < WATER_LEVEL) f111 = mean;
  }

  double k000 = f000 * (p1 - p) * (t1 - t) * (r1 - r);
  double k100 = f100 * (p - p0) * (t1 - t) * (r1 - r);
  double k010 = f010 * (p1 - p) * (t - t0) * (r1 - r);
  double k110 = f110 * (p - p0) * (t - t0) * (r1 - r);
  double k001 = f001 * (p1 - p) * (t1 - t) * (r - r0);
  double k101 = f101 * (p - p0) * (t1 - t) * (r - r0);
  double k011 = f011 * (p1 - p) * (t - t0) * (r - r0);
  double k111 = f111 * (p - p0) * (t - t0) * (r - r0);

  return (k000 + k100 + k010 + k110 +
          k001 + k101 + k011 + k111) / (dp * dt * dr);
}

static void expandCrust (struct SphericalBoundaries *sb,
                         double r1, double r2, double r,
                         unsigned np, unsigned nt,
                         double Theta[nt], double Phi[np],
                         unsigned np_b, unsigned nt_b, unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b], double M[np][nt],
                         unsigned nl,
                         struct MeanModel Mm[nl], bool dvv)
{
  /* Expands crust from block model */
  if (r >= r1 && r <= r2)

    for (unsigned i = 0; i < nt; i++)

      for (unsigned j = 0; j < np; j++)
      {
        double v = block2Crust (sb, r1, r2, r, Theta[i], Phi[j],
                                np_b, nt_b, nr_b, Bm);

        M[j][i] = dvv && v ? 100 * log (v / gMeanModel (r, nl, Mm)) : v;
      }
}

static void expandMantel (double r1, double r2, double r,
                          unsigned nt, unsigned np,
                          double Theta[nt], double Phi[np],
                          unsigned ns, unsigned dg,
                          unsigned nnt, unsigned N,
                          unsigned nlg, double T[nnt],
                          double A[ns][nlg], double B[ns][nlg],
                          double M[np][nt], unsigned nl,
                          struct MeanModel Mm[nl], bool dvv)
{
  /* Expands mantle from B-splines and spherical harmonics coefficients */
  if (r >= r1 && r <= r2)
  {
    double dr_1 = (NPT - 1) / (r2 - r1);
    unsigned ri = v2Index (r, r1, dr_1);

    double (*Rb)[ns] = malloc (sizeof (double[NPT][ns]));

    radialBasis (r1, r2, NPT, dg, nnt, ns, T, Rb);

    long double *nF = malloc (sizeof (long double[nlg]));

    nmlFactors (N, nlg, nF);

    long double  (*P)[N + 2] = malloc (sizeof (long double[nt][N + 2]));
    long double (*nP)[N + 2] = malloc (sizeof (long double[nt][N + 2]));

    double (*C)[N + 1] = malloc (sizeof (double[np][N + 1]));
    double (*S)[N + 1] = malloc (sizeof (double[np][N + 1]));

    azimuthalBasis (np, Phi, N, C, S);

    for (unsigned n = 0; n <= N; n++)
    {
      polarBasis (nt, Theta, n, nlg, N, nF, P, nP);

      for (unsigned m = 0; m <= n; m++)

        for (unsigned ti = 0; ti < nt; ti++)

          for (unsigned pi = 0; pi < np; pi++)
          {
            double Cmn = nP[ti][m + 1] * C[pi][m];
            double Smn = nP[ti][m + 1] * S[pi][m];

            unsigned mni = mN2I (m, n);

            for (unsigned s = 0; s < ns; s++)

              if (Rb[ri][s]) M[pi][ti] += Rb[ri][s] * (Cmn * A[s][mni]
                                                    +  Smn * B[s][mni]);
          }
    }

    if (dvv)
    {
      for (unsigned ti = 0; ti < nt; ti++)

        for (unsigned pi = 0; pi < np; pi++)

          M[pi][ti] = 100 * log (M[pi][ti] / gMeanModel (r, nl, Mm));
    }

    free (Rb);
    free (nF);
    free (P);
    free (nP);
    free (C);
    free (S);
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n BKMNS2LL"

                    "\n\n USAGE"
                    "\n    ./bin/bkmns2ll PARAMETER DEPTH RESOLUTION NMAX"

                    "\n\n EXAMPLE"
                    "\n    ./bin/bkmns2ll vsv 100 0.5 10"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)"
                    "\n    DEPTH                  - depth in which the depth slice will be created"
                    "\n    RESOLUTION             - spatial distance (in degrees) between both latitude and longitude"
                    "\n                             grid points of the depth slice"
                    "\n    NMAX (optional)        - maximum degree of the spherical harmonics expansion"

                    "\n\n DESCRIPTION"
                    "\n    Reads the desired model parameter, the depth, and the horizontal resolution from the command"
                    "\n    line, and creates a vertical cross-section of the model up to the requested spherical harmonics"
                    "\n    degree. In case you don't provide NMAX, all the coefficients, the routine will expand all the"
                    "\n    coefficients. If you want the perturbations instead of the absolute values, just add a 'd' at"
                    "\n    beginning of the parameter code (e.g., dvs, drho, etc). The routine writes the output to a file"
                    "\n    called PARAMETER_DEPTH_DS.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 4 && argc != 5)
  {
    fprintf (stderr, "\n Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  char *prm = argv[1];

  bool dvv = prm[0] == 'd' ? true : false;

  if (dvv) prm = &prm[1];

  unsigned nmax = UINT_MAX;

  if (argc == 5) nmax = atoi (argv[4]);

  double r  = depth2R (atof (argv[2]));
  double dh = atof (argv[3]);

  unsigned np = (unsigned) (360.0 / dh + 1);
  unsigned nt = (unsigned) (180.0 / dh + 1);

  struct SphericalBoundaries sb;

  sb.rmin =  CMB_R;
  sb.rmax =  TOP_R;
  sb.tmin =  0;
  sb.tmax =  PI;
  sb.pmin = -PI;
  sb.pmax =  PI;

  if (r < sb.rmin || r > sb.rmax)
  {
    fprintf (stderr, "\n Error: depth should be between %.1lf and %.1lf km....\n",
             r2Depth (sb.rmin), r2Depth (sb.rmax)); return 1;
  }

  fprintf (stderr, "\nReading mean model...\n");

  unsigned nl = 0;

  if (checkMeanModelHeaderIO (readMeanModelHeader (prm, &nl))) return 1;

  struct MeanModel Mm[nl];

  if (checkMeanModelIO (readMeanModel (prm, nl, Mm))) return 1;

  double r1min = MOHO_R;
  double r1max = TOP_R;

  unsigned np_b, nt_b, nr_b;

  if (checkBlockModelHeaderIO (readBlockModelHeader (prm,
                                                     &np_b, &nt_b, &nr_b))) return 1;

  double (*Bm)[nt_b][nr_b] = malloc (sizeof (double[np_b][nt_b][nr_b]));

  fprintf (stderr, "\nReading block model from %s.bin file...\n", prm);

  if (checkBlockModelIO (readBlockModel (prm, np_b, nt_b, nr_b, Bm))) return 1;

  fprintf (stderr, "\nCreating output grid...\n");

  double Theta[nt], Phi[np];

  createOutputGridLL (nt, np, Theta, Phi);

  double (*M)[nt] = malloc (sizeof (double[np][nt]));

  initialize2DArray (np, nt, M);

  fprintf (stderr, "Zone 1 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r1max), r2Depth (r1min));
  fprintf (stderr, "Expanding zone 1...");

  expandCrust (&sb, r1min, r1max, r, np, nt, Theta, Phi, np_b, nt_b, nr_b, Bm, M, nl, Mm, dvv);

  free (Bm);

  unsigned nn2;
  unsigned dg2;

  double r2min;
  double r2max;

  fprintf (stderr, "\n\nReading B-splines information from 'knots_Z2.dat' file...\n");

  if (checkBSplinesHeaderIO (readBSplinesHeader (&r2min, &r2max,
                                                 &nn2, &dg2, 2))) return 1;

  unsigned nnt2 = nnAndDg2Nnt (nn2, dg2);

  double T2[nnt2];

  fprintf (stderr, "Creating knots array...\n");

  if (checkKnotsIO (readKnots (r2min, r2max, dg2, 2, nnt2, T2))) return 1;

  unsigned N2;
  unsigned ns2;

  if (checkCoefficientsHeaderIO (readCoefficientsHeader (2, prm,
                                                         &N2, &ns2))) return 1;

  unsigned nlg2 = N2Nlg (N2);

  double (*A2)[nlg2] = malloc (sizeof (double[ns2][nlg2]));
  double (*B2)[nlg2] = malloc (sizeof (double[ns2][nlg2]));

  fprintf (stderr, "Reading coefficients from 'mns_Z2_%s.dat' file...\n", prm);

  if (checkCoefficientsIO (readCoefficients (2, prm, N2, ns2,
                                             nlg2, A2, B2))) return 1;

  if (N2 > nmax) N2 = nmax;

  fprintf (stderr, "Zone 2 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r2max), r2Depth (r2min));
  fprintf (stderr, "Expanding zone 2...");

  expandMantel (r2min, r2max, r, nt, np, Theta, Phi, ns2, dg2,
                nnt2, N2, nlg2, T2, A2, B2, M, nl, Mm, dvv);

  free (A2);
  free (B2);

  unsigned nn3;
  unsigned dg3;

  double r3min;
  double r3max;

  fprintf (stderr, "\n\nReading B-splines information from 'knots_Z3.dat' file...\n");

  if (checkBSplinesHeaderIO (readBSplinesHeader (&r3min, &r3max,
                                                 &nn3, &dg3, 3))) return 1;

  unsigned nnt3 = nnAndDg2Nnt (nn3, dg3);

  double T3[nnt3];

  fprintf (stderr, "Creating knots array...\n");

  if (checkKnotsIO (readKnots (r3min, r3max, dg3, 3, nnt3, T3))) return 1;

  unsigned N3;
  unsigned ns3;

  if (checkCoefficientsHeaderIO (readCoefficientsHeader (3, prm,
                                                         &N3, &ns3))) return 1;

  unsigned nlg3 = N2Nlg (N3);

  double (*A3)[nlg3] = malloc (sizeof (double[ns3][nlg3]));
  double (*B3)[nlg3] = malloc (sizeof (double[ns3][nlg3]));

  fprintf (stderr, "Reading coefficients from 'mns_Z3_%s.dat' file...\n", prm);

  if (checkCoefficientsIO (readCoefficients (3, prm, N3, ns3,
                                             nlg3, A3, B3))) return 1;

  if (N3 > nmax) N3 = nmax;

  fprintf (stderr, "Zone 3 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r3max), r2Depth (r3min));
  fprintf (stderr, "Expanding zone 3...");

  expandMantel (r3min, r3max, r, nt, np, Theta, Phi, ns3, dg3,
                nnt3, N3, nlg3, T3, A3, B3, M, nl, Mm, dvv);

  free (A3);
  free (B3);

  unsigned nn4;
  unsigned dg4;

  double r4min;
  double r4max;

  fprintf (stderr, "\n\nReading B-splines information from 'knots_Z4.dat' file...\n");

  if (checkBSplinesHeaderIO (readBSplinesHeader (&r4min, &r4max,
                                                 &nn4, &dg4, 4))) return 1;

  unsigned nnt4 = nnAndDg2Nnt (nn4, dg4);

  double T4[nnt4];

  fprintf (stderr, "Creating knots array...\n");

  if (checkKnotsIO (readKnots (r4min, r4max, dg4, 4, nnt4, T4))) return 1;

  unsigned N4;
  unsigned ns4;

  if (checkCoefficientsHeaderIO (readCoefficientsHeader (4, prm,
                                                         &N4, &ns4))) return 1;

  unsigned nlg4 = N2Nlg (N4);

  double (*A4)[nlg4] = malloc (sizeof (double[ns4][nlg4]));
  double (*B4)[nlg4] = malloc (sizeof (double[ns4][nlg4]));

  fprintf (stderr, "Reading coefficients from 'mns_Z4_%s.dat' file...\n", prm);

  if (checkCoefficientsIO (readCoefficients (4, prm, N4, ns4,
                                             nlg4, A4, B4))) return 1;

  if (N4 > nmax) N4 = nmax;

  fprintf (stderr, "Zone 4 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r4max), r2Depth (r4min));
  fprintf (stderr, "Expanding zone 4...");

  expandMantel (r4min, r4max, r, nt, np, Theta, Phi, ns4, dg4,
                nnt4, N4, nlg4, T4, A4, B4, M, nl, Mm, dvv);

  free (A4);
  free (B4);

  fprintf (stderr, "\n\nWriting expansion...\n");

  if (checkExpansionIO (writeExpansionLL (prm, dvv, &sb, r, np, nt, M))) return 1;

  free (M);

  fprintf (stderr, "Done!\n");

  return 0;
}

