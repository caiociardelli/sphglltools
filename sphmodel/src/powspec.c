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

 POWSPEC

 USAGE
   ./bin/powspec PARAMETER NMAX

 EXAMPLE
   ./bin/powspec vsv 60

 COMMAND-LINE ARGUMENTS
   PARAMETER              - model parameter for which the power spectrum will be created
                            (vph, rho, eta, vsv, etc.)
   NMAX                   - maximum degree of the spherical harmonics expansion

 DESCRIPTION
   Reads the parameter and the maximum degree of the expansion from the command line and creates
   both the one- and the two-dimensional power spectra. The routine writes two output files called
   PARAMETER_pwspc.dat and  PARAMETER_pwspc2D.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "spharm.h"
#include "io.h"
#include "constants.h"

static void powerSpectrum1D (unsigned ns, unsigned nmax, unsigned nlg,
                             double Bm[ns][ns],
                             double A[ns][nlg], double B[ns][nlg],
                             double Pws1D[nmax + 1])
{
  /* Computes 1D power spectrum */
  for (unsigned n = 0; n <= nmax; n++)
  {
    double nf = 1.0 / (2 * n + 1);

    Pws1D[n] = 0.0;

    for (unsigned s1 = 0; s1 < ns; s1++)

      for (unsigned s2 = 0; s2 < ns; s2++)
      {
        double sum = 0.0;

        for (unsigned m = 0; m <= n; m++)
        {
          unsigned mni = mN2I (m, n);

          sum += A[s1][mni] * A[s2][mni] +
                 B[s1][mni] * B[s2][mni];
        }

        Pws1D[n] += Bm[s1][s2] * sum;
      }

    Pws1D[n] *= nf;
  }
}

static void powerSpectrum2D (unsigned Nr, unsigned ns, unsigned nmax,
                             unsigned nlg, double R[Nr][ns],
                             double A[ns][nlg], double B[ns][nlg],
                             double Pws2D[Nr][nmax + 1])
{
  /* Computes 2D power spectrum */
  for (unsigned n = 0; n <= nmax; n++)
  {
    double nf = 1.0 / (2 * n + 1);

    for (unsigned i = 0; i < Nr; i++)
    {
      Pws2D[i][n] = 0.0;

      for (unsigned s1 = 0; s1 < ns; s1++)

        for (unsigned s2 = 0; s2 < ns; s2++)
        {
          double sum = 0.0;

          for (unsigned m = 0; m <= n; m++)
          {
            unsigned mni = mN2I (m, n);

            sum += A[s1][mni] * A[s2][mni] +
                   B[s1][mni] * B[s2][mni];
          }

          Pws2D[i][n] += R[i][s1] * R[i][s2] * sum;
        }

      Pws2D[i][n] *= nf;
    }
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n POWSPEC"

                    "\n\n USAGE"
                    "\n    ./bin/powspec PARAMETER NMAX"

                    "\n\n EXAMPLE"
                    "\n    ./bin/powspec vsv 60"

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter for which the power spectrum will be created"
                    "\n                             (vph, rho, eta, vsv, etc.)"
                    "\n    NMAX                   - maximum degree of the spherical harmonics expansion"

                    "\n\n DESCRIPTION"
                    "\n    Reads the parameter and the maximum degree of the expansion from the command line and creates"
                    "\n    both the one- and the two-dimensional power spectra. The routine writes two output files called"
                    "\n    PARAMETER_pwspc.dat and  PARAMETER_pwspc2D.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 3)
  {
    fprintf (stderr, "\n Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  char *prm = argv[1];
  unsigned nmax = atoi (argv[2]);

  fprintf (stderr, "\nStarting power spectrum calculation from Zone 2 to 4...\n");

  unsigned nn2;
  unsigned dg2;

  double r2min;
  double r2max;

  fprintf (stderr, "\nReading B-splines information from 'knots_Z2.dat' file...\n");

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

  fprintf (stderr, "Zone 2 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r2max), r2Depth (r2min));

  double (*Rb2)[ns2] = malloc (sizeof (double[NPT][ns2]));

  fprintf (stderr, "Creating radial basis...\n");

  radialBasis (r2min, r2max, NPT, dg2, nnt2, ns2, T2, Rb2);

  fprintf (stderr, "Computing power spectrum for Z2...\n");

  double Pws2D2[NPT][nmax + 1];

  powerSpectrum2D (NPT, ns2, nmax, nlg2, Rb2, A2, B2, Pws2D2);

  unsigned nn3;
  unsigned dg3;

  double r3min;
  double r3max;

  fprintf (stderr, "\nReading B-splines information from 'knots_Z3.dat' file...\n");

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

  fprintf (stderr, "Zone 3 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r3max), r2Depth (r3min));

  double (*Rb3)[ns3] = malloc (sizeof (double[NPT][ns3]));

  fprintf (stderr, "Creating radial basis...\n");

  radialBasis (r3min, r3max, NPT, dg3, nnt3, ns3, T3, Rb3);

  fprintf (stderr, "Computing power spectrum for Z3...\n");

  double Pws2D3[NPT][nmax + 1];

  powerSpectrum2D (NPT, ns3, nmax, nlg3, Rb3, A3, B3, Pws2D3);

  unsigned nn4;
  unsigned dg4;

  double r4min;
  double r4max;

  fprintf (stderr, "\nReading B-splines information from 'knots_Z4.dat' file...\n");

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

  fprintf (stderr, "Zone 4 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (r4max), r2Depth (r4min));

  double (*Rb4)[ns4] = malloc (sizeof (double[NPT][ns4]));

  fprintf (stderr, "Creating radial basis...\n");

  radialBasis (r4min, r4max, NPT, dg4, nnt4, ns4, T4, Rb4);

  fprintf (stderr, "Computing power spectrum for Z4...\n");

  double Pws2D4[NPT][nmax + 1];

  powerSpectrum2D (NPT, ns4, nmax, nlg4, Rb4, A4, B4, Pws2D4);

  double Bm2[ns2][ns2];

  fprintf (stderr, "\nCreating basis matrix for Z2...\n");

  basisMatrix (r2min, r2max, NPT, dg2, ns2, Rb2, Bm2);

  free (Rb2);

  double Bm3[ns3][ns3];

  fprintf (stderr, "Creating basis matrix for Z3...\n");

  basisMatrix (r3min, r3max, NPT, dg3, ns3, Rb3, Bm3);

  free (Rb3);

  double Bm4[ns4][ns4];

  fprintf (stderr, "Creating basis matrix for Z4...\n");

  basisMatrix (r4min, r4max, NPT, dg4, ns4, Rb4, Bm4);

  free (Rb4);

  fprintf (stderr, "\nComputing 1D power spectrum for the whole mantle...\n");

  double Pws1D2[nmax + 1];
  double Pws1D3[nmax + 1];
  double Pws1D4[nmax + 1];

  powerSpectrum1D (ns2, nmax, nlg2, Bm2, A2, B2, Pws1D2);
  powerSpectrum1D (ns3, nmax, nlg3, Bm3, A3, B3, Pws1D3);
  powerSpectrum1D (ns4, nmax, nlg4, Bm4, A4, B4, Pws1D4);

  free (A2); free (B2);
  free (A3); free (B3);
  free (A4); free (B4);

  fprintf (stderr, "\nWriting 1D power spectrum...");

  if (checkPowSpecIO (writePowSpec1D (prm, nmax,
                                      Pws1D2, Pws1D3, Pws1D3))) return 1;

  fprintf (stderr, "\nWriting 2D power spectrum...\n");

  if (checkPowSpecIO (writePowSpec2D (prm, nmax, NPT,
                                      Pws2D2, Pws2D3, Pws2D4))) return 1;

  fprintf (stderr, "Done!\n");

  return 0;
}

