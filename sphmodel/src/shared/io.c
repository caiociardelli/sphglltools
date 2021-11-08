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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"
#include "coordinates.h"
#include "constants.h"

unsigned readMeanModelHeader (char *prm, unsigned *nl)
{
  /* Reads the header of the mean model file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "1D_mean/%s.dat", prm) < 14) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  int f1 = fscanf (file, "%*s");
  int f2 = fscanf (file, "%u", nl);

  fclose (file);

  if (f1)      return 3;
  if (f2 != 1) return 4;

  return 0;
}

unsigned readMeanModel (char *prm, unsigned nl,
                        struct MeanModel Mm[nl])
{
  /* Reads the mean model file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "1D_mean/%s.dat", prm) < 14) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  if (fscanf (file, "%*[^\n]\n") != 0) return 3;
  if (fscanf (file, "%*[^\n]\n") != 0) return 3;

  double depth, vam, vgm;

  for (unsigned l = 0; l < nl; l++)
  {
    if (fscanf (file, "%lf %lf %lf", &depth, &vam, &vgm) != 3) return 3;

    Mm[nl - l - 1].r   = depth2R (depth);
    Mm[nl - l - 1].vam = vam;
    Mm[nl - l - 1].vgm = vgm;
  }

  fclose (file);

  return 0;
}

static inline char machineEndianess (void)
{
  /* Returns the machine endianness */
  int c = 1; return (*(char*) &c == 1) ? '<' : '>';
}

static unsigned uintSwap (unsigned input)
{
  /* Swaps byte order for unsigned
     integer numbers */
  int output;

  char *in  = (char*) &input;
  char *out = (char*) &output;

  out[0] = in[3];
  out[1] = in[2];
  out[2] = in[1];
  out[3] = in[0];

  return output;
}

static float floatSwap (float input)
{
  /* Swaps byte order for floating-point
     numbers */
  float output;

  char *in  = (char*) &input;
  char *out = (char*) &output;

  out[0] = in[3];
  out[1] = in[2];
  out[2] = in[1];
  out[3] = in[0];

  return output;
}

unsigned readBlockModelHeader (char *prm,
                               unsigned *np_b,
                               unsigned *nt_b,
                               unsigned *nr_b)
{
  /* Reads the crustal block model header */
  char name[MAX_STRING_LEN];

  if (sprintf (name, "crust/%s.bin", prm) < 12) return 1;

  FILE *input = fopen (name, "rb");

  if (input == NULL) return 2;

  char e;

  if (fread (&e, sizeof (char), 1, input) != 1) return 3;

  if (fread (np_b, sizeof (unsigned), 1, input) != 1) return 4;
  if (fread (nt_b, sizeof (unsigned), 1, input) != 1) return 4;
  if (fread (nr_b, sizeof (unsigned), 1, input) != 1) return 4;

  if (e != machineEndianess ())
  {
    *np_b = uintSwap (*np_b);
    *nt_b = uintSwap (*nt_b);
    *nr_b = uintSwap (*nr_b);
  }

  return 0;
}

unsigned readBlockModel (char *prm,
                         unsigned np_b,
                         unsigned nt_b,
                         unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b])
{
  /* Reads the crustal block model */
  unsigned n = np_b * nt_b * nr_b;

  char name[MAX_STRING_LEN];

  if (sprintf (name, "crust/%s.bin", prm) < 12) return 1;

  FILE *input = fopen (name, "rb");

  if (input == NULL) return 2;

  char e;

  if (fread (&e, sizeof (char), 1, input) != 1) return 3;

  unsigned junk;

  if (fread (&junk, sizeof (unsigned), 1, input) != 1) return 4;
  if (fread (&junk, sizeof (unsigned), 1, input) != 1) return 4;
  if (fread (&junk, sizeof (unsigned), 1, input) != 1) return 4;

  float (*tBm)[nt_b][nr_b] = malloc (sizeof (float[np_b][nt_b][nr_b]));

  if (fread (tBm, sizeof (float), n, input) != n) return 5;

  for (unsigned i = 0; i < np_b; i++)

    for (unsigned j = 0; j < nt_b; j++)

      for (unsigned k = 0; k < nr_b; k++)

        Bm[i][j][k] = (e == machineEndianess ()) ? tBm[i][j][k] : floatSwap (tBm[i][j][k]);

  free (tBm);

  return 0;
}

unsigned readBSplinesHeader (double *rmin, double *rmax,
                             unsigned *nn, unsigned *dg,
                             unsigned zone)
{
  /* Reads the header of the B-splines file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "mantle/knots_Z%u.dat", zone) < 19) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  int f1 = fscanf (file, "%*[^\n]\n");
  int f2 = fscanf (file, "%u", nn);
  int f3 = fscanf (file, "%u", dg);
  int f4 = fscanf (file, "%lf", rmin);
  int f5 = fscanf (file, "%lf", rmax);

  fclose (file);

  if (f1)      return 3;
  if (f2 != 1) return 4;
  if (f3 != 1) return 5;
  if (f4 != 1) return 6;
  if (f5 != 1) return 7;

  return 0;
}

unsigned readKnots (double r1, double r2, unsigned dg,
                    unsigned zone, unsigned nnt, double T[nnt])
{
  /* Reads the B-splines knots */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "mantle/knots_Z%u.dat", zone) < 19) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  for (unsigned i = 0; i <= dg; i++)
  {
    T[i] = r1;
    T[nnt - i - 1] = r2;
  }

  int f1 = fscanf (file, "%*[^\n]\n");
  int f2 = fscanf (file, "%*[^\n]\n");
  int f3 = fscanf (file, "%*[^\n]\n");
  int f4 = fscanf (file, "%*[^\n]\n");
  int f5 = fscanf (file, "%*[^\n]\n");

  if (f1 || f2 || f3 || f4 || f5)
  {
    fclose (file); return 3;
  }

  for (unsigned i = dg + 1; i < nnt - dg - 1; i++)

    if (fscanf (file, "%lf", &T[i]) != 1)
    {
      fclose (file); return 4;
    }

  fclose (file);

  return 0;
}

unsigned readCoefficientsHeader (unsigned zone, char *prm,
                                 unsigned *N, unsigned *ns)
{
  /* Reads coefficients file header */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "mantle/mns_Z%u_%s.dat", zone, prm) < 20) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  char string[MAX_STRING_LEN];

  for (unsigned i = 0; i < 21; i++)
  {
    int f = fscanf (file, "%s", string);

    if (f != 1)
    {
      fclose (file);

      return 3;
    }

    if (i == 17) *N  = atoi (string);
    if (i == 20) *ns = atoi (string);
  }

  fclose (file);

  return 0;
}

unsigned readCoefficients (unsigned zone, char *prm, unsigned N,
                           unsigned ns, unsigned nlg,
                           double A[ns][nlg], double B[ns][nlg])
{
  /* Reads coefficients file header */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "mantle/mns_Z%u_%s.dat", zone, prm) < 20) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  if (fscanf (file, "%*[^\n]\n"))
  {
    fclose (file);

    return 3;
  }

  char junk1[MAX_STRING_LEN];
  char junk2[MAX_STRING_LEN];

  for (unsigned s = 0; s < ns; s++)

    for (unsigned n = 0; n <= N; n++)

      for (unsigned m = 0; m <= n; m++)
      {
        unsigned mni = mN2I (m, n);

        int fA = fscanf (file, "%s %s %lf", junk1, junk2, &A[s][mni]);
        int fB = fscanf (file, "%s %s %lf", junk1, junk2, &B[s][mni]);

        if (fA != 3 || fB != 3)
        {
          fclose (file);

          return 4;
        }
      }

  fclose (file);

  return 0;
}

unsigned readNumberOfPoints (unsigned *np, char *filename)
{
  /* Reads number of points in the
     input file */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  for (int c = getc (file); c != EOF; c = getc (file))

    if (c == '\n') *np += 1;

  fclose (file);

  return 0;
}

unsigned readCoordinates (unsigned np, char *filename,
                          double *rmin, double *rmax,
                          double R[np], double Theta[np], double Phi[np])
{
  /* Reads the coordinates of the points */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  unsigned l = 0;

  double latitude, longitude, depth;

  for (unsigned i = 0; i < np; i++)
  {
    l += fscanf (file, "%lf %lf %lf", &latitude, &longitude, &depth);

    R[i]     = depth2R (depth);
    Theta[i] = degree2Rad (90 - latitude);
    Phi[i]   = degree2Rad (longitude);

    if (R[i] < *rmin) *rmin = R[i];
    if (R[i] > *rmax) *rmax = R[i];
  }

  fclose (file);

  if (l != 3 * np) return 2;

  return 0;
}

void createProfilePf (double r1, double r2,
                      unsigned nr, double R[nr])
{
  /* Creates grid for 1D profile */
  double dr = (r2 - r1) / (nr - 1);

  for (unsigned i = 0; i < nr; i++)

    R[i] = r1 + i * dr;
}

void createOutputGridLL (unsigned nt, unsigned np,
                         double Theta[nt], double Phi[np])
{
  /* Creates grid for depth slice */
  double dt =     PI / (nt - 1);
  double dp = 2 * PI / (np - 1);

  double t = 0;

  for (unsigned i = 0; i < nt; i++)
  {
    Theta[i] = t; t += dt;
  }

  double p = -PI;

  for (unsigned i = 0; i < np; i++)
  {
    Phi[i] = p; p += dp;
  }
}

static void crossProduct (double x1, double y1, double z1,
                          double x2, double y2, double z2,
                          double *u, double *v, double *w)
{
  /* Computes cross product between pair of points in 3D */
  *u = y1 * z2 - z1 * y2;
  *v = z1 * x2 - x1 * z2;
  *w = x1 * y2 - y1 * x2;

  double norm = sqrt (square (*u) + square (*v) + square (*w));

  *u = *u / norm; *v = *v / norm; *w = *w / norm;
}

void createOutputGridDD (double r1, double r2,
                         double t1, double t2,
                         double p1, double p2,
                         unsigned nr, unsigned nd,
                         double R[nr], double Delta[nd],
                         double Theta[nd], double Phi[nd])
{
  /* Creates grid for vertical cross section */
  double dr = (r2 - r1) / (nr - 1);

  for (unsigned i = 0; i < nr; i++)

    R[i] = r1 + i * dr;

  double dd = vincenty (t1, p1, t2, p2) / (nd - 1);

  double x1, y1, z1;
  double x2, y2, z2;

  rThetaPhi2XYZ (1, t1, p1, &x1, &y1, &z1);
  rThetaPhi2XYZ (1, t2, p2, &x2, &y2, &z2);

  double u, v, w;

  crossProduct (x1, y1, z1, x2, y2, z2, &u, &v, &w);

  for (unsigned i = 0; i < nd; i++)
  {
    Delta[i] = rad2Degree (i * dd);

    double sind = sin (i * dd);
    double cosd = cos (i * dd);

    double x = x1 * (cosd + square (u) * (1 - cosd))
             + y1 * (u * v * (1 - cosd) - w * sind)
             + z1 * (u * w * (1 - cosd) + v * sind);

    double y = x1 * (v * u * (1 - cosd) + w * sind)
             + y1 * (cosd + square (v) * (1 - cosd))
             + z1 * (v * w * (1 - cosd) - u * sind);

    double z = x1 * (w * u * (1 - cosd) - v * sind)
             + y1 * (w * v * (1 - cosd) + u * sind)
             + z1 * (cosd + square (w) * (1 - cosd));

    double r;

    xYZ2RThetaPhi (x, y, z, &r, &Theta[i], &Phi[i]);
  }
}

unsigned writeExpansionLL (char *prm, bool dvv,
                           struct SphericalBoundaries *sb,
                           double r,
                           unsigned np, unsigned nt,
                           double M[np][nt])
{
  /* Writes depth slice */
  char name[MAX_STRING_LEN];

  double depth = r2Depth (r);

  if (dvv)
  {
    if (sprintf (name, "d%s_%g_DS.dat", prm, depth) < 12) return 1;
  }

  else
  {
    if (sprintf (name, "%s_%g_DS.dat", prm, depth) < 11) return 1;
  }

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %g %u %u\n", depth, nt, np);

  if (dvv)
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)      dVpVp\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dVpvVpv\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dVphVph\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)      dVsVs\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dVsvVsv\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dVshVsh\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dEtaEta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)     dRhoRho\n");
  }

  else
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vp (km/s)\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vpv (km/s)\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vph (km/s)\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vs (km/s)\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vsv (km/s)\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Vsh (km/s)\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)       Eta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#latitude (degrees)   longitude (degrees)    Rho (g/cm^3)\n");
  }

  double dp = (sb->pmax - sb->pmin) / (np - 1);
  double dt = (sb->tmax - sb->tmin) / (nt - 1);

  double t = sb->tmax;

  for (int ti = nt - 1; ti >= 0; ti--)
  {
    double p = sb->pmin;

    for (unsigned pi = 0; pi < np; pi++)
    {
      fprintf (file, "%15.7lf %21.7lf %18E\n", 90 - rad2Degree (t),
                                               rad2Degree (p),
                                               M[pi][ti]);

      p += dp;
    }

    t -= dt;
  }

  fclose (file);

  return 0;
}

unsigned writeExpansionDD (char *argv[], char *prm, bool dvv,
                           unsigned nr, unsigned nd,
                           double R[nr],
                           double Delta[nd],
                           double M[nr][nd])
{
  /* Writes vertical cross section */
  char name[MAX_STRING_LEN];

  if (dvv)
  {
    if (sprintf (name, "d%s_VCS.dat", prm) < 11) return 1;
  }

  else
  {
    if (sprintf (name, "%s_VCS.dat", prm) < 10) return 1;
  }

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n", argv[4], argv[5], argv[6], argv[7]);

  if (dvv)
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#radius (km)   delta (degrees)      dVpVp\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dVpvVpv\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dVphVph\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#radius (km)   delta (degrees)      dVsVs\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dVsvVsv\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dVshVsh\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dEtaEta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#radius (km)   delta (degrees)     dRhoRho\n");
  }

  else
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vp (km/s)\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vpv (km/s)\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vph (km/s)\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vs (km/s)\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vsv (km/s)\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Vsh (km/s)\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#radius (km)   delta (degrees)       Eta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#radius (km)   delta (degrees)    Rho (g/cm^3)\n");
  }

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3lf %14.3lf %19E\n", R[i] * EARTH_R,
                                               Delta[j], M[i][j]);

  return 0;
}

unsigned writeExpansionPf (char *argv[], bool dvv,
                           unsigned nr,
                           double R[nr], double M[nr])
{
  /* Writes 1D profile */
  char *prm = argv[1];
  char name[MAX_STRING_LEN];

  if (dvv)
  {
    if (sprintf (name, "%s_PF.dat", prm) < 10) return 1;
  }

  else
  {
    if (sprintf (name, "%s_PF.dat", prm) < 9) return 1;
  }

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[4], argv[5]);

  if (dvv)
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#depth (km)       dVpVp\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#depth (km)      dVpvVpv\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#depth (km)      dVphVph\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#depth (km)       dVsVs\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#depth (km)      dVsvVsv\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#depth (km)      dVshVsh\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#depth (km)      dEtaEta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#depth (km)      dRhoRho\n");
  }

  else
  {
    if (strcmp (prm, "vp") == 0)

      fprintf (file, "#depth (km)     Vp (km/s)\n");

    else if (strcmp (prm, "vpv") == 0)

      fprintf (file, "#depth (km)     Vpv (km/s)\n");

    else if (strcmp (prm, "vph") == 0)

      fprintf (file, "#depth (km)     Vph (km/s)\n");

    else if (strcmp (prm, "vs") == 0)

      fprintf (file, "#depth (km)     Vs (km/s)\n");

    else if (strcmp (prm, "vsv") == 0)

      fprintf (file, "#depth (km)     Vsv (km/s)\n");

    else if (strcmp (prm, "vsh") == 0)

      fprintf (file, "#depth (km)     Vsh (km/s)\n");

    else if (strcmp (prm, "eta") == 0)

      fprintf (file, "#depth (km)        Eta\n");

    else if (strcmp (prm, "rho") == 0)

      fprintf (file, "#depth (km)     Rho (g/cm^3)\n");
  }

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3lf %17E\n", r2Depth (R[i]), M[i]);

  return 0;
}

unsigned writeModel1D (char *prm, unsigned nr,
                       double R[nr], struct MeanModel Mo[nr])
{
  /* Writes mean model */
  char name[MAX_STRING_LEN];

  FILE *file;

  if (sprintf (name, "%s.dat", prm) < 6) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);

  if (strcmp (prm, "vp") == 0)

    fprintf (file, "#depth (km)    Vp_arithmetic (km/s)    Vp_geometric (km/s)\n");

  else if (strcmp (prm, "vpv") == 0)

    fprintf (file, "#depth (km)   Vpv_arithmetic (km/s)   Vpv_geometric (km/s)\n");

  else if (strcmp (prm, "vph") == 0)

    fprintf (file, "#depth (km)   Vph_arithmetic (km/s)   Vph_geometric (km/s)\n");

  else if (strcmp (prm, "vs") == 0)

    fprintf (file, "#depth (km)    Vs_arithmetic (km/s)    Vs_geometric (km/s)\n");

  else if (strcmp (prm, "vsv") == 0)

    fprintf (file, "#depth (km)   Vsv_arithmetic (km/s)   Vsv_geometric (km/s)\n");

  else if (strcmp (prm, "vsh") == 0)

    fprintf (file, "#depth (km)   Vsh_arithmetic (km/s)   Vsh_geometric (km/s)\n");

  else if (strcmp (prm, "eta") == 0)

    fprintf (file, "#depth (km)      Eta_arithmetic          Eta_geometric\n");

  else if (strcmp (prm, "rho") == 0)

    fprintf (file, "#depth (km)   Rho_arithmetic (g/cm^3)   Rho_geometric (g/cm^3)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1lf %20.4E %23.4E\n", r2Depth (R[i]),
                                             Mo[i].vam, Mo[i].vgm);

  fclose (file);

  return 0;
}

void writeExpansionPath (unsigned np, double M[np])
{
  /* Writes values along the path to the standard output */
  for (unsigned i = 0; i < np; i++)

    fprintf (stdout, "%12E\n", M[i]);
}

unsigned writePowSpec1D (char *prm, unsigned nmax,
                         double Pws1D2[nmax + 1],
                         double Pws1D3[nmax + 1],
                         double Pws1D4[nmax + 1])
{
  /* Writes power spectra */
  char name[MAX_STRING_LEN];

  if (sprintf (name, "%s_pwspc.dat", prm) < 12) return 1;

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#Nmax: %u\n", nmax);
  fprintf (file, "#n  Pws (%s)\n", prm);

  for (unsigned n = 0; n <= nmax; n++)

    fprintf (file, "%u %E\n", n, Pws1D2[n] +
                                 Pws1D3[n] +
                                 Pws1D4[n]);

  fclose (file);

  return 0;
}

static inline unsigned mapIndex (unsigned i, unsigned nr, unsigned Nr)
{
  /* Maps local index to global */
  return (unsigned) (((double) (nr - i) / nr) * (Nr - 1) + 0.5);
}

unsigned writePowSpec2D (char *prm, unsigned nmax,
                         unsigned Nr,
                         double Pws2D2[Nr][nmax + 1],
                         double Pws2D3[Nr][nmax + 1],
                         double Pws2D4[Nr][nmax + 1])
{
  /* Writes power spectra per depth */
  char name[MAX_STRING_LEN];

  if (sprintf (name, "%s_pwspc2D.dat", prm) < 14) return 1;

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndepth Nmax: %u %u\n", 2812, nmax);
  fprintf (file, "#depth (km) ");

  for (unsigned n = 0; n <= nmax; n++)

    fprintf (file, "N%u ", n);

  fprintf (file, "\n");

  for (unsigned i = 0; i < 330; i++)
  {
    double depth = 80 + i;
    unsigned index = mapIndex (i, 330, Nr);

    fprintf (file, "%6.0lf   ", depth);

    for (unsigned n = 0; n <= nmax; n++)

      fprintf (file, "%E ", Pws2D2[index][n]);

    fprintf (file, "\n");
  }

  for (unsigned i = 0; i < 240; i++)
  {
    double depth = 410 + i;
    unsigned index = mapIndex (i, 240, Nr);

    fprintf (file, "%g ", depth);

    for (unsigned n = 0; n <= nmax; n++)

      fprintf (file, "%E ", Pws2D3[index][n]);

    fprintf (file, "\n");
  }

  for (unsigned i = 0; i <= 2241; i++)
  {
    double depth = 650 + i;
    unsigned index = mapIndex (i, 2241, Nr);

    fprintf (file, "%g ", depth);

    for (unsigned n = 0; n <= nmax; n++)

      fprintf (file, "%E ", Pws2D4[index][n]);

    fprintf (file, "\n");
  }

  fclose (file);

  return 0;
}

unsigned checkBlockModelHeaderIO (unsigned rvalue)
{
  /* Checks IO of the block model */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "\n Error: could not read in the endianess of the file...!\n");
    break;

    case 4:
      fprintf (stderr, "\n Error: could not read block model dimensions!\n");
    break;
  }

  return rvalue;
}

unsigned checkBlockModelIO (unsigned rvalue)
{
  /* Checks IO of the block model */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "\n Error: could not read in the endianess of the file...!\n");
    break;

    case 4:
      fprintf (stderr, "\n Error: could not read block model dimensions!\n");
    break;

    case 5:
      fprintf (stderr, "\n Error: could not read block model file!\n");
    break;
  }

  return rvalue;
}

unsigned checkBSplinesHeaderIO (unsigned rvalue)
{
  /* Checks IO of the B-splines file header */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file header!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read number of knots!\n");
    break;

    case 5:
      fprintf (stderr, "Error: could not read b-splines degree!\n");
    break;

    case 6:
      fprintf (stderr, "Error: could not read minimum radius!\n");
    break;

    case 7:
      fprintf (stderr, "Error: could not read maximum radius!\n");
    break;
  }

  return rvalue;
}

unsigned checkKnotsIO (unsigned rvalue)
{
  /* Checks IO of the B-splines file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read header information!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read knots positions!\n");
    break;
  }

  return rvalue;
}

unsigned checkMeanModelHeaderIO (unsigned rvalue)
{
  /* Checks IO of the mean model file header */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file header!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read number of points!\n");
    break;
  }

  return rvalue;
}

unsigned checkMeanModelIO (unsigned rvalue)
{
  /* Checks IO of the mean model file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not open file!\n");
    break;

    case 3:
      fprintf (stderr, "\n Error: could not read file!\n");
    break;
  }

  return rvalue;
}

unsigned checkCoordinatesIO (unsigned rvalue)
{
  /* Checks IO of the file containing the coordinates */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\nError opening the file...\n");
    break;

    case 2:
      fprintf (stderr, "\nError reading the file...\n");
    break;
  }

  return rvalue;
}


unsigned checkCoefficientsHeaderIO (unsigned rvalue)
{
  /* Checks IO of the coefficients file header */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "\n Error: could not read header information!\n");
    break;
  }

  return rvalue;
}

unsigned checkCoefficientsIO (unsigned rvalue)
{
  /* Checks IO of the coefficients file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "\n Error: could not read header information!\n");
    break;

    case 4:
      fprintf (stderr, "\n Error: could not read coefficients!\n");
    break;
  }

  return rvalue;
}

unsigned checkExpansionIO (unsigned rvalue)
{
  /* Checks IO of the output file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not open file!\n");
    break;
  }

  return rvalue;
}

unsigned checkPowSpecIO (unsigned rvalue)
{
  /* Checks IO of the output file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\n Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "\n Error: could not open file!\n");
    break;
  }

  return rvalue;
}

void helpMenuLL (void)
{
  /* Prints the help menu */
  char *help_menu = "\n BKMNS2LL"

                    "\n\n USAGE"
                    "\n    ./bkmns2ll PARAMETER DEPTH RESOLUTION NMAX"

                    "\n\n EXAMPLE"
                    "\n    ./bkmns2ll vs 10 0.5 40"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vp, vs, rho, vpv, vph, vsv, vsh or eta)"
                    "\n    DEPTH                  - depth in which the depth slice will be created"
                    "\n    RESOLUTION             - spatial distance (in degrees) between both latitude and longitude"
                    "\n                             grid points of the depth slice"
                    "\n    NMAX (optional)        - maximum degree of the spherical harmonics expansion"

                    "\n\n DESCRIPTION"
                    "\n    Creates a depth slice of the model up to the requested spherical harmonics degree given the"
                    "\n    desired model parameter, depth and horizontal resolution. In case NMAX is not provided, all"
                    "\n    the coefficients are expanded. If you want the perturbations instead of the absolute values,"
                    "\n    just add a 'd' at beginning of the parameter code (e.g dvs, drho, etc). The output is writen"
                    "\n    to a file called PARAMETER_DEPTH_DS.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

void helpMenuDD (void)
{
  char *help_menu = "\n BKMNS2DD"

                    "\n\n USAGE"
                    "\n    ./bkmns2dd PARAMETER MIN_RADIUS MAX_RADIUS LAT1 LON1 LAT2 LON2 RADIAL_RESOLUTION HORIZONTAL_RESOLUTION NMAX"

                    "\n\n EXAMPLE"
                    "\n    ./bkmns2dd dvs 80 2891 -1.4 -25.1 23.7 51.3 2.0 0.1 40"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vp, vs, rho, vpv, vph, vsv, vsh or eta)"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    LAT1, LON1             - latitude and longitude of the firts point"
                    "\n    LAT2, LON2             - latitude and longitude of the second point"
                    "\n    RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial direction"
                    "\n    HORIZONTAL_RESOLUTION  - spatial distance (in degrees) between grid points along the horizontal direction"
                    "\n    NMAX (optional)        - maximum degree of the spherical harmonics expansion"

                    "\n\n DESCRIPTION"
                    "\n    Creates a vertical cross section of the model up to the requested spherical harmonics degree"
                    "\n    given the desired model parameter, minimum depth, maximum depth, coordinates of the two points"
                    "\n    defining the great circle path and the spatial resolutions along the radial direction and the"
                    "\n    great circle. In case NMAX is not provided, all the coefficients are expanded. If you want the"
                    "\n    perturbations instead of the absolute values, just add a 'd' at beginning of the parameter code"
                    "\n    (e.g dvs, drho, etc). The output is writen to a file called PARAMETER_VCS.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

void helpMenuPf (void)
{
  char *help_menu = "\n BKMNS2PF"

                    "\n\n USAGE"
                    "\n    ./bkmns2pf PARAMETER MIN_DEPTH MAX_DEPTH LAT LON RADIAL_RESOLUTION NMAX"

                    "\n\n EXAMPLE"
                    "\n    ./bkmns2pf vs 10 2891 0.8 -24.2 1.0 40"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vp, vs, rho, vpv, vph, vsv, vsh or eta)"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    LAT LON                - latitude and longitude of the point"
                    "\n    RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial direction"
                    "\n    NMAX (optional)        - maximum degree of the spherical harmonics expansion"

                    "\n\n DESCRIPTION"
                    "\n    Creates a 1D vertical profile of the model up to the requested spherical harmonics degree given"
                    "\n    the desired model parameter, minimum depth, maximum depth and the coordinates of the point in"
                    "\n    which you want the profile. In case NMAX is not provided, all the coefficients are expanded."
                    "\n    If you want the perturbations instead of the absolute values, just add a 'd' at beginning of"
                    "\n    the parameter code (e.g dvs, drho, etc). The output is writen to a file called"
                    "\n    PARAMETER_PF.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

