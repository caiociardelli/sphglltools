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

 GLL2MEAN

 USAGE
   mpiexec -n 36 bin/gll2mean MIN_DEPTH MAX_DEPTH HORIZONTAL_RESOLUTION RADIAL_RESOLUTION
                              INPUT_DIRECTORY OUTPUT_DIRECTORY

 EXAMPLE
   mpiexec -n 36 bin/gll2mean -5 2891 0.5 1.0 data/INPUT/ .

 COMMAND-LINE ARGUMENTS
   MIN_DEPTH              - minimum depth
   MAX_DEPTH              - maximum depth
   HORIZONTAL_RESOLUTION  - spatial distance (in degrees) between grid points along the
                            horizontal direction
   RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial
                            direction
   INPUT_DIRECTORY        - directory containing the input files
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads the minimum and maximum depths, the horizontal and vertical grid spacings, and the
   input and output directory names from the command line and computes the mean models for the
   desired parameters (defined in 'config.h'). The routine writes the output to files called
   PARAMETER.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "structs.h"
#include "metrics.h"
#include "io.h"
#include "boundaries.h"
#include "progress.h"
#include "constants.h"
#include "config.h"

static long double iDWCorrection (long double u, long double v, struct FaceNode Fn[NPf])
{
  /* Computes a correction for a face surface using IDW */
  long double s = 0.L, ws = 0.L;

  for (unsigned i = 0; i < NPf; i++)
  {
    long double k  = squaredDistance2D (u, Fn[i].u, v, Fn[i].v);
    long double wi = 1.L / (WATER_LEVEL_1 + powl (k, 1.5L));

    s += wi * Fn[i].w; ws += wi;
  }

  return s / ws;
}

static bool isInFaceShadow (struct Point *p, struct Edges *e)
{
  /* Checks if a point is inside the face shadow */
  long double e1_u = getCoordinate (p, e->e1.corners.o);
  long double e1_v = getCoordinate (p, e->e1.o);
  long double e2_u = getCoordinate (p, e->e2.corners.o);
  long double e2_v = getCoordinate (p, e->e2.o);
  long double e3_u = getCoordinate (p, e->e3.corners.o);
  long double e3_v = getCoordinate (p, e->e3.o);
  long double e4_u = getCoordinate (p, e->e4.corners.o);
  long double e4_v = getCoordinate (p, e->e4.o);

  long double c_v;

  c_v = isInEdgeShadow (e1_u, &e->e1.corners)
      ? edge (e1_u, &e->e1.inner)
      : edge (e1_u, &e->e1.outer);

  if ((e->e1.p == '>' && e1_v + e->ti < c_v) ||
      (e->e1.p == '<' && e1_v - e->ti > c_v)) return false;

  c_v = isInEdgeShadow (e2_u, &e->e2.corners)
      ? edge (e2_u, &e->e2.inner)
      : edge (e2_u, &e->e2.outer);

  if ((e->e2.p == '>' && e2_v + e->ti < c_v) ||
      (e->e2.p == '<' && e2_v - e->ti > c_v)) return false;

  c_v = isInEdgeShadow (e3_u, &e->e3.corners)
      ? edge (e3_u, &e->e3.inner)
      : edge (e3_u, &e->e3.outer);

  if ((e->e3.p == '>' && e3_v + e->tj < c_v) ||
      (e->e3.p == '<' && e3_v - e->tj > c_v)) return false;

  c_v = isInEdgeShadow (e4_u, &e->e4.corners)
      ? edge (e4_u, &e->e4.inner)
      : edge (e4_u, &e->e4.outer);

  if ((e->e4.p == '>' && e4_v + e->tj < c_v) ||
      (e->e4.p == '<' && e4_v - e->tj > c_v)) return false;

  return true;
}

static bool isInSpecElm (struct Point *p, struct Boundaries *b, struct Faces *f)
{
  /* Checks if a point is inside a spectral element */
  if (p->x     < b->xmin) return false;
  if (p->x     > b->xmax) return false;
  if (p->y     < b->ymin) return false;
  if (p->y     > b->ymax) return false;
  if (p->z     < b->zmin) return false;
  if (p->z     > b->zmax) return false;
  if (p->r     < b->rmin) return false;
  if (p->r     > b->rmax) return false;
  if (p->theta < b->tmin) return false;
  if (p->theta > b->tmax) return false;
  if (p->phi   < b->pmin) return false;
  if (p->phi   > b->pmax) return false;

  long double f1_u = 0.L, f1_v = 0.L, f1_w = 0.L;
  long double f2_u = 0.L, f2_v = 0.L, f2_w = 0.L;
  long double f3_u = 0.L, f3_v = 0.L, f3_w = 0.L;
  long double f4_u = 0.L, f4_v = 0.L, f4_w = 0.L;
  long double f5_u = 0.L, f5_v = 0.L, f5_w = 0.L;
  long double f6_u = 0.L, f6_v = 0.L, f6_w = 0.L;

  switch (f->f1.o)
  {
    case 'x':
      f1_u = p->y;
      f1_v = p->z;
      f1_w = p->x;
    break;

    case 'y':
      f1_u = p->x;
      f1_v = p->z;
      f1_w = p->y;
    break;

    case 'z':
      f1_u = p->x;
      f1_v = p->y;
      f1_w = p->z;
    break;
  }

  switch (f->f2.o)
  {
    case 'x':
      f2_u = p->y;
      f2_v = p->z;
      f2_w = p->x;
    break;

    case 'y':
      f2_u = p->x;
      f2_v = p->z;
      f2_w = p->y;
    break;

    case 'z':
      f2_u = p->x;
      f2_v = p->y;
      f2_w = p->z;
    break;
  }

  switch (f->f3.o)
  {
    case 'x':
      f3_u = p->y;
      f3_v = p->z;
      f3_w = p->x;
    break;

    case 'y':
      f3_u = p->x;
      f3_v = p->z;
      f3_w = p->y;
    break;

    case 'z':
      f3_u = p->x;
      f3_v = p->y;
      f3_w = p->z;
    break;
  }

  switch (f->f4.o)
  {
    case 'x':
      f4_u = p->y;
      f4_v = p->z;
      f4_w = p->x;
    break;

    case 'y':
      f4_u = p->x;
      f4_v = p->z;
      f4_w = p->y;
    break;

    case 'z':
      f4_u = p->x;
      f4_v = p->y;
      f4_w = p->z;
    break;
  }

  switch (f->f5.o)
  {
    case 'x':
      f5_u = p->y;
      f5_v = p->z;
      f5_w = p->x;
    break;

    case 'y':
      f5_u = p->x;
      f5_v = p->z;
      f5_w = p->y;
    break;

    case 'z':
      f5_u = p->x;
      f5_v = p->y;
      f5_w = p->z;
    break;
  }

  switch (f->f6.o)
  {
    case 'x':
      f6_u = p->y;
      f6_v = p->z;
      f6_w = p->x;
    break;

    case 'y':
      f6_u = p->x;
      f6_v = p->z;
      f6_w = p->y;
    break;

    case 'z':
      f6_u = p->x;
      f6_v = p->y;
      f6_w = p->z;
    break;
  }

  long double tlr, s_w;

  tlr = f->ti + f->f1.sigma;

  s_w = isInFaceShadow (p, &f->f1.edges)
      ? face (f1_u, f1_v, &f->f1.inner) + iDWCorrection (f1_u, f1_v, f->f1.Fn)
      : face (f1_u, f1_v, &f->f1.outer);

  if ((f->f1.p == '>' && f1_w + tlr < s_w) ||
      (f->f1.p == '<' && f1_w - tlr > s_w)) return false;

  tlr = f->ti + f->f2.sigma;

  s_w = isInFaceShadow (p, &f->f2.edges)
      ? face (f2_u, f2_v, &f->f2.inner) + iDWCorrection (f2_u, f2_v, f->f2.Fn)
      : face (f2_u, f2_v, &f->f2.outer);

  if ((f->f2.p == '>' && f2_w + tlr < s_w) ||
      (f->f2.p == '<' && f2_w - tlr > s_w)) return false;

  tlr = f->tj + f->f3.sigma;

  s_w = isInFaceShadow (p, &f->f3.edges)
      ? face (f3_u, f3_v, &f->f3.inner) + iDWCorrection (f3_u, f3_v, f->f3.Fn)
      : face (f3_u, f3_v, &f->f3.outer);

  if ((f->f3.p == '>' && f3_w + tlr < s_w) ||
      (f->f3.p == '<' && f3_w - tlr > s_w)) return false;

  tlr = f->tj + f->f4.sigma;

  s_w = isInFaceShadow (p, &f->f4.edges)
      ? face (f4_u, f4_v, &f->f4.inner) + iDWCorrection (f4_u, f4_v, f->f4.Fn)
      : face (f4_u, f4_v, &f->f4.outer);

  if ((f->f4.p == '>' && f4_w + tlr < s_w) ||
      (f->f4.p == '<' && f4_w - tlr > s_w)) return false;

  tlr = f->tk + f->f5.sigma;

  s_w = isInFaceShadow (p, &f->f5.edges)
      ? face (f5_u, f5_v, &f->f5.inner) + iDWCorrection (f5_u, f5_v, f->f5.Fn)
      : face (f5_u, f5_v, &f->f5.outer);

  if ((f->f5.p == '>' && f5_w + tlr < s_w) ||
      (f->f5.p == '<' && f5_w - tlr > s_w)) return false;

  tlr = f->tk + f->f6.sigma;

  s_w = isInFaceShadow (p, &f->f6.edges)
      ? face (f6_u, f6_v, &f->f6.inner) + iDWCorrection (f6_u, f6_v, f->f6.Fn)
      : face (f6_u, f6_v, &f->f6.outer);

  if ((f->f6.p == '>' && f6_w + tlr < s_w) ||
      (f->f6.p == '<' && f6_w - tlr > s_w)) return false;

  return true;
}

static bool findSpecElm (struct Point *p, unsigned nel,
                         struct Point Ti[nel][NX][NY][NZ],
                         struct Boundaries b[nel], struct Faces f[nel],
                         unsigned *element)
{
  /* Finds the spectral element that contains a given point */
  unsigned hnx = NX / 2, hny = NY / 2, hnz = NZ / 2;
  long double mnd = INFINITY; bool found = false;

  for (unsigned el = 0; el < nel; el++)

    if (isInSpecElm (p, &b[el], &f[el]))
    {
      long double d = squaredDistance3D (p, &Ti[el][hnx][hny][hnz]);

      if (d < mnd)
      {
        *element = el;

        mnd   = d;
        found = true;
      }
    }

  return found;
}

static void findClosestPoint (struct Point *p, unsigned nel, unsigned el,
                              struct Point T[nel][NX][NY][NZ],
                              unsigned *ii, unsigned *jj, unsigned *kk)
{
  /* Finds the closest point inside the spectral element */
  long double mnd = INFINITY;

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)
      {
        long double d = squaredDistance3D (p, &T[el][i][j][k]);

        if (d < mnd)
        {
          *ii = i; *jj = j; *kk = k; mnd = d;
        }
      }
}

static unsigned findFirstGuesses (int ic, unsigned nel,
                                  unsigned np_l, unsigned nt, unsigned nr,
                                  struct Boundaries b[nel],
                                  struct Faces f[nel],
                                  struct Point Ti[nel][NX][NY][NZ],
                                  long double R[nr],
                                  long double Theta[nt],
                                  long double Phi[np_l],
                                  struct Indexes Eijk[np_l][nt][nr])
{
  /* Finds initial guesses for the interpolation points */
  clock_t starttime = clock ();

  for (unsigned k = 0; k < nr; k++)
  {
    long double r = R[k];

    for (unsigned j = 0; j < nt; j++)
    {
      long double theta = Theta[j];

      for (unsigned i = 0; i < np_l; i++)
      {
        long double phi = Phi[i], x, y, z;

        rThetaPhi2XYZl (r, theta, phi, &x, &y, &z);

        struct Point p;

        p.r     = r;
        p.theta = theta;
        p.phi   = phi;
        p.x     = x;
        p.y     = y;
        p.z     = z;

        struct Indexes eijk;

        if (findSpecElm (&p, nel, Ti, b, f, &eijk.element))
        {
          findClosestPoint (&p, nel, eijk.element, Ti, &eijk.i, &eijk.j, &eijk.k);

          eijk.outside_mesh = false;

          Eijk[i][j][k] = eijk;
        }

        else if (p.r > MIN_SURFACE_R) Eijk[i][j][k].outside_mesh = true;

        else return 1;
      }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (ic == 0) progressBar (k, 1, nr, starttime);
  }

  return 0;
}

static unsigned computeShape (unsigned nel, struct CartesianPoint spf[nel][NN],
                              struct Point T[nel][NX][NY][NZ])
{
  /* Computes shape functions for all spectrals elements */
  long double xi[]    = {-1.L,  0.L,  1.L, -1.L,  0.L,  1.L, -1.L,  0.L,  1.L,
                         -1.L,  0.L,  1.L, -1.L,  0.L,  1.L, -1.L,  0.L,  1.L,
                         -1.L,  0.L,  1.L, -1.L,  0.L,  1.L, -1.L,  0.L,  1.L};
  long double eta[]   = {-1.L, -1.L, -1.L,  0.L,  0.L,  0.L,  1.L,  1.L,  1.L,
                         -1.L, -1.L, -1.L,  0.L,  0.L,  0.L,  1.L,  1.L,  1.L,
                         -1.L, -1.L, -1.L,  0.L,  0.L,  0.L,  1.L,  1.L,  1.L};
  long double gamma[] = {-1.L, -1.L, -1.L, -1.L, -1.L, -1.L, -1.L, -1.L, -1.L,
                          0.L,  0.L,  0.L,  0.L,  0.L,  0.L,  0.L,  0.L,  0.L,
                          1.L,  1.L,  1.L,  1.L,  1.L,  1.L,  1.L,  1.L,  1.L};

  long double M[NN][NN];

  for (unsigned r = 0; r < NN; r++)
  {
    unsigned c = 0;

    for (unsigned k = 0; k < NNe; k++)

      for (unsigned j = 0; j < NNe; j++)

        for (unsigned i = 0; i < NNe; i++)

          M[r][c++] = cNNN (xi[r], eta[r], gamma[r], i, j, k);
  }

  if (gaussJordanl (NN, M)) return 1;

  unsigned ni[] = {0, 2, 4}; /* Assumes that NX = NY = NZ = 5 */

  for (unsigned el = 0; el < nel; el++)
  {
    unsigned r = 0;

    long double vx[NN], vy[NN], vz[NN];

    for (unsigned k = 0; k < NNe; k++)

      for (unsigned j = 0; j < NNe; j++)

        for (unsigned i = 0; i < NNe; i++)
        {
          struct Point p = T[el][ni[i]][ni[j]][ni[k]];

          vx[r] = p.x;
          vy[r] = p.y;
          vz[r] = p.z;

          r++;
        }

    for (unsigned r = 0; r < NN; r++)
    {
      spf[el][r].x = 0; spf[el][r].y = 0; spf[el][r].z = 0;

      for (unsigned c = 0; c < NN; c++)
      {
        spf[el][r].x += M[r][c] * vx[c];
        spf[el][r].y += M[r][c] * vy[c];
        spf[el][r].z += M[r][c] * vz[c];
      }
    }
  }

  return 0;
}

static void shape (long double xi, long double eta, long double gamma,
                   unsigned nel, unsigned el,
                   struct CartesianPoint spf[nel][NN],
                   struct CartesianPoint *p)
{
  /* Maps a point from the standard base to the new base */
  unsigned r = 0;

  p->x = 0.L; p->y = 0.L; p->z = 0.L;

  for (unsigned k = 0; k < NNe; k++)

    for (unsigned j = 0; j < NNe; j++)

      for (unsigned i = 0; i < NNe; i++)
      {
        long double c = cNNN (xi, eta, gamma, i, j, k);

        p->x += spf[el][r].x * c;
        p->y += spf[el][r].y * c;
        p->z += spf[el][r].z * c;

        r++;
      }
}

static void jacobian (long double xi, long double eta, long double gamma,
                      unsigned nel, unsigned el,
                      struct CartesianPoint spf[nel][NN],
                      long double J[3][3])
{
  /* Computes the Jacobian matrix */
  unsigned r = 0;

  J[0][0] = 0.L; J[0][1] = 0.L; J[0][2] = 0.L;
  J[1][0] = 0.L; J[1][1] = 0.L; J[1][2] = 0.L;
  J[2][0] = 0.L; J[2][1] = 0.L; J[2][2] = 0.L;

  for (unsigned k = 0; k < 3; k++)

    for (unsigned j = 0; j < 3; j++)

      for (unsigned i = 0; i < 3; i++)
      {
        long double spx = spf[el][r].x;
        long double spy = spf[el][r].y;
        long double spz = spf[el][r].z;

        if (i > 0)
        {
          long double ci = ((long double) i) * cNNN (xi, eta, gamma, i - 1, j, k);

          J[0][0] += ci * spx;
          J[1][0] += ci * spy;
          J[2][0] += ci * spz;
        }

        if (j > 0)
        {
          long double cj = ((long double) j) * cNNN (xi, eta, gamma, i, j - 1, k);

          J[0][1] += cj * spx;
          J[1][1] += cj * spy;
          J[2][1] += cj * spz;
        }

        if (k > 0)
        {
          long double ck = ((long double) k) * cNNN (xi, eta, gamma, i, j, k - 1);

          J[0][2] += ck * spx;
          J[1][2] += ck * spy;
          J[2][2] += ck * spz;
        }

        r++;
      }
}

static unsigned newton (long double *xi, long double *eta, long double *gamma,
                        unsigned nel, unsigned el,
                        struct Point *tp,
                        struct CartesianPoint spf[nel][NN])
{
  /* Unmaps the coordinates back to the standard base using Newton's method */
  long double xin    = *xi;
  long double etan   = *eta;
  long double gamman = *gamma;

  for (unsigned n = 0; n < MAX_NITER; n++)
  {
    struct CartesianPoint p;

    shape (xin, etan, gamman, nel, el, spf, &p);

    long double xn = p.x - tp->x;
    long double yn = p.y - tp->y;
    long double zn = p.z - tp->z;

    if (norm (xn, yn, zn) < EPSILON)
    {
      *xi = xin; *eta = etan; *gamma = gamman; break;
    }

    long double J[3][3];

    jacobian (xin, etan, gamman, nel, el, spf, J);

    if (gaussJordanl (3, J)) return 1;

    xin    -= J[0][0] * xn + J[0][1] * yn + J[0][2] * zn;
    etan   -= J[1][0] * xn + J[1][1] * yn + J[1][2] * zn;
    gamman -= J[2][0] * xn + J[2][1] * yn + J[2][2] * zn;
  }

  return 0;
}

static unsigned computeInterpolator (unsigned nel, struct Parameters a[nel][NX * NY * NZ],
                                     long double X[NX], long double Y[NY], long double Z[NZ],
                                     struct Parameters M[nel][NX][NY][NZ])
{
  /* Computes the interpolator for all spectral elements */
  unsigned N = NX * NY * NZ;

  unsigned r = 0;

  long double xi[N], eta[N], gamma[N];

  for (unsigned k = 0; k < NZ; k++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned i = 0; i < NX; i++)
      {
        xi[r]    = X[i];
        eta[r]   = Y[j];
        gamma[r] = Z[k];

        r++;
      }

  long double C[N][N];

  for (unsigned r = 0; r < N; r++)
  {
    unsigned c = 0;

    for (unsigned k = 0; k < NZ; k++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned i = 0; i < NX; i++)

          C[r][c++] = cNNN (xi[r], eta[r], gamma[r], i, j, k);
  }

  if (gaussJordanl (N, C)) return 1;

  for (unsigned el = 0; el < nel; el++)
  {
    unsigned r = 0;

    struct Parameters v[N];

    for (unsigned k = 0; k < NZ; k++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned i = 0; i < NX; i++)
        {
          #if VP
          v[r].vp  = M[el][i][j][k].vp;
          #endif

          #if VS
          v[r].vs  = M[el][i][j][k].vs;
          #endif

          #if RHO
          v[r].rho = M[el][i][j][k].rho;
          #endif

          #if VPV
          v[r].vpv = M[el][i][j][k].vpv;
          #endif

          #if VPH
          v[r].vph = M[el][i][j][k].vph;
          #endif

          #if VSV
          v[r].vsv = M[el][i][j][k].vsv;
          #endif

          #if VSH
          v[r].vsh = M[el][i][j][k].vsh;
          #endif

          #if ETA
          v[r].eta = M[el][i][j][k].eta;
          #endif

          #if QMU
          v[r].qmu = M[el][i][j][k].qmu;
          #endif

          #if GCP
          v[r].gcp = M[el][i][j][k].gcp;
          #endif

          #if GSP
          v[r].gsp = M[el][i][j][k].gsp;
          #endif

          #if MU0
          v[r].mu0 = M[el][i][j][k].mu0;
          #endif

          #if APK
          v[r].apk = M[el][i][j][k].apk;
          #endif

          #if BTK
          v[r].btk = M[el][i][j][k].btk;
          #endif

          #if RHK
          v[r].rhk = M[el][i][j][k].rhk;
          #endif

          #if BCK
          v[r].bck = M[el][i][j][k].bck;
          #endif

          #if BBK
          v[r].bbk = M[el][i][j][k].bbk;
          #endif

          #if BVK
          v[r].bvk = M[el][i][j][k].bvk;
          #endif

          #if BHK
          v[r].bhk = M[el][i][j][k].bhk;
          #endif

          #if ETK
          v[r].etk = M[el][i][j][k].etk;
          #endif

          #if GCK
          v[r].gck = M[el][i][j][k].gck;
          #endif

          #if GSK
          v[r].gsk = M[el][i][j][k].gsk;
          #endif

          #if HSK
          v[r].hsk = M[el][i][j][k].hsk;
          #endif

          r++;
        }

    for (unsigned r = 0; r < N; r++)
    {
      #if VP
      a[el][r].vp  = 0.L;
      #endif

      #if VS
      a[el][r].vs  = 0.L;
      #endif

      #if RHO
      a[el][r].rho = 0.L;
      #endif

      #if VPV
      a[el][r].vpv = 0.L;
      #endif

      #if VPH
      a[el][r].vph = 0.L;
      #endif

      #if VSV
      a[el][r].vsv = 0.L;
      #endif

      #if VSH
      a[el][r].vsh = 0.L;
      #endif

      #if ETA
      a[el][r].eta = 0.L;
      #endif

      #if QMU
      a[el][r].qmu = 0.L;
      #endif

      #if GCP
      a[el][r].gcp = 0.L;
      #endif

      #if GSP
      a[el][r].gsp = 0.L;
      #endif

      #if MU0
      a[el][r].mu0 = 0.L;
      #endif

      #if APK
      a[el][r].apk = 0.L;
      #endif

      #if BTK
      a[el][r].btk = 0.L;
      #endif

      #if RHK
      a[el][r].rhk = 0.L;
      #endif

      #if BCK
      a[el][r].bck = 0.L;
      #endif

      #if BBK
      a[el][r].bbk = 0.L;
      #endif

      #if BVK
      a[el][r].bvk = 0.L;
      #endif

      #if BHK
      a[el][r].bhk = 0.L;
      #endif

      #if ETK
      a[el][r].etk = 0.L;
      #endif

      #if GCK
      a[el][r].gck = 0.L;
      #endif

      #if GSK
      a[el][r].gsk = 0.L;
      #endif

      #if HSK
      a[el][r].hsk = 0.L;
      #endif

      for (unsigned c = 0; c < N; c++)
      {
        #if VP
        a[el][r].vp  += C[r][c] * v[c].vp;
        #endif

        #if VS
        a[el][r].vs  += C[r][c] * v[c].vs;
        #endif

        #if RHO
        a[el][r].rho += C[r][c] * v[c].rho;
        #endif

        #if VPV
        a[el][r].vpv += C[r][c] * v[c].vpv;
        #endif

        #if VPH
        a[el][r].vph += C[r][c] * v[c].vph;
        #endif

        #if VSV
        a[el][r].vsv += C[r][c] * v[c].vsv;
        #endif

        #if VSH
        a[el][r].vsh += C[r][c] * v[c].vsh;
        #endif

        #if ETA
        a[el][r].eta += C[r][c] * v[c].eta;
        #endif

        #if QMU
        a[el][r].qmu += C[r][c] * v[c].qmu;
        #endif

        #if GCP
        a[el][r].gcp += C[r][c] * v[c].gcp;
        #endif

        #if GSP
        a[el][r].gsp += C[r][c] * v[c].gsp;
        #endif

        #if MU0
        a[el][r].mu0 += C[r][c] * v[c].mu0;
        #endif

        #if APK
        a[el][r].apk += C[r][c] * v[c].apk;
        #endif

        #if BTK
        a[el][r].btk += C[r][c] * v[c].btk;
        #endif

        #if RHK
        a[el][r].rhk += C[r][c] * v[c].rhk;
        #endif

        #if BCK
        a[el][r].bck += C[r][c] * v[c].bck;
        #endif

        #if BBK
        a[el][r].bbk += C[r][c] * v[c].bbk;
        #endif

        #if BVK
        a[el][r].bvk += C[r][c] * v[c].bvk;
        #endif

        #if BHK
        a[el][r].bhk += C[r][c] * v[c].bhk;
        #endif

        #if ETK
        a[el][r].etk += C[r][c] * v[c].etk;
        #endif

        #if GCK
        a[el][r].gck += C[r][c] * v[c].gck;
        #endif

        #if GSK
        a[el][r].gsk += C[r][c] * v[c].gsk;
        #endif

        #if HSK
        a[el][r].hsk += C[r][c] * v[c].hsk;
        #endif
      }
    }
  }

  return 0;
}

static void surface (struct Parameters *m, unsigned nel, unsigned el,
                     long double xi, long double eta, long double gamma,
                     struct Parameters a[nel][NX * NY * NZ])
{
  /* Returns interpolated value */
  unsigned r = 0;

  #if VP
  m->vp  = 0.L;
  #endif

  #if VS
  m->vs  = 0.L;
  #endif

  #if RHO
  m->rho = 0.L;
  #endif

  #if VPV
  m->vpv = 0.L;
  #endif

  #if VPH
  m->vph = 0.L;
  #endif

  #if VSV
  m->vsv = 0.L;
  #endif

  #if VSH
  m->vsh = 0.L;
  #endif

  #if ETA
  m->eta = 0.L;
  #endif

  #if QMU
  m->qmu = 0.L;
  #endif

  #if GCP
  m->gcp = 0.L;
  #endif

  #if GSP
  m->gsp = 0.L;
  #endif

  #if MU0
  m->mu0 = 0.L;
  #endif

  #if APK
  m->apk = 0.L;
  #endif

  #if BTK
  m->btk = 0.L;
  #endif

  #if RHK
  m->rhk = 0.L;
  #endif

  #if BCK
  m->bck = 0.L;
  #endif

  #if BBK
  m->bbk = 0.L;
  #endif

  #if BVK
  m->bvk = 0.L;
  #endif

  #if BHK
  m->bhk = 0.L;
  #endif

  #if ETK
  m->etk = 0.L;
  #endif

  #if GCK
  m->gck = 0.L;
  #endif

  #if GSK
  m->gsk = 0.L;
  #endif

  #if HSK
  m->hsk = 0.L;
  #endif

  for (unsigned k = 0; k < NZ; k++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned i = 0; i < NX; i++)
      {
        long double cn = cNNN (xi, eta, gamma, i, j, k);

        #if VP
        m->vp  += a[el][r].vp  * cn;
        #endif

        #if VS
        m->vs  += a[el][r].vs  * cn;
        #endif

        #if RHO
        m->rho += a[el][r].rho * cn;
        #endif

        #if VPV
        m->vpv += a[el][r].vpv * cn;
        #endif

        #if VPH
        m->vph += a[el][r].vph * cn;
        #endif

        #if VSV
        m->vsv += a[el][r].vsv * cn;
        #endif

        #if VSH
        m->vsh += a[el][r].vsh * cn;
        #endif

        #if ETA
        m->eta += a[el][r].eta * cn;
        #endif

        #if QMU
        m->qmu += a[el][r].qmu * cn;
        #endif

        #if GCP
        m->gcp += a[el][r].gcp * cn;
        #endif

        #if GSP
        m->gsp += a[el][r].gsp * cn;
        #endif

        #if MU0
        m->mu0 += a[el][r].mu0 * cn;
        #endif

        #if APK
        m->apk += a[el][r].apk * cn;
        #endif

        #if BTK
        m->btk += a[el][r].btk * cn;
        #endif

        #if RHK
        m->rhk += a[el][r].rhk * cn;
        #endif

        #if BCK
        m->bck += a[el][r].bck * cn;
        #endif

        #if BBK
        m->bbk += a[el][r].bbk * cn;
        #endif

        #if BVK
        m->bvk += a[el][r].bvk * cn;
        #endif

        #if BHK
        m->bhk += a[el][r].bhk * cn;
        #endif

        #if ETK
        m->etk += a[el][r].etk * cn;
        #endif

        #if GCK
        m->gck += a[el][r].gck * cn;
        #endif

        #if GSK
        m->gsk += a[el][r].gsk * cn;
        #endif

        #if HSK
        m->hsk += a[el][r].hsk * cn;
        #endif

        r++;
      }
}

static void computeMeans (struct Means *mo,
                          unsigned np, unsigned nt,
                          struct Parameters M[np][nt])
{
  /* Computes arithmetic and geometric means */
  #if VP
  unsigned nvp = 0;
  mo->a.vp = 0.L;
  mo->g.vp = 0.L;
  #endif

  #if VS
  unsigned nvs = 0;
  mo->a.vs = 0.L;
  mo->g.vs = 0.L;
  #endif

  #if RHO
  unsigned nrho = 0;
  mo->a.rho = 0.L;
  mo->g.rho = 0.L;
  #endif

  #if VPV
  unsigned nvpv = 0;
  mo->a.vpv = 0.L;
  mo->g.vpv = 0.L;
  #endif

  #if VPH
  unsigned nvph = 0;
  mo->a.vph = 0.L;
  mo->g.vph = 0.L;
  #endif

  #if VSV
  unsigned nvsv = 0;
  mo->a.vsv = 0.L;
  mo->g.vsv = 0.L;
  #endif

  #if VSH
  unsigned nvsh = 0;
  mo->a.vsh = 0.L;
  mo->g.vsh = 0.L;
  #endif

  #if ETA
  unsigned neta = 0;
  mo->a.eta = 0.L;
  mo->g.eta = 0.L;
  #endif

  #if QMU
  unsigned nqmu = 0;
  mo->a.qmu = 0.L;
  mo->g.qmu = 0.L;
  #endif

  #if GCP
  unsigned ngcp = 0;
  mo->a.gcp = 0.L;
  mo->g.gcp = 0.L;
  #endif

  #if GSP
  unsigned ngsp = 0;
  mo->a.gsp = 0.L;
  mo->g.gsp = 0.L;
  #endif

  #if MU0
  unsigned nmu0 = 0;
  mo->a.mu0 = 0.L;
  mo->g.mu0 = 0.L;
  #endif

  #if APK
  unsigned napk = 0;
  mo->a.apk = 0.L;
  mo->g.apk = 0.L;
  #endif

  #if BTK
  unsigned nbtk = 0;
  mo->a.btk = 0.L;
  mo->g.btk = 0.L;
  #endif

  #if RHK
  unsigned nrhk = 0;
  mo->a.rhk = 0.L;
  mo->g.rhk = 0.L;
  #endif

  #if BCK
  unsigned nbck = 0;
  mo->a.bck = 0.L;
  mo->g.bck = 0.L;
  #endif

  #if BBK
  unsigned nbbk = 0;
  mo->a.bbk = 0.L;
  mo->g.bbk = 0.L;
  #endif

  #if BVK
  unsigned nbvk = 0;
  mo->a.bvk = 0.L;
  mo->g.bvk = 0.L;
  #endif

  #if BHK
  unsigned nbhk = 0;
  mo->a.bhk = 0.L;
  mo->g.bhk = 0.L;
  #endif

  #if ETK
  unsigned netk = 0;
  mo->a.etk = 0.L;
  mo->g.etk = 0.L;
  #endif

  #if GCK
  unsigned ngck = 0;
  mo->a.gck = 0.L;
  mo->g.gck = 0.L;
  #endif

  #if GSK
  unsigned ngsk = 0;
  mo->a.gsk = 0.L;
  mo->g.gsk = 0.L;
  #endif

  #if HSK
  unsigned nhsk = 0;
  mo->a.hsk = 0.L;
  mo->g.hsk = 0.L;
  #endif

  for (unsigned i = 0; i < np; i++)

    for (unsigned j = 0; j < nt; j++)
    {
      struct Parameters m = M[i][j];

      #if VP
      if (m.vp > TOLERANCE)
      {
        mo->a.vp += m.vp;
        mo->g.vp += log (m.vp);

        nvp++;
      }
      #endif

      #if VS
      if (m.vs > TOLERANCE)
      {
        mo->a.vs += m.vs;
        mo->g.vs += log (m.vs);

        nvs++;
      }
      #endif

      #if RHO
      if (m.rho > TOLERANCE)
      {
        mo->a.rho += m.rho;
        mo->g.rho += log (m.rho);

        nrho++;
      }
      #endif

      #if VPV
      if (m.vpv > TOLERANCE)
      {
        mo->a.vpv += m.vpv;
        mo->g.vpv += log (m.vpv);

        nvpv++;
      }
      #endif

      #if VPH
      if (m.vph > TOLERANCE)
      {
        mo->a.vph += m.vph;
        mo->g.vph += log (m.vph);

        nvph++;
      }
      #endif

      #if VSV
      if (m.vsv > TOLERANCE)
      {
        mo->a.vsv += m.vsv;
        mo->g.vsv += log (m.vsv);

        nvsv++;
      }
      #endif

      #if VSH
      if (m.vsh > TOLERANCE)
      {
        mo->a.vsh += m.vsh;
        mo->g.vsh += log (m.vsh);

        nvsh++;
      }
      #endif

      #if ETA
      if (m.eta > TOLERANCE)
      {
        mo->a.eta += m.eta;
        mo->g.eta += log (m.eta);

        neta++;
      }
      #endif

      #if QMU
      if (m.qmu > TOLERANCE)
      {
        mo->a.qmu += m.rho;
        mo->g.qmu += log (m.rho);

        nqmu++;
      }
      #endif

      #if GCP
      if (m.gcp > TOLERANCE)
      {
        mo->a.gcp += m.gcp;
        mo->g.gcp += log (m.gcp);

        ngcp++;
      }
      #endif

      #if GSP
      if (m.gsp > TOLERANCE)
      {
        mo->a.gsp += m.gsp;
        mo->g.gsp += log (m.gsp);

        ngsp++;
      }
      #endif

      #if MU0
      if (m.mu0 > TOLERANCE)
      {
        mo->a.mu0 += m.mu0;
        mo->g.mu0 += log (m.mu0);

        nmu0++;
      }
      #endif

      #if APK
      if (m.apk > TOLERANCE)
      {
        mo->a.apk += m.apk;
        mo->g.apk += log (m.apk);

        napk++;
      }
      #endif

      #if BTK
      if (m.btk > TOLERANCE)
      {
        mo->a.btk += m.btk;
        mo->g.btk += log (m.btk);

        nbtk++;
      }
      #endif

      #if RHK
      if (m.rhk > TOLERANCE)
      {
        mo->a.rhk += m.rhk;
        mo->g.rhk += log (m.rhk);

        nrhk++;
      }
      #endif

      #if BCK
      if (m.bck > TOLERANCE)
      {
        mo->a.bck += m.bck;
        mo->g.bck += log (m.bck);

        nbck++;
      }
      #endif

      #if BBK
      if (m.bbk > TOLERANCE)
      {
        mo->a.bbk += m.bbk;
        mo->g.bbk += log (m.bbk);

        nbbk++;
      }
      #endif

      #if BVK
      if (m.bvk > TOLERANCE)
      {
        mo->a.bvk += m.bvk;
        mo->g.bvk += log (m.bvk);

        nbvk++;
      }
      #endif

      #if BHK
      if (m.bhk > TOLERANCE)
      {
        mo->a.bhk += m.bhk;
        mo->g.bhk += log (m.bhk);

        nbhk++;
      }
      #endif

      #if ETK
      if (m.etk > TOLERANCE)
      {
        mo->a.etk += m.etk;
        mo->g.etk += log (m.etk);

        netk++;
      }
      #endif

      #if GCK
      if (m.gck > TOLERANCE)
      {
        mo->a.gck += m.gck;
        mo->g.gck += log (m.gck);

        ngck++;
      }
      #endif

      #if GSK
      if (m.gsk > TOLERANCE)
      {
        mo->a.gsk += m.gsk;
        mo->g.gsk += log (m.gsk);

        ngsk++;
      }
      #endif

      #if HSK
      if (m.hsk > TOLERANCE)
      {
        mo->a.hsk += m.hsk;
        mo->g.hsk += log (m.hsk);

        nhsk++;
      }
      #endif
    }

  #if VP
  if (nvp)
  {
    mo->a.vp /= nvp;
    mo->g.vp /= nvp;
    mo->g.vp = exp (mo->g.vp);
  }
  #endif

  #if VS
  if (nvs)
  {
    mo->a.vs /= nvs;
    mo->g.vs /= nvs;
    mo->g.vs = exp (mo->g.vs);
  }
  #endif

  #if RHO
  if (nrho)
  {
    mo->a.rho /= nrho;
    mo->g.rho /= nrho;
    mo->g.rho = exp (mo->g.rho);
  }
  #endif

  #if VPV
  if (nvpv)
  {
    mo->a.vpv /= nvpv;
    mo->g.vpv /= nvpv;
    mo->g.vpv = exp (mo->g.vpv);
  }
  #endif

  #if VPH
  if (nvph)
  {
    mo->a.vph /= nvph;
    mo->g.vph /= nvph;
    mo->g.vph = exp (mo->g.vph);
  }
  #endif

  #if VSV
  if (nvsv)
  {
    mo->a.vsv /= nvsv;
    mo->g.vsv /= nvsv;
    mo->g.vsv = exp (mo->g.vsv);
  }
  #endif

  #if VSH
  if (nvsh)
  {
    mo->a.vsh /= nvsh;
    mo->g.vsh /= nvsh;
    mo->g.vsh = exp (mo->g.vsh);
  }
  #endif

  #if ETA
  if (neta)
  {
    mo->a.eta /= neta;
    mo->g.eta /= neta;
    mo->g.eta = exp (mo->g.eta);
  }
  #endif

  #if QMU
  if (nqmu)
  {
    mo->a.qmu /= nqmu;
    mo->g.qmu /= nqmu;
    mo->g.qmu = exp (mo->g.qmu);
  }
  #endif

  #if GCP
  if (ngcp)
  {
    mo->a.gcp /= ngcp;
    mo->g.gcp /= ngcp;
    mo->g.gcp = exp (mo->g.gcp);
  }
  #endif

  #if GSP
  if (ngsp)
  {
    mo->a.gsp /= ngsp;
    mo->g.gsp /= ngsp;
    mo->g.gsp = exp (mo->g.gsp);
  }
  #endif

  #if MU0
  if (nmu0)
  {
    mo->a.mu0 /= nmu0;
    mo->g.mu0 /= nmu0;
    mo->g.mu0 = exp (mo->g.mu0);
  }
  #endif

  #if APK
  if (napk)
  {
    mo->a.apk /= napk;
    mo->g.apk /= napk;
    mo->g.apk = exp (mo->g.apk);
  }
  #endif

  #if BTK
  if (nbtk)
  {
    mo->a.btk /= nbtk;
    mo->g.btk /= nbtk;
    mo->g.btk = exp (mo->g.btk);
  }
  #endif

  #if RHK
  if (nrhk)
  {
    mo->a.rhk /= nrhk;
    mo->g.rhk /= nrhk;
    mo->g.rhk = exp (mo->g.rhk);
  }
  #endif

  #if BCK
  if (nbck)
  {
    mo->a.bck /= nbck;
    mo->g.bck /= nbck;
    mo->g.bck = exp (mo->g.bck);
  }
  #endif

  #if BBK
  if (nbbk)
  {
    mo->a.bbk /= nbbk;
    mo->g.bbk /= nbbk;
    mo->g.bbk = exp (mo->g.bbk);
  }
  #endif

  #if BVK
  if (nbvk)
  {
    mo->a.bvk /= nbvk;
    mo->g.bvk /= nbvk;
    mo->g.bvk = exp (mo->g.bvk);
  }
  #endif

  #if BHK
  if (nbhk)
  {
    mo->a.bhk /= nbhk;
    mo->g.bhk /= nbhk;
    mo->g.bhk = exp (mo->g.bhk);
  }
  #endif

  #if ETK
  if (netk)
  {
    mo->a.etk /= netk;
    mo->g.etk /= netk;
    mo->g.etk = exp (mo->g.etk);
  }
  #endif

  #if GCK
  if (ngck)
  {
    mo->a.gck /= ngck;
    mo->g.gck /= ngck;
    mo->g.gck = exp (mo->g.gck);
  }
  #endif

  #if GSK
  if (ngsk)
  {
    mo->a.gsk /= ngsk;
    mo->g.gsk /= ngsk;
    mo->g.gsk = exp (mo->g.gsk);
  }
  #endif

  #if HSK
  if (nhsk)
  {
    mo->a.hsk /= nhsk;
    mo->g.hsk /= nhsk;
    mo->g.hsk = exp (mo->g.hsk);
  }
  #endif
}

static unsigned optimize (int ic, int nc, unsigned nel,
                          unsigned np, unsigned np_l, unsigned nt, unsigned nr,
                          long double X[NX], long double Y[NY], long double Z[NZ],
                          struct Indexes Eijk[np_l][nt][nr],
                          long double R[nr],
                          long double Theta[nt],
                          long double Phi[np_l],
                          struct CartesianPoint spf[nel][NN],
                          MPI_Datatype PARAMETERS,
                          struct Means Mo[nr],
                          struct Parameters a[nel][NX * NY * NZ])
{
  /* Optimizes initial guesses to improve the interpolation quality */
  int kv[nc], kkv[nc];

  np_l = computeGatheringArrays (ic, nc, np, nt, 1, kv, kkv);

  struct Parameters M_l[np_l][nt], M[np][nt];

  clock_t starttime = clock ();

  for (unsigned k = 0; k < nr; k++)
  {
    long double r = R[k];

    for (unsigned j = 0; j < nt; j++)
    {
      long double theta = Theta[j];

      for (unsigned i = 0; i < np_l; i++)
      {
        struct Indexes eijk = Eijk[i][j][k];

        if (eijk.outside_mesh)
        {
          #if VP
          M_l[i][j].vp  = 0.L;
          #endif

          #if VS
          M_l[i][j].vs  = 0.L;
          #endif

          #if RHO
          M_l[i][j].rho = 0.L;
          #endif

          #if VPV
          M_l[i][j].vpv = 0.L;
          #endif

          #if VPH
          M_l[i][j].vph = 0.L;
          #endif

          #if VSV
          M_l[i][j].vsv = 0.L;
          #endif

          #if VSH
          M_l[i][j].vsh = 0.L;
          #endif

          #if ETA
          M_l[i][j].eta = 0.L;
          #endif

          #if QMU
          M_l[i][j].qmu = 0.L;
          #endif

          #if GCP
          M_l[i][j].gcp = 0.L;
          #endif

          #if GSP
          M_l[i][j].gsp = 0.L;
          #endif

          #if MU0
          M_l[i][j].mu0 = 0.L;
          #endif

          #if APK
          M_l[i][j].apk = 0.L;
          #endif

          #if BTK
          M_l[i][j].btk = 0.L;
          #endif

          #if RHK
          M_l[i][j].rhk = 0.L;
          #endif

          #if BCK
          M_l[i][j].bck = 0.L;
          #endif

          #if BBK
          M_l[i][j].bbk = 0.L;
          #endif

          #if BVK
          M_l[i][j].bvk = 0.L;
          #endif

          #if BHK
          M_l[i][j].bhk = 0.L;
          #endif

          #if ETK
          M_l[i][j].etk = 0.L;
          #endif

          #if GCK
          M_l[i][j].gck = 0.L;
          #endif

          #if GSK
          M_l[i][j].gsk = 0.L;
          #endif

          #if HSK
          M_l[i][j].hsk = 0.L;
          #endif

          continue;
        }

        long double phi = Phi[i], x, y, z;

        rThetaPhi2XYZl (r, theta, phi, &x, &y, &z);

        struct Point p;

        p.r     = r;
        p.theta = theta;
        p.phi   = phi;
        p.x     = x;
        p.y     = y;
        p.z     = z;

        long double xi = X[eijk.i], eta = Y[eijk.j], gamma = Z[eijk.k];

        if (newton (&xi, &eta, &gamma, nel, eijk.element, &p, spf))
        {
          fprintf (stderr, "Singular Jacobian. No solution found.\n");

          return 1;
        }

        if (xi    < MIN_XI)    xi    = MIN_XI;
        if (xi    > MAX_XI)    xi    = MAX_XI;
        if (eta   < MIN_ETA)   eta   = MIN_ETA;
        if (eta   > MAX_ETA)   eta   = MAX_ETA;
        if (gamma < MIN_GAMMA) gamma = MIN_GAMMA;
        if (gamma > MAX_GAMMA) gamma = MAX_GAMMA;

        surface (&M_l[i][j], nel, eijk.element, xi, eta, gamma, a);
      }
    }

    MPI_Gatherv (M_l, kv[ic], PARAMETERS, M, kv, kkv, PARAMETERS,
                                                      0,
                                                      MPI_COMM_WORLD);

    computeMeans (&Mo[k], np, nt, M);

    MPI_Barrier (MPI_COMM_WORLD);

    if (ic == 0) progressBar (k, 1, nr, starttime);
  }

  return 0;
}

static void helpMenu (void)
{
  char *help_menu = "\n GLL2MEAN"

                    "\n\n USAGE"
                    "\n    mpiexec -n 36 bin/gll2mean MIN_DEPTH MAX_DEPTH HORIZONTAL_RESOLUTION RADIAL_RESOLUTION"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 36 bin/gll2mean -5 2891 0.5 1.0 data/INPUT/ ."

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    HORIZONTAL_RESOLUTION  - spatial distance (in degrees) between grid points along the horizontal"
                    "\n                             direction"
                    "\n    RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial direction"
                    "\n    INPUT_DIRECTORY        - directory containing the input files"
                    "\n                             grid points of the depth slice"
                    "\n    OUTPUT_DIRECTORY       - directory where the routine will write the output files"

                    "\n\n DESCRIPTION"
                    "\n    Reads the minimum and maximum depths, the horizontal and vertical grid spacings, and the"
                    "\n    input and output directory names from the command line and computes the mean models for the"
                    "\n    desired parameters (defined in 'config.h'). The routine writes the output to files called"
                    "\n    PARAMETER.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  int ic, nc;

  MPI_Init (NULL, NULL);

  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  if (ic == 0 && argc != 7)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command line...\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  long double dmin = atof (argv[1]);
  long double dmax = atof (argv[2]);
  long double dh   = atof (argv[3]);
  long double dv   = atof (argv[4]);

  unsigned np = (unsigned) (360.L / dh + 1.5);
  unsigned nt = (unsigned) (180.L / dh + 1.5);
  unsigned nr = (unsigned) (fabsl (dmax - dmin) / dv + 1.5);

  long double r1 = depth2Rl (dmax);
  long double r2 = depth2Rl (dmin);

  struct Boundaries gb;

  int k[nc], kk[nc];

  unsigned np_l = computeGatheringArrays (ic, nc, np, nt, nr, k, kk);

  long double R[nr], Theta[nt], Phi[np_l];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nCreating output grid...");

  createOutputGrid1D (ic, nc, k, kk, r1, r2, np, np_l, nt, nr, &gb, R, Theta, Phi);

  unsigned nel = 0;

  struct ElNode *lle;
  initializeElNode (&lle);

  struct PmNode *llm;
  initializePmNode (&llm);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nScanning input mesh and model files...\n\n");

  if (checkTopoAndModelIO (scanMeshAndModel (ic, argv[5], NEL, NG, &gb, lle, llm,
                                             &nel))) MPI_Abort (MPI_COMM_WORLD, 1);

  unsigned nels; MPI_Reduce(&nel, &nels, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nAn average of %u espectral elements were selected "
                                "in each slice.\n\n", nels / nc);

  struct Point Ti[nel][NX][NY][NZ];
  struct Parameters Mi[nel][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Converting mesh and model linked lists to arrays...\n");

  toArrayElAndPmNodes (nel, lle, llm, Ti, Mi);

  struct Boundaries b[nel];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing spectral elements approximate boundaries...\n");

  getBoundaries (nel, Ti, b);

  struct Faces f[nel];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing spectral elements faces...\n");

  if (computeFaces (nel, Ti, f))
  {
    fprintf (stderr, "\n\nError: ill-conditioned face or edge!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing faces corrections and uncertainties...\n");

  if (facesCorrections (nel, Ti, f))
  {
    fprintf (stderr, "Error: face uncertainty is too large!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  struct Indexes Eijk[np_l][nt][nr];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nLooking for the closest spectral element "
                                "and internal point...\n\n");

  if (findFirstGuesses (ic, nel, np_l, nt, nr, b, f, Ti, R, Theta, Phi, Eijk))
  {
    fprintf (stderr, "\n\nError: no spectral element contains that point!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  struct CartesianPoint spf[nel][NN];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\n\nComputing shape functions for each "
                                "spectral element...\n");

  if (computeShape (nel, spf, Ti))
  {
    fprintf (stderr, "Error: matrix is singular!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  long double X[NX], Y[NY], Z[NZ];
  long double W[NX][NY][NZ];

  gLLNodesAndWeights3Dl (X, Y, Z, W);

  struct Parameters a[nel][NX * NY * NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Computing spectral elements Lagrange polynomials...\n");

  if (computeInterpolator (nel, a, X, Y, Z, Mi))
  {
    fprintf (stderr, "Error: matrix is singular!\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  struct Means Mo[nr];

  MPI_Datatype PARAMETERS, types[] = {
                                      #if VP
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if VS
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if RHO
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if VPV
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if VPH
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if VSV
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if VSH
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if ETA
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if QMU
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if GCP
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if GSP
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if MU0
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if APK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if BTK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if RHK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if BCK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if BBK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if BVK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if BHK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if ETK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if GCK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if GSK
                                      MPI_LONG_DOUBLE,
                                      #endif

                                      #if HSK
                                      MPI_LONG_DOUBLE
                                      #endif
                                      };

  int blocks[] = {
                  #if VP
                  1,
                  #endif

                  #if VS
                  1,
                  #endif

                  #if RHO
                  1,
                  #endif

                  #if VPV
                  1,
                  #endif

                  #if VPH
                  1,
                  #endif

                  #if VSV
                  1,
                  #endif

                  #if VSH
                  1,
                  #endif

                  #if ETA
                  1,
                  #endif

                  #if QMU
                  1,
                  #endif

                  #if GCP
                  1,
                  #endif

                  #if GSP
                  1,
                  #endif

                  #if MU0
                  1,
                  #endif

                  #if APK
                  1,
                  #endif

                  #if BTK
                  1,
                  #endif

                  #if RHK
                  1,
                  #endif

                  #if BCK
                  1,
                  #endif

                  #if BBK
                  1,
                  #endif

                  #if BVK
                  1,
                  #endif

                  #if BHK
                  1,
                  #endif

                  #if ETK
                  1,
                  #endif

                  #if GCK
                  1,
                  #endif

                  #if GSK
                  1,
                  #endif

                  #if HSK
                  1
                  #endif
                  };

  MPI_Aint displ[] = {
                      #if VP
                      offsetof (struct Parameters, vp),
                      #endif

                      #if VS
                      offsetof (struct Parameters, vs),
                      #endif

                      #if RHO
                      offsetof (struct Parameters, rho),
                      #endif

                      #if VPV
                      offsetof (struct Parameters, vpv),
                      #endif

                      #if VPH
                      offsetof (struct Parameters, vph),
                      #endif

                      #if VSV
                      offsetof (struct Parameters, vsv),
                      #endif

                      #if VSH
                      offsetof (struct Parameters, vsh),
                      #endif

                      #if ETA
                      offsetof (struct Parameters, eta),
                      #endif

                      #if QMU
                      offsetof (struct Parameters, qmu),
                      #endif

                      #if GCP
                      offsetof (struct Parameters, gcp),
                      #endif

                      #if GSP
                      offsetof (struct Parameters, gsp),
                      #endif

                      #if MU0
                      offsetof (struct Parameters, mu0),
                      #endif

                      #if APK
                      offsetof (struct Parameters, apk),
                      #endif

                      #if BTK
                      offsetof (struct Parameters, btk),
                      #endif

                      #if RHK
                      offsetof (struct Parameters, rhk),
                      #endif

                      #if BCK
                      offsetof (struct Parameters, bck),
                      #endif

                      #if BBK
                      offsetof (struct Parameters, bbk),
                      #endif

                      #if BVK
                      offsetof (struct Parameters, bvk),
                      #endif

                      #if BHK
                      offsetof (struct Parameters, bhk),
                      #endif

                      #if ETK
                      offsetof (struct Parameters, etk),
                      #endif

                      #if GCK
                      offsetof (struct Parameters, gck),
                      #endif

                      #if GSK
                      offsetof (struct Parameters, gsk),
                      #endif

                      #if HSK
                      offsetof (struct Parameters, hsk)
                      #endif
                      };

  unsigned count = 0;

  #if VP
  count++;
  #endif

  #if VS
  count++;
  #endif

  #if RHO
  count++;
  #endif

  #if VPV
  count++;
  #endif

  #if VPH
  count++;
  #endif

  #if VSV
  count++;
  #endif

  #if VSH
  count++;
  #endif

  #if ETA
  count++;
  #endif

  #if QMU
  count++;
  #endif

  #if GCP
  count++;
  #endif

  #if GSP
  count++;
  #endif

  #if MU0
  count++;
  #endif

  #if APK
  count++;
  #endif

  #if BTK
  count++;
  #endif

  #if RHK
  count++;
  #endif

  #if BCK
  count++;
  #endif

  #if BBK
  count++;
  #endif

  #if BVK
  count++;
  #endif

  #if BHK
  count++;
  #endif

  #if ETK
  count++;
  #endif

  #if GCK
  count++;
  #endif

  #if GSK
  count++;
  #endif

  #if HSK
  count++;
  #endif

  MPI_Type_create_struct (count, blocks, displ, types, &PARAMETERS);
  MPI_Type_commit (&PARAMETERS);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nOptimizing initial guesses and "
                                "interpolating model parameters...\n\n");

  if (optimize (ic, nc, nel, np, np_l, nt, nr, X, Y, Z,
                Eijk, R, Theta, Phi, spf, PARAMETERS, Mo, a))
  {
    fprintf (stderr, "Singular Jacobian. No solution found.\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Type_free (&PARAMETERS);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\n\nWriting 1D mean model...\n");

    if (checkModelIO (writeModel1D (argv[6], nr,
                                    R, Mo))) MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Done!\n");

  MPI_Finalize ();

  return 0;
}

