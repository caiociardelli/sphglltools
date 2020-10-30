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

 GLL2PF

 USAGE
   mpiexec -n 1 bin/gll2pf MIN_DEPTH MAX_DEPTH LAT LON RESOLUTION INPUT_DIRECTORY
                           OUTPUT_DIRECTORY

 EXAMPLE
   mpiexec -n 1 bin/gll2pf 0 2891 -28 140 1.0 data/INPUT/ .

 COMMAND LINE ARGUMENTS
   MIN_DEPTH              - minimum depth
   MAX_DEPTH              - maximum depth
   LAT LON                - latitude and longitude of the point
   RESOLUTION             - spatial distance (in km) between grid points along the radial
                            direction
   INPUT_DIRECTORY        - directory containing the input files
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads the minimum and maximum depths, the coordinates of the point in which you want the
   profile, the radial grid spacing, and the input and output directory names from the command
   line and creates a one-dimensional profile of the model for the desired model parameters
   (defined in 'config.h'). The routine writes the output to files called PARAMETER_PF.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"
#include "metrics.h"
#include "coordinates.h"
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
                                  unsigned nr,
                                  struct Boundaries b[nel],
                                  struct Faces f[nel],
                                  struct Point Ti[nel][NX][NY][NZ],
                                  struct Point To[nr],
                                  struct Indexes Eijk[nr])
{
  /* Finds initial guesses for the interpolation points */
  clock_t starttime = clock ();

  for (unsigned i = 0; i < nr; i++)
  {
    struct Point p = To[i];
    struct Indexes eijk;

    if (findSpecElm (&p, nel, Ti, b, f, &eijk.element))
    {
      findClosestPoint (&p, nel, eijk.element, Ti, &eijk.i, &eijk.j, &eijk.k);

      eijk.outside_mesh = false;

      Eijk[i] = eijk;
    }

    else Eijk[i].outside_mesh = true;

    if (ic == 0) progressBar (i, 10, nr, starttime);
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

        r++;
      }
}

static unsigned optimize (int ic, unsigned nel,
                          unsigned nr,
                          long double X[NX], long double Y[NY], long double Z[NZ],
                          struct Indexes Eijk[nr],
                          struct Point To[nr],
                          struct CartesianPoint spf[nel][NN],
                          struct Parameters Mo[nr],
                          struct Parameters a[nel][NX * NY * NZ])
{
  /* Optimizes initial guesses to improve the interpolation quality */
  clock_t starttime = clock ();

  for (unsigned i = 0; i < nr; i++)
  {
    struct Indexes eijk = Eijk[i];

    if (eijk.outside_mesh)
    {
      #if VP
      Mo[i].vp  = 0.L;
      #endif

      #if VS
      Mo[i].vs  = 0.L;
      #endif

      #if RHO
      Mo[i].rho = 0.L;
      #endif

      #if VPV
      Mo[i].vpv = 0.L;
      #endif

      #if VPH
      Mo[i].vph = 0.L;
      #endif

      #if VSV
      Mo[i].vsv = 0.L;
      #endif

      #if VSH
      Mo[i].vsh = 0.L;
      #endif

      #if ETA
      Mo[i].eta = 0.L;
      #endif

      #if QMU
      Mo[i].qmu = 0.L;
      #endif
    }

    else
    {
      long double xi = X[eijk.i], eta = Y[eijk.j], gamma = Z[eijk.k];

      if (newton (&xi, &eta, &gamma, nel, eijk.element, &To[i], spf))
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

      surface (&Mo[i], nel, eijk.element, xi, eta, gamma, a);
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (ic == 0) progressBar (i, 10, nr, starttime);
  }

  return 0;
}

static void helpMenu (void)
{
  char *help_menu = "\n GLL2PF"

                    "\n\n USAGE"
                    "\n    mpiexec -n 1 bin/gll2pf MIN_DEPTH MAX_DEPTH LAT LON RESOLUTION INPUT_DIRECTORY"
                    "\n                            OUTPUT_DIRECTORY"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 1 bin/gll2pf 0 2891 -28 140 1.0 data/INPUT/ ."

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    LAT LON                - latitude and longitude of the point"
                    "\n    RESOLUTION             - spatial distance (in km) between grid points along the radial"
                    "\n                             direction"
                    "\n    INPUT_DIRECTORY        - directory containing the input files"
                    "\n    OUTPUT_DIRECTORY       - directory where the routine will write the output files"

                    "\n\n DESCRIPTION"
                    "\n    Reads the minimum and maximum depths, the coordinates of the point in which you want the"
                    "\n    profile, the radial grid spacing, and the input and output directory names from the command"
                    "\n    line and creates a one-dimensional profile of the model for the desired model parameters"
                    "\n    (defined in 'config.h'). The routine writes the output to files called PARAMETER_PF.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  int ic, nc;

  MPI_Init (NULL, NULL);

  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  if (ic == 0 && argc != 8)
  {
    fprintf (stderr, "Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  long double dmin = atof (argv[1]);
  long double dmax = atof (argv[2]);
  long double t    = degree2Radl (90.L - atof (argv[3]));
  long double p    = degree2Radl (atof (argv[4]));
  long double dr   = atof (argv[5]);

  unsigned nr = (unsigned) (fabsl (dmax - dmin) / dr + 1);

  long double r1 = depth2Rl (dmax);
  long double r2 = depth2Rl (dmin);

  struct Boundaries gb;

  if (ic == 0 && (r1 < CMB_R || r1 > MAX_SURFACE_R ||
                  r2 < CMB_R || r2 > MAX_SURFACE_R))
  {
    fprintf (stderr, "\nError: depths should be between %.1Lf and %.1Lf km....\n",
             r2Depthl (MAX_SURFACE_R), r2Depthl (CMB_R));

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  long double R[nr]; struct Point To[nr];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nCreating output grid for vertical cross section...");

  createOutputGridPf (r1, r2, t, p, nr, R, &gb, To);

  unsigned nel = 0;

  struct ElNode *lle;
  initializeElNode (&lle);

  struct PmNode *llm;
  initializePmNode (&llm);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nScanning input mesh and model files...\n\n");

  if (checkTopoAndModelIO (scanMeshAndModel (ic, argv[6], NEL, NG, &gb, lle, llm,
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

  struct Indexes Eijk[nr];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nLooking for the closest spectral element "
                                "and internal point...\n\n");

  if (findFirstGuesses (ic, nel, nr, b, f, Ti, To, Eijk))
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

  struct Parameters Mo[nr];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nOptimizing initial guesses and "
                                "interpolating model parameters...\n\n");

  if (optimize (ic, nel, nr, X, Y, Z, Eijk, To, spf, Mo, a))
  {
    fprintf (stderr, "Singular Jacobian. No solution found.\n");

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0)
  {
    fprintf (stderr, "\n\nWriting profile...\n");

    if (checkModelIO (writeModelPf (argv, nr, R, Mo))) MPI_Abort (MPI_COMM_WORLD, 1);
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Done!\n");

  MPI_Finalize ();

  return 0;
}

