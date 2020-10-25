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

#include <stdbool.h>
#include <math.h>
#include "exmath.h"
#include "structs.h"
#include "coordinates.h"
#include "metrics.h"
#include "boundaries.h"
#include "constants.h"

static inline long double min (long double a, long double b, long double c)
{
  /* Returns the minimum of thee values */
  if (a < b && a < c) return a; else if (b < c) return b; return c;
}

static inline long double crossProduct (long double u1, long double u2,
                                        long double v1, long double v2)
{
  /* Computes the cross-product of two 2D arrays */
  return u1 * v2 - u2 * v1;
}

static char orientation (long double sax, long double say, long double saz,
                         long double ncx, long double ncy, long double ncz,
                         long double srx, long double sry, long double srz)
{
  /* Computes the best orientation for a face of a spectral element */
  long double wl = min (sax, say, saz) + WATER_LEVEL_2;

  sax += wl;
  say += wl;
  saz += wl;

  wl = min (srx, sry, srz) + WATER_LEVEL_2;

  long double nx = ncx * (srx + wl);
  long double ny = ncy * (sry + wl);
  long double nz = ncz * (srz + wl);

  if (nx * say < ny * sax && nx * saz < nz * sax) return 'x';

  else if (ny * saz < nz * say) return 'y';

  return 'z';
}

static inline char position (long double p1, long double p2)
{
  /* Returns the appropriate character according to the
     relative position of two points */
  return p1 > p2 ? '>' : '<';
}

void getBoundaries (unsigned nel, struct Point T[nel][NX][NY][NZ],
                    struct Boundaries b[nel])
{
  /* Gets the boundaries of all spectral elements */
  for (unsigned el = 0; el < nel; el++)
  {
    b[el].xmin =  INFINITY;
    b[el].xmax = -INFINITY;
    b[el].ymin =  INFINITY;
    b[el].ymax = -INFINITY;
    b[el].zmin =  INFINITY;
    b[el].zmax = -INFINITY;
    b[el].rmin =  INFINITY;
    b[el].rmax = -INFINITY;
    b[el].tmin =  INFINITY;
    b[el].tmax = -INFINITY;
    b[el].pmin =  INFINITY;
    b[el].pmax = -INFINITY;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          long double x = T[el][i][j][k].x;
          long double y = T[el][i][j][k].y;
          long double z = T[el][i][j][k].z;

          if (x < b[el].xmin) b[el].xmin = x;
          if (x > b[el].xmax) b[el].xmax = x;
          if (y < b[el].ymin) b[el].ymin = y;
          if (y > b[el].ymax) b[el].ymax = y;
          if (z < b[el].zmin) b[el].zmin = z;
          if (z > b[el].zmax) b[el].zmax = z;

          long double r, theta, phi;

          xYZ2RThetaPhil (x, y, z, &r, &theta, &phi);

          if (r     < b[el].rmin) b[el].rmin = r;
          if (r     > b[el].rmax) b[el].rmax = r;
          if (theta < b[el].tmin) b[el].tmin = theta;
          if (theta > b[el].tmax) b[el].tmax = theta;
          if (phi   < b[el].pmin) b[el].pmin = phi;
          if (phi   > b[el].pmax) b[el].pmax = phi;
        }

    long double tol;

    tol = BOUNDARY_RATIO_1 * (b[el].xmax - b[el].xmin);

    b[el].xmin -= tol;
    b[el].xmax += tol;

    tol = BOUNDARY_RATIO_1 * (b[el].ymax - b[el].ymin);

    b[el].ymin -= tol;
    b[el].ymax += tol;

    tol = BOUNDARY_RATIO_1 * (b[el].zmax - b[el].zmin);

    b[el].zmin -= tol;
    b[el].zmax += tol;

    tol = BOUNDARY_RATIO_1 * (b[el].rmax - b[el].rmin);

    b[el].rmin -= tol;
    b[el].rmax += tol;

    tol = BOUNDARY_RATIO_1 * (b[el].tmax - b[el].tmin);

    b[el].tmin -= tol;
    b[el].tmax += tol;

    tol = BOUNDARY_RATIO_1 * (b[el].pmax - b[el].pmin);

    b[el].pmin -= tol;
    b[el].pmax += tol;
  }
}

static long double fitStraightLine (long double u0, long double u1,
                                    long double v0, long double v1,
                                    struct Curve *c)
{
  /* Fits a straight line given two points */
  long double k1 = u1 - u0;

  if (fabsl (k1) < TOLERANCE) return 1;

  long double k2 = v1 - v0;
  long double k3 = k2 / k1;

  c->a   = v0 - k3 * u0;
  c->au  = k3;
  c->au2 = 0;

  return fabsl (c->a) + fabsl (c->au);
}

static long double fitParabola (long double u0, long double u1, long double u2,
                                long double v0, long double v1, long double v2,
                                struct Curve *c)
{
  /* Fits a parabola given three points */
  long double k1 = (u0 - u1) * (u0 - u2) * (u1 - u2);

  if (fabsl (k1) < TOLERANCE) return INFINITY;

  long double k2 = 1.L / k1;

  c->a   = k2 * (u1 * u2 * (u1 - u2) * v0 +
                 u2 * u0 * (u2 - u0) * v1 +
                 u0 * u1 * (u0 - u1) * v2);
  c->au  = k2 * (squarel (u2) * (v0 - v1) +
                 squarel (u1) * (v2 - v0) +
                 squarel (u0) * (v1 - v2));
  c->au2 = k2 * (u2 * (v1 - v0) +
                 u1 * (v0 - v2) +
                 u0 * (v2 - v1));

  return fabsl (c->a) + fabsl (c->au) + fabsl (c->au2);
}

static inline void computeCorners (long double u0, long double u1,
                                   struct Corners *c, char o)
{
  /* Conputes the corners of an edge of a spectral element */
  c->o = o; if (u0 < u1) {c->c1 = u0; c->c2 = u1;} else {c->c1 = u1; c->c2 = u0;};
}

static unsigned computeEdges (struct Point ep[NNe], struct Edge *e,
                              struct Point cp, char o)
{
  /* Computes the edges of a face of a spectral element */
  struct Curve cx, cy, cz;

  long double ncx, ncy, ncz;

  switch (o)
  {
    case 'x':
      ncx = INFINITY;
      ncy = fitParabola (ep[0].z, ep[1].z, ep[2].z,
                         ep[0].y, ep[1].y, ep[2].y, &cy);
      ncz = fitParabola (ep[0].y, ep[1].y, ep[2].y,
                         ep[0].z, ep[1].z, ep[2].z, &cz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1)
      {
        ncx = INFINITY;
        ncy = fitStraightLine (ep[0].z, ep[2].z, ep[0].y, ep[2].y, &cy);
        ncz = fitStraightLine (ep[0].y, ep[2].y, ep[0].z, ep[2].z, &cz);
      }

      if (min (ncx, ncy, ncz) > MAX_NORM_1) return 1;

      if (ncy < ncz)
      {
        e->o = 'y';
        e->p = position (cp.y, edge (cp.z, &cy));
        e->inner = cy;
        fitStraightLine (ep[0].z, ep[2].z, ep[0].y, ep[2].y, &e->outer);
        computeCorners (ep[0].z, ep[2].z, &e->corners, 'z');
      }

      else
      {
        e->o = 'z';
        e->p = position (cp.z, edge (cp.y, &cz));
        e->inner = cz;
        fitStraightLine (ep[0].y, ep[2].y, ep[0].z, ep[2].z, &e->outer);
        computeCorners (ep[0].y, ep[2].y, &e->corners, 'y');
      }
    break;

    case 'y':
      ncx = fitParabola (ep[0].z, ep[1].z, ep[2].z,
                         ep[0].x, ep[1].x, ep[2].x, &cx);
      ncy = INFINITY;
      ncz = fitParabola (ep[0].x, ep[1].x, ep[2].x,
                         ep[0].z, ep[1].z, ep[2].z, &cz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1)
      {
        ncx = fitStraightLine (ep[0].z, ep[2].z, ep[0].x, ep[2].x, &cx);
        ncy = INFINITY;
        ncz = fitStraightLine (ep[0].x, ep[2].x, ep[0].z, ep[2].z, &cz);
      }

      if (min (ncx, ncy, ncz) > MAX_NORM_1) return 1;

      if (ncx < ncz)
      {
        e->o = 'x';
        e->p = position (cp.x, edge (cp.z, &cx));
        e->inner = cx;
        fitStraightLine (ep[0].z, ep[2].z, ep[0].x, ep[2].x, &e->outer);
        computeCorners (ep[0].z, ep[2].z, &e->corners, 'z');
      }

      else
      {
        e->o = 'z';
        e->p = position (cp.z, edge (cp.x, &cz));
        e->inner = cz;
        fitStraightLine (ep[0].x, ep[2].x, ep[0].z, ep[2].z, &e->outer);
        computeCorners (ep[0].x, ep[2].x, &e->corners, 'x');
      }
    break;

    case 'z':
      ncx = fitParabola (ep[0].y, ep[1].y, ep[2].y,
                         ep[0].x, ep[1].x, ep[2].x, &cx);
      ncy = fitParabola (ep[0].x, ep[1].x, ep[2].x,
                         ep[0].y, ep[1].y, ep[2].y, &cy);
      ncz = INFINITY;

      if (min (ncx, ncy, ncz) > MAX_NORM_1)
      {
        ncx = fitStraightLine (ep[0].y, ep[2].y, ep[0].x, ep[2].x, &cx);
        ncy = fitStraightLine (ep[0].x, ep[2].x, ep[0].y, ep[2].y, &cy);
        ncz = INFINITY;
      }

      if (min (ncx, ncy, ncz) > MAX_NORM_1) return 1;

      if (ncx < ncy)
      {
        e->o = 'x';
        e->p = position (cp.x, edge (cp.y, &cx));
        e->inner = cx;
        fitStraightLine (ep[0].y, ep[2].y, ep[0].x, ep[2].x, &e->outer);
        computeCorners (ep[0].y, ep[2].y, &e->corners, 'y');
      }

      else
      {
        e->o = 'y';
        e->p = position (cp.y, edge (cp.x, &cy));
        e->inner = cy;
        fitStraightLine (ep[0].x, ep[2].x, ep[0].y, ep[2].y, &e->outer);
        computeCorners (ep[0].x, ep[2].x, &e->corners, 'x');
      }
    break;
  }

  e->corners.ti = BOUNDARY_RATIO_3 * distance1D (e->corners.c1, e->corners.c2);

  return 0;
}

static long double fitPlaneSurface (long double C_s[NCF_p][NCF_p], long double c_s[NCF_p],
                                    struct Surface *s)
{
  /* Fits a plane to a face of a spectral element using Least Squares */
  s->a     = 0.L;
  s->au    = 0.L;
  s->au2   = 0.L;
  s->av    = 0.L;
  s->auv   = 0.L;
  s->au2v  = 0.L;
  s->av2   = 0.L;
  s->auv2  = 0.L;
  s->au2v2 = 0.L;

  if (gaussJordanl (NCF_p, C_s)) return INFINITY;

  else

    for (unsigned i = 0; i < NCF_p; i++)
    {
      s->a   += C_s[0][i] * c_s[i];
      s->au  += C_s[1][i] * c_s[i];
      s->av  += C_s[2][i] * c_s[i];
    }

  return fabsl (s->a) + fabsl (s->au) + fabsl (s->av);
}

static long double fitSimpleSurface (long double C_s[SSM][SSM], long double c_s[SSM],
                                     struct Surface *s)
{
  /* Fits a quadratic surface of minimum curvature to a face of a spectral element
     using Lagrange Multipliers */
  s->a     = 0.L;
  s->au    = 0.L;
  s->au2   = 0.L;
  s->av    = 0.L;
  s->auv   = 0.L;
  s->au2v  = 0.L;
  s->av2   = 0.L;
  s->auv2  = 0.L;
  s->au2v2 = 0.L;

  if (gaussJordanl (SSM, C_s)) return INFINITY;

  for (unsigned i = 0; i < SSM; i++)
  {
    s->a   += C_s[0][i] * c_s[i];
    s->au  += C_s[1][i] * c_s[i];
    s->au2 += C_s[2][i] * c_s[i];
    s->av  += C_s[3][i] * c_s[i];
    s->auv += C_s[4][i] * c_s[i];
    s->av2 += C_s[5][i] * c_s[i];
  }

  return fabsl (s->a)  + fabsl (s->au)  + fabsl (s->au2) +
         fabsl (s->av) + fabsl (s->auv) + fabsl (s->av2);
}

static long double fitSurface (long double C[NNf][NNf], long double c[NNf],
                               struct Surface *s)
{
  /* Fits a fourth-degree surface to a face of a spectral element */
  s->a     = 0.L;
  s->au    = 0.L;
  s->au2   = 0.L;
  s->av    = 0.L;
  s->auv   = 0.L;
  s->au2v  = 0.L;
  s->av2   = 0.L;
  s->auv2  = 0.L;
  s->au2v2 = 0.L;

  if (gaussJordanl (NNf, C)) return INFINITY;

  for (unsigned i = 0; i < NNf; i++)
  {
    s->a     += C[0][i] * c[i];
    s->au    += C[1][i] * c[i];
    s->au2   += C[2][i] * c[i];
    s->av    += C[3][i] * c[i];
    s->auv   += C[4][i] * c[i];
    s->au2v  += C[5][i] * c[i];
    s->av2   += C[6][i] * c[i];
    s->auv2  += C[7][i] * c[i];
    s->au2v2 += C[8][i] * c[i];
  }

  return fabsl (s->a)   + fabsl (s->au)   + fabsl (s->au2)  +
         fabsl (s->av)  + fabsl (s->auv)  + fabsl (s->au2v) +
         fabsl (s->av2) + fabsl (s->auv2) + fabsl (s->au2v2);
}

static long double surfaceResiduals (char o, struct Point fn[NPf], struct Surface *s)
{
  /* Computes the residuals of a surface fitted to a face of a spectral element */
  long double sr = 0;

  switch (o)
  {
    case 'x':
      for (unsigned i = 0; i < NPf; i++)

        sr += fabsl (fn[i].x - face (fn[i].y, fn[i].z, s));
    break;

    case 'y':
      for (unsigned i = 0; i < NPf; i++)

        sr += fabsl (fn[i].y - face (fn[i].x, fn[i].z, s));
    break;

    case 'z':
      for (unsigned i = 0; i < NPf; i++)

        sr += fabsl (fn[i].z - face (fn[i].x, fn[i].y, s));
    break;
  }

  return sr > MAX_NORM_2 ? INFINITY : sr;
}

static long double meanCurvature (long double u, long double v, struct Surface *s)
{
  /* Computes the mean curvature for a single point in the surface */
  long double dw_du = s->au +
                      s->auv * v + 2 * s->au2 * u +
                      2 * s->au2v * u * v + s->auv2 * squarel (v) +
                      2 * s->au2v2 * u * squarel (v);

  long double dw_dv = s->av +
                      s->auv * u + 2 * s->av2 * v +
                      s->au2v * squarel (u) + 2 * s->auv2 * u * v +
                      2 * s->au2v2 * squarel (u) * v;

  long double d2w_du2  = 2 * s->au2 +
                         2 * s->au2v * v +
                         2 * s->au2v2 * squarel (v);

  long double d2w_dudv = s->auv +
                         2 * s->au2v * u + 2 * s->auv2 * v +
                         4 * s->au2v2 * u * v;

  long double d2w_dv2  = 2 * s->av2 +
                         2 * s->auv2 * u +
                         2 * s->au2v2 * squarel (u);

  return 0.5 * (d2w_du2 * (1 + squarel (dw_dv))
             - 2 * d2w_dudv * dw_du * dw_dv
             + d2w_dv2 * (1 + squarel (dw_du)))
             / powl (1 + squarel (dw_du) + squarel (dw_dv), 1.5);
}

static long double maxCurvature (char o, struct Point fn[NPf], struct Surface *s)
{
  /* Computes the maximum curvature of a surface
     fitted to a face of a spectral element */
  long double sc = -INFINITY;

  switch (o)
  {
    case 'x':
      for (unsigned i = 0; i < NPf; i++)
      {
        long double H = fabsl (meanCurvature (fn[i].y, fn[i].z, s));

        if (H > sc) sc = H;
      }
    break;

    case 'y':
      for (unsigned i = 0; i < NPf; i++)
      {
        long double H = fabsl (meanCurvature (fn[i].x, fn[i].z, s));

        if (H > sc) sc = H;
      }
    break;

    case 'z':
      for (unsigned i = 0; i < NPf; i++)
      {
        long double H = fabsl (meanCurvature (fn[i].x, fn[i].y, s));

        if (H > sc) sc = H;
      }
    break;
  }

  return sc;
}


unsigned computeFaces (unsigned nel, struct Point T[nel][NX][NY][NZ],
                       struct Faces f[nel])
{
  /* Computes all faces of all spectral element */
  unsigned hnx = NX / 2, hny = NY / 2, hnz = NZ / 2;

  long double Cx1[NNf][NNf], Cy1[NNf][NNf], Cz1[NNf][NNf];
  long double Cx2[NNf][NNf], Cy2[NNf][NNf], Cz2[NNf][NNf];
  long double Cx3[NNf][NNf], Cy3[NNf][NNf], Cz3[NNf][NNf];
  long double Cx4[NNf][NNf], Cy4[NNf][NNf], Cz4[NNf][NNf];
  long double Cx5[NNf][NNf], Cy5[NNf][NNf], Cz5[NNf][NNf];
  long double Cx6[NNf][NNf], Cy6[NNf][NNf], Cz6[NNf][NNf];

  long double Cx1_s[SSM][SSM], Cy1_s[SSM][SSM], Cz1_s[SSM][SSM];
  long double Cx2_s[SSM][SSM], Cy2_s[SSM][SSM], Cz2_s[SSM][SSM];
  long double Cx3_s[SSM][SSM], Cy3_s[SSM][SSM], Cz3_s[SSM][SSM];
  long double Cx4_s[SSM][SSM], Cy4_s[SSM][SSM], Cz4_s[SSM][SSM];
  long double Cx5_s[SSM][SSM], Cy5_s[SSM][SSM], Cz5_s[SSM][SSM];
  long double Cx6_s[SSM][SSM], Cy6_s[SSM][SSM], Cz6_s[SSM][SSM];

  long double Cx1_p[NCF_p][NCF_p], Cy1_p[NCF_p][NCF_p], Cz1_p[NCF_p][NCF_p];
  long double Cx2_p[NCF_p][NCF_p], Cy2_p[NCF_p][NCF_p], Cz2_p[NCF_p][NCF_p];
  long double Cx3_p[NCF_p][NCF_p], Cy3_p[NCF_p][NCF_p], Cz3_p[NCF_p][NCF_p];
  long double Cx4_p[NCF_p][NCF_p], Cy4_p[NCF_p][NCF_p], Cz4_p[NCF_p][NCF_p];
  long double Cx5_p[NCF_p][NCF_p], Cy5_p[NCF_p][NCF_p], Cz5_p[NCF_p][NCF_p];
  long double Cx6_p[NCF_p][NCF_p], Cy6_p[NCF_p][NCF_p], Cz6_p[NCF_p][NCF_p];

  long double cx1[NNf], cy1[NNf], cz1[NNf];
  long double cx2[NNf], cy2[NNf], cz2[NNf];
  long double cx3[NNf], cy3[NNf], cz3[NNf];
  long double cx4[NNf], cy4[NNf], cz4[NNf];
  long double cx5[NNf], cy5[NNf], cz5[NNf];
  long double cx6[NNf], cy6[NNf], cz6[NNf];

  long double cx1_s[SSM], cy1_s[SSM], cz1_s[SSM];
  long double cx2_s[SSM], cy2_s[SSM], cz2_s[SSM];
  long double cx3_s[SSM], cy3_s[SSM], cz3_s[SSM];
  long double cx4_s[SSM], cy4_s[SSM], cz4_s[SSM];
  long double cx5_s[SSM], cy5_s[SSM], cz5_s[SSM];
  long double cx6_s[SSM], cy6_s[SSM], cz6_s[SSM];

  long double cx1_p[NCF_p], cy1_p[NCF_p], cz1_p[NCF_p];
  long double cx2_p[NCF_p], cy2_p[NCF_p], cz2_p[NCF_p];
  long double cx3_p[NCF_p], cy3_p[NCF_p], cz3_p[NCF_p];
  long double cx4_p[NCF_p], cy4_p[NCF_p], cz4_p[NCF_p];
  long double cx5_p[NCF_p], cy5_p[NCF_p], cz5_p[NCF_p];
  long double cx6_p[NCF_p], cy6_p[NCF_p], cz6_p[NCF_p];

  struct Point f1e1[NNe], f1e2[NNe], f1e3[NNe], f1e4[NNe];
  struct Point f2e1[NNe], f2e2[NNe], f2e3[NNe], f2e4[NNe];
  struct Point f3e1[NNe], f3e2[NNe], f3e3[NNe], f3e4[NNe];
  struct Point f4e1[NNe], f4e2[NNe], f4e3[NNe], f4e4[NNe];
  struct Point f5e1[NNe], f5e2[NNe], f5e3[NNe], f5e4[NNe];
  struct Point f6e1[NNe], f6e2[NNe], f6e3[NNe], f6e4[NNe];

  for (unsigned el = 0; el < nel; el++)
  {
    struct Point cpf1 = T[el][0][hny][hnz];
    struct Point cpf2 = T[el][NX - 1][hny][hnz];
    struct Point cpf3 = T[el][hnx][0][hnz];
    struct Point cpf4 = T[el][hnx][NY - 1][hnz];
    struct Point cpf5 = T[el][hnx][hny][0];
    struct Point cpf6 = T[el][hnx][hny][NZ - 1];

    f[el].ti = BOUNDARY_RATIO_2 * distance3D (&cpf1, &cpf2);
    f[el].tj = BOUNDARY_RATIO_2 * distance3D (&cpf3, &cpf4);
    f[el].tk = BOUNDARY_RATIO_2 * distance3D (&cpf5, &cpf6);

    for (unsigned i = 0, j = 0; i < NNe; i++, j += 2)
    {
      f1e1[i] = T[el][0][0][j];
      f1e2[i] = T[el][0][NY - 1][j];
      f1e3[i] = T[el][0][j][0];
      f1e4[i] = T[el][0][j][NZ - 1];

      f2e1[i] = T[el][NX - 1][0][j];
      f2e2[i] = T[el][NX - 1][NY - 1][j];
      f2e3[i] = T[el][NX - 1][j][0];
      f2e4[i] = T[el][NX - 1][j][NZ - 1];

      f3e1[i] = T[el][0][0][j];
      f3e2[i] = T[el][NX - 1][0][j];
      f3e3[i] = T[el][j][0][0];
      f3e4[i] = T[el][j][0][NZ - 1];

      f4e1[i] = T[el][0][NY - 1][j];
      f4e2[i] = T[el][NX - 1][NY - 1][j];
      f4e3[i] = T[el][j][NY - 1][0];
      f4e4[i] = T[el][j][NY - 1][NZ - 1];

      f5e1[i] = T[el][0][j][0];
      f5e2[i] = T[el][NX - 1][j][0];
      f5e3[i] = T[el][j][0][0];
      f5e4[i] = T[el][j][NY - 1][0];

      f6e1[i] = T[el][0][j][NZ - 1];
      f6e2[i] = T[el][NX - 1][j][NZ - 1];
      f6e3[i] = T[el][j][0][NZ - 1];
      f6e4[i] = T[el][j][NY - 1][NZ - 1];
    }

    for (unsigned l = 0, i = 0, j = 0; l < NNf; l++)
    {
      unsigned c = 0;

      struct Point p1 = T[el][0][i][j];
      struct Point p2 = T[el][NX - 1][i][j];
      struct Point p3 = T[el][i][0][j];
      struct Point p4 = T[el][i][NY - 1][j];
      struct Point p5 = T[el][i][j][0];
      struct Point p6 = T[el][i][j][NZ - 1];

      for (unsigned e2 = 0; e2 < NNe; e2++)

        for (unsigned e1 = 0; e1 < NNe; e1++)
        {
          Cx1[l][c] = cNN (p1.y, p1.z, e1, e2);
          Cy1[l][c] = cNN (p1.x, p1.z, e1, e2);
          Cz1[l][c] = cNN (p1.x, p1.y, e1, e2);

          Cx2[l][c] = cNN (p2.y, p2.z, e1, e2);
          Cy2[l][c] = cNN (p2.x, p2.z, e1, e2);
          Cz2[l][c] = cNN (p2.x, p2.y, e1, e2);

          Cx3[l][c] = cNN (p3.y, p3.z, e1, e2);
          Cy3[l][c] = cNN (p3.x, p3.z, e1, e2);
          Cz3[l][c] = cNN (p3.x, p3.y, e1, e2);

          Cx4[l][c] = cNN (p4.y, p4.z, e1, e2);
          Cy4[l][c] = cNN (p4.x, p4.z, e1, e2);
          Cz4[l][c] = cNN (p4.x, p4.y, e1, e2);

          Cx5[l][c] = cNN (p5.y, p5.z, e1, e2);
          Cy5[l][c] = cNN (p5.x, p5.z, e1, e2);
          Cz5[l][c] = cNN (p5.x, p5.y, e1, e2);

          Cx6[l][c] = cNN (p6.y, p6.z, e1, e2);
          Cy6[l][c] = cNN (p6.x, p6.z, e1, e2);
          Cz6[l][c] = cNN (p6.x, p6.y, e1, e2);

          c++;
        }

      cx1[l] = p1.x; cy1[l] = p1.y; cz1[l] = p1.z;
      cx2[l] = p2.x; cy2[l] = p2.y; cz2[l] = p2.z;
      cx3[l] = p3.x; cy3[l] = p3.y; cz3[l] = p3.z;
      cx4[l] = p4.x; cy4[l] = p4.y; cz4[l] = p4.z;
      cx5[l] = p5.x; cy5[l] = p5.y; cz5[l] = p5.z;
      cx6[l] = p6.x; cy6[l] = p6.y; cz6[l] = p6.z;

      i += 2;

      if (i >= NPe) /* Assumes that NX = NY = NZ */
      {
        i = 0; j += 2;
      }
    }

    for (unsigned l = 0; l < NCF_p; l++)

      for (unsigned c = 0; c < NCF_p; c++)
      {
        Cx1_p[l][c] = 0.L; Cy1_p[l][c] = 0.L; Cz1_p[l][c] = 0.L;
        Cx2_p[l][c] = 0.L; Cy2_p[l][c] = 0.L; Cz2_p[l][c] = 0.L;
        Cx3_p[l][c] = 0.L; Cy3_p[l][c] = 0.L; Cz3_p[l][c] = 0.L;
        Cx4_p[l][c] = 0.L; Cy4_p[l][c] = 0.L; Cz4_p[l][c] = 0.L;
        Cx5_p[l][c] = 0.L; Cy5_p[l][c] = 0.L; Cz5_p[l][c] = 0.L;
        Cx6_p[l][c] = 0.L; Cy6_p[l][c] = 0.L; Cz6_p[l][c] = 0.L;

        cx1_p[l] = 0.L; cy1_p[l] = 0.L; cz1_p[l] = 0.L;
        cx2_p[l] = 0.L; cy2_p[l] = 0.L; cz2_p[l] = 0.L;
        cx3_p[l] = 0.L; cy3_p[l] = 0.L; cz3_p[l] = 0.L;
        cx4_p[l] = 0.L; cy4_p[l] = 0.L; cz4_p[l] = 0.L;
        cx5_p[l] = 0.L; cy5_p[l] = 0.L; cz5_p[l] = 0.L;
        cx6_p[l] = 0.L; cy6_p[l] = 0.L; cz6_p[l] = 0.L;
      }

    for (unsigned l = 0; l < SSM; l++)

      for (unsigned c = 0; c < SSM; c++)
      {
        Cx1_s[l][c] = 0.L; Cy1_s[l][c] = 0.L; Cz1_s[l][c] = 0.L;
        Cx2_s[l][c] = 0.L; Cy2_s[l][c] = 0.L; Cz2_s[l][c] = 0.L;
        Cx3_s[l][c] = 0.L; Cy3_s[l][c] = 0.L; Cz3_s[l][c] = 0.L;
        Cx4_s[l][c] = 0.L; Cy4_s[l][c] = 0.L; Cz4_s[l][c] = 0.L;
        Cx5_s[l][c] = 0.L; Cy5_s[l][c] = 0.L; Cz5_s[l][c] = 0.L;
        Cx6_s[l][c] = 0.L; Cy6_s[l][c] = 0.L; Cz6_s[l][c] = 0.L;

        cx1_s[l] = 0.L; cy1_s[l] = 0.L; cz1_s[l] = 0.L;
        cx2_s[l] = 0.L; cy2_s[l] = 0.L; cz2_s[l] = 0.L;
        cx3_s[l] = 0.L; cy3_s[l] = 0.L; cz3_s[l] = 0.L;
        cx4_s[l] = 0.L; cy4_s[l] = 0.L; cz4_s[l] = 0.L;
        cx5_s[l] = 0.L; cy5_s[l] = 0.L; cz5_s[l] = 0.L;
        cx6_s[l] = 0.L; cy6_s[l] = 0.L; cz6_s[l] = 0.L;
      }

    Cx1_s[0][0] = 1.L; Cy1_s[0][0] = 1.L; Cz1_s[0][0] = 1.L;
    Cx1_s[1][1] = 1.L; Cy1_s[1][1] = 1.L; Cz1_s[1][1] = 1.L;
    Cx1_s[2][2] = 1.L; Cy1_s[2][2] = 1.L; Cz1_s[2][2] = 1.L;

    Cx2_s[0][0] = 1.L; Cy2_s[0][0] = 1.L; Cz2_s[0][0] = 1.L;
    Cx2_s[1][1] = 1.L; Cy2_s[1][1] = 1.L; Cz2_s[1][1] = 1.L;
    Cx2_s[2][2] = 1.L; Cy2_s[2][2] = 1.L; Cz2_s[2][2] = 1.L;

    Cx3_s[0][0] = 1.L; Cy3_s[0][0] = 1.L; Cz3_s[0][0] = 1.L;
    Cx3_s[1][1] = 1.L; Cy3_s[1][1] = 1.L; Cz3_s[1][1] = 1.L;
    Cx3_s[2][2] = 1.L; Cy3_s[2][2] = 1.L; Cz3_s[2][2] = 1.L;

    Cx4_s[0][0] = 1.L; Cy4_s[0][0] = 1.L; Cz4_s[0][0] = 1.L;
    Cx4_s[1][1] = 1.L; Cy4_s[1][1] = 1.L; Cz4_s[1][1] = 1.L;
    Cx4_s[2][2] = 1.L; Cy4_s[2][2] = 1.L; Cz4_s[2][2] = 1.L;

    Cx5_s[0][0] = 1.L; Cy5_s[0][0] = 1.L; Cz5_s[0][0] = 1.L;
    Cx5_s[1][1] = 1.L; Cy5_s[1][1] = 1.L; Cz5_s[1][1] = 1.L;
    Cx5_s[2][2] = 1.L; Cy5_s[2][2] = 1.L; Cz5_s[2][2] = 1.L;

    Cx6_s[0][0] = 1.L; Cy6_s[0][0] = 1.L; Cz6_s[0][0] = 1.L;
    Cx6_s[1][1] = 1.L; Cy6_s[1][1] = 1.L; Cz6_s[1][1] = 1.L;
    Cx6_s[2][2] = 1.L; Cy6_s[2][2] = 1.L; Cz6_s[2][2] = 1.L;

    unsigned l = NCF_s;

    for (unsigned j = 0; j < NPe; j += 2)

      for (unsigned i = 0; i < NPe; i += 2)
      {
        struct Point p1 = T[el][0][i][j];
        struct Point p2 = T[el][NX - 1][i][j];
        struct Point p3 = T[el][i][0][j];
        struct Point p4 = T[el][i][NY - 1][j];
        struct Point p5 = T[el][i][j][0];
        struct Point p6 = T[el][i][j][NZ - 1];

        Cx1_p[0][0] += 1.L; Cx1_p[0][1] += p1.y;        Cx1_p[0][2] += p1.z;
        Cx1_p[1][0] += p1.y; Cx1_p[1][1] += p1.y * p1.y; Cx1_p[1][2] += p1.y * p1.z;
        Cx1_p[2][0] += p1.z; Cx1_p[2][1] += p1.y * p1.z; Cx1_p[2][2] += p1.z * p1.z;
        Cy1_p[0][0] += 1.L; Cy1_p[0][1] += p1.x;        Cy1_p[0][2] += p1.z;
        Cy1_p[1][0] += p1.x; Cy1_p[1][1] += p1.x * p1.x; Cy1_p[1][2] += p1.x * p1.z;
        Cy1_p[2][0] += p1.z; Cy1_p[2][1] += p1.x * p1.z; Cy1_p[2][2] += p1.z * p1.z;
        Cz1_p[0][0] += 1.L; Cz1_p[0][1] += p1.x;        Cz1_p[0][2] += p1.y;
        Cz1_p[1][0] += p1.x; Cz1_p[1][1] += p1.x * p1.x; Cz1_p[1][2] += p1.x * p1.y;
        Cz1_p[2][0] += p1.y; Cz1_p[2][1] += p1.x * p1.y; Cz1_p[2][2] += p1.y * p1.y;

        Cx2_p[0][0] += 1.L; Cx2_p[0][1] += p2.y;        Cx2_p[0][2] += p2.z;
        Cx2_p[1][0] += p2.y; Cx2_p[1][1] += p2.y * p2.y; Cx2_p[1][2] += p2.y * p2.z;
        Cx2_p[2][0] += p2.z; Cx2_p[2][1] += p2.y * p2.z; Cx2_p[2][2] += p2.z * p2.z;
        Cy2_p[0][0] += 1.L; Cy2_p[0][1] += p2.x;        Cy2_p[0][2] += p2.z;
        Cy2_p[1][0] += p2.x; Cy2_p[1][1] += p2.x * p2.x; Cy2_p[1][2] += p2.x * p2.z;
        Cy2_p[2][0] += p2.z; Cy2_p[2][1] += p2.x * p2.z; Cy2_p[2][2] += p2.z * p2.z;
        Cz2_p[0][0] += 1.L; Cz2_p[0][1] += p2.x;        Cz2_p[0][2] += p2.y;
        Cz2_p[1][0] += p2.x; Cz2_p[1][1] += p2.x * p2.x; Cz2_p[1][2] += p2.x * p2.y;
        Cz2_p[2][0] += p2.y; Cz2_p[2][1] += p2.x * p2.y; Cz2_p[2][2] += p2.y * p2.y;

        Cx3_p[0][0] += 1.L; Cx3_p[0][1] += p3.y;        Cx3_p[0][2] += p3.z;
        Cx3_p[1][0] += p3.y; Cx3_p[1][1] += p3.y * p3.y; Cx3_p[1][2] += p3.y * p3.z;
        Cx3_p[2][0] += p3.z; Cx3_p[2][1] += p3.y * p3.z; Cx3_p[2][2] += p3.z * p3.z;
        Cy3_p[0][0] += 1.L; Cy3_p[0][1] += p3.x;        Cy3_p[0][2] += p3.z;
        Cy3_p[1][0] += p3.x; Cy3_p[1][1] += p3.x * p3.x; Cy3_p[1][2] += p3.x * p3.z;
        Cy3_p[2][0] += p3.z; Cy3_p[2][1] += p3.x * p3.z; Cy3_p[2][2] += p3.z * p3.z;
        Cz3_p[0][0] += 1.L; Cz3_p[0][1] += p3.x;        Cz3_p[0][2] += p3.y;
        Cz3_p[1][0] += p3.x; Cz3_p[1][1] += p3.x * p3.x; Cz3_p[1][2] += p3.x * p3.y;
        Cz3_p[2][0] += p3.y; Cz3_p[2][1] += p3.x * p3.y; Cz3_p[2][2] += p3.y * p3.y;

        Cx4_p[0][0] += 1.L; Cx4_p[0][1] += p4.y;        Cx4_p[0][2] += p4.z;
        Cx4_p[1][0] += p4.y; Cx4_p[1][1] += p4.y * p4.y; Cx4_p[1][2] += p4.y * p4.z;
        Cx4_p[2][0] += p4.z; Cx4_p[2][1] += p4.y * p4.z; Cx4_p[2][2] += p4.z * p4.z;
        Cy4_p[0][0] += 1.L; Cy4_p[0][1] += p4.x;        Cy4_p[0][2] += p4.z;
        Cy4_p[1][0] += p4.x; Cy4_p[1][1] += p4.x * p4.x; Cy4_p[1][2] += p4.x * p4.z;
        Cy4_p[2][0] += p4.z; Cy4_p[2][1] += p4.x * p4.z; Cy4_p[2][2] += p4.z * p4.z;
        Cz4_p[0][0] += 1.L; Cz4_p[0][1] += p4.x;        Cz4_p[0][2] += p4.y;
        Cz4_p[1][0] += p4.x; Cz4_p[1][1] += p4.x * p4.x; Cz4_p[1][2] += p4.x * p4.y;
        Cz4_p[2][0] += p4.y; Cz4_p[2][1] += p4.x * p4.y; Cz4_p[2][2] += p4.y * p4.y;

        Cx5_p[0][0] += 1.L; Cx5_p[0][1] += p5.y;        Cx5_p[0][2] += p5.z;
        Cx5_p[1][0] += p5.y; Cx5_p[1][1] += p5.y * p5.y; Cx5_p[1][2] += p5.y * p5.z;
        Cx5_p[2][0] += p5.z; Cx5_p[2][1] += p5.y * p5.z; Cx5_p[2][2] += p5.z * p5.z;
        Cy5_p[0][0] += 1.L; Cy5_p[0][1] += p5.x;        Cy5_p[0][2] += p5.z;
        Cy5_p[1][0] += p5.x; Cy5_p[1][1] += p5.x * p5.x; Cy5_p[1][2] += p5.x * p5.z;
        Cy5_p[2][0] += p5.z; Cy5_p[2][1] += p5.x * p5.z; Cy5_p[2][2] += p5.z * p5.z;
        Cz5_p[0][0] += 1.L; Cz5_p[0][1] += p5.x;        Cz5_p[0][2] += p5.y;
        Cz5_p[1][0] += p5.x; Cz5_p[1][1] += p5.x * p5.x; Cz5_p[1][2] += p5.x * p5.y;
        Cz5_p[2][0] += p5.y; Cz5_p[2][1] += p5.x * p5.y; Cz5_p[2][2] += p5.y * p5.y;

        Cx6_p[0][0] += 1.L; Cx6_p[0][1] += p6.y;        Cx6_p[0][2] += p6.z;
        Cx6_p[1][0] += p6.y; Cx6_p[1][1] += p6.y * p6.y; Cx6_p[1][2] += p6.y * p6.z;
        Cx6_p[2][0] += p6.z; Cx6_p[2][1] += p6.y * p6.z; Cx6_p[2][2] += p6.z * p6.z;
        Cy6_p[0][0] += 1.L; Cy6_p[0][1] += p6.x;        Cy6_p[0][2] += p6.z;
        Cy6_p[1][0] += p6.x; Cy6_p[1][1] += p6.x * p6.x; Cy6_p[1][2] += p6.x * p6.z;
        Cy6_p[2][0] += p6.z; Cy6_p[2][1] += p6.x * p6.z; Cy6_p[2][2] += p6.z * p6.z;
        Cz6_p[0][0] += 1.L; Cz6_p[0][1] += p6.x;        Cz6_p[0][2] += p6.y;
        Cz6_p[1][0] += p6.x; Cz6_p[1][1] += p6.x * p6.x; Cz6_p[1][2] += p6.x * p6.y;
        Cz6_p[2][0] += p6.y; Cz6_p[2][1] += p6.x * p6.y; Cz6_p[2][2] += p6.y * p6.y;

        cx1_p[0] += p1.x; cx1_p[1] += p1.x * p1.y; cx1_p[2] += p1.x * p1.z;
        cy1_p[0] += p1.y; cy1_p[1] += p1.x * p1.y; cy1_p[2] += p1.y * p1.z;
        cz1_p[0] += p1.z; cz1_p[1] += p1.x * p1.z; cz1_p[2] += p1.y * p1.z;

        cx2_p[0] += p2.x; cx2_p[1] += p2.x * p2.y; cx2_p[2] += p2.x * p2.z;
        cy2_p[0] += p2.y; cy2_p[1] += p2.x * p2.y; cy2_p[2] += p2.y * p2.z;
        cz2_p[0] += p2.z; cz2_p[1] += p2.x * p2.z; cz2_p[2] += p2.y * p2.z;

        cx3_p[0] += p3.x; cx3_p[1] += p3.x * p3.y; cx3_p[2] += p3.x * p3.z;
        cy3_p[0] += p3.y; cy3_p[1] += p3.x * p3.y; cy3_p[2] += p3.y * p3.z;
        cz3_p[0] += p3.z; cz3_p[1] += p3.x * p3.z; cz3_p[2] += p3.y * p3.z;

        cx4_p[0] += p4.x; cx4_p[1] += p4.x * p4.y; cx4_p[2] += p4.x * p4.z;
        cy4_p[0] += p4.y; cy4_p[1] += p4.x * p4.y; cy4_p[2] += p4.y * p4.z;
        cz4_p[0] += p4.z; cz4_p[1] += p4.x * p4.z; cz4_p[2] += p4.y * p4.z;

        cx5_p[0] += p5.x; cx5_p[1] += p5.x * p5.y; cx5_p[2] += p5.x * p5.z;
        cy5_p[0] += p5.y; cy5_p[1] += p5.x * p5.y; cy5_p[2] += p5.y * p5.z;
        cz5_p[0] += p5.z; cz5_p[1] += p5.x * p5.z; cz5_p[2] += p5.y * p5.z;

        cx6_p[0] += p6.x; cx6_p[1] += p6.x * p6.y; cx6_p[2] += p6.x * p6.z;
        cy6_p[0] += p6.y; cy6_p[1] += p6.x * p6.y; cy6_p[2] += p6.y * p6.z;
        cz6_p[0] += p6.z; cz6_p[1] += p6.x * p6.z; cz6_p[2] += p6.y * p6.z;

        if ((i == 2 && j != 2) || (j == 2 && i != 2)) continue;

        Cx1_s[l][0] = 1.L; Cx1_s[l][1] = p1.y;        Cx1_s[l][2] = p1.y * p1.y;
        Cx1_s[l][3] = p1.z; Cx1_s[l][4] = p1.y * p1.z; Cx1_s[l][5] = p1.z * p1.z;
        Cy1_s[l][0] = 1.L; Cy1_s[l][1] = p1.x;        Cy1_s[l][2] = p1.x * p1.x;
        Cy1_s[l][3] = p1.z; Cy1_s[l][4] = p1.x * p1.z; Cy1_s[l][5] = p1.z * p1.z;
        Cz1_s[l][0] = 1.L; Cz1_s[l][1] = p1.x;        Cz1_s[l][2] = p1.x * p1.x;
        Cz1_s[l][3] = p1.y; Cz1_s[l][4] = p1.x * p1.y; Cz1_s[l][5] = p1.y * p1.y;

        Cx2_s[l][0] = 1.L; Cx2_s[l][1] = p2.y;        Cx2_s[l][2] = p2.y * p2.y;
        Cx2_s[l][3] = p2.z; Cx2_s[l][4] = p2.y * p2.z; Cx2_s[l][5] = p2.z * p2.z;
        Cy2_s[l][0] = 1.L; Cy2_s[l][1] = p2.x;        Cy2_s[l][2] = p2.x * p2.x;
        Cy2_s[l][3] = p2.z; Cy2_s[l][4] = p2.x * p2.z; Cy2_s[l][5] = p2.z * p2.z;
        Cz2_s[l][0] = 1.L; Cz2_s[l][1] = p2.x;        Cz2_s[l][2] = p2.x * p2.x;
        Cz2_s[l][3] = p2.y; Cz2_s[l][4] = p2.x * p2.y; Cz2_s[l][5] = p2.y * p2.y;

        Cx3_s[l][0] = 1.L; Cx3_s[l][1] = p3.y;        Cx3_s[l][2] = p3.y * p3.y;
        Cx3_s[l][3] = p3.z; Cx3_s[l][4] = p3.y * p3.z; Cx3_s[l][5] = p3.z * p3.z;
        Cy3_s[l][0] = 1.L; Cy3_s[l][1] = p3.x;        Cy3_s[l][2] = p3.x * p3.x;
        Cy3_s[l][3] = p3.z; Cy3_s[l][4] = p3.x * p3.z; Cy3_s[l][5] = p3.z * p3.z;
        Cz3_s[l][0] = 1.L; Cz3_s[l][1] = p3.x;        Cz3_s[l][2] = p3.x * p3.x;
        Cz3_s[l][3] = p3.y; Cz3_s[l][4] = p3.x * p3.y; Cz3_s[l][5] = p3.y * p3.y;

        Cx4_s[l][0] = 1.L; Cx4_s[l][1] = p4.y;        Cx4_s[l][2] = p4.y * p4.y;
        Cx4_s[l][3] = p4.z; Cx4_s[l][4] = p4.y * p4.z; Cx4_s[l][5] = p4.z * p4.z;
        Cy4_s[l][0] = 1.L; Cy4_s[l][1] = p4.x;        Cy4_s[l][2] = p4.x * p4.x;
        Cy4_s[l][3] = p4.z; Cy4_s[l][4] = p4.x * p4.z; Cy4_s[l][5] = p4.z * p4.z;
        Cz4_s[l][0] = 1.L; Cz4_s[l][1] = p4.x;        Cz4_s[l][2] = p4.x * p4.x;
        Cz4_s[l][3] = p4.y; Cz4_s[l][4] = p4.x * p4.y; Cz4_s[l][5] = p4.y * p4.y;

        Cx5_s[l][0] = 1.L; Cx5_s[l][1] = p5.y;        Cx5_s[l][2] = p5.y * p5.y;
        Cx5_s[l][3] = p5.z; Cx5_s[l][4] = p5.y * p5.z; Cx5_s[l][5] = p5.z * p5.z;
        Cy5_s[l][0] = 1.L; Cy5_s[l][1] = p5.x;        Cy5_s[l][2] = p5.x * p5.x;
        Cy5_s[l][3] = p5.z; Cy5_s[l][4] = p5.x * p5.z; Cy5_s[l][5] = p5.z * p5.z;
        Cz5_s[l][0] = 1.L; Cz5_s[l][1] = p5.x;        Cz5_s[l][2] = p5.x * p5.x;
        Cz5_s[l][3] = p5.y; Cz5_s[l][4] = p5.x * p5.y; Cz5_s[l][5] = p5.y * p5.y;

        Cx6_s[l][0] = 1.L; Cx6_s[l][1] = p6.y;        Cx6_s[l][2] = p6.y * p6.y;
        Cx6_s[l][3] = p6.z; Cx6_s[l][4] = p6.y * p6.z; Cx6_s[l][5] = p6.z * p6.z;
        Cy6_s[l][0] = 1.L; Cy6_s[l][1] = p6.x;        Cy6_s[l][2] = p6.x * p6.x;
        Cy6_s[l][3] = p6.z; Cy6_s[l][4] = p6.x * p6.z; Cy6_s[l][5] = p6.z * p6.z;
        Cz6_s[l][0] = 1.L; Cz6_s[l][1] = p6.x;        Cz6_s[l][2] = p6.x * p6.x;
        Cz6_s[l][3] = p6.y; Cz6_s[l][4] = p6.x * p6.y; Cz6_s[l][5] = p6.y * p6.y;

        for (unsigned ii = 0; ii < NCF_s; ii++)
        {
          Cx1_s[ii][l] = Cx1_s[l][ii];
          Cy1_s[ii][l] = Cy1_s[l][ii];
          Cz1_s[ii][l] = Cz1_s[l][ii];

          Cx2_s[ii][l] = Cx2_s[l][ii];
          Cy2_s[ii][l] = Cy2_s[l][ii];
          Cz2_s[ii][l] = Cz2_s[l][ii];

          Cx3_s[ii][l] = Cx3_s[l][ii];
          Cy3_s[ii][l] = Cy3_s[l][ii];
          Cz3_s[ii][l] = Cz3_s[l][ii];

          Cx4_s[ii][l] = Cx4_s[l][ii];
          Cy4_s[ii][l] = Cy4_s[l][ii];
          Cz4_s[ii][l] = Cz4_s[l][ii];

          Cx5_s[ii][l] = Cx5_s[l][ii];
          Cy5_s[ii][l] = Cy5_s[l][ii];
          Cz5_s[ii][l] = Cz5_s[l][ii];

          Cx6_s[ii][l] = Cx6_s[l][ii];
          Cy6_s[ii][l] = Cy6_s[l][ii];
          Cz6_s[ii][l] = Cz6_s[l][ii];
        }

        cx1_s[l] = p1.x; cy1_s[l] = p1.y; cz1_s[l] = p1.z;
        cx2_s[l] = p2.x; cy2_s[l] = p2.y; cz2_s[l] = p2.z;
        cx3_s[l] = p3.x; cy3_s[l] = p3.y; cz3_s[l] = p3.z;
        cx4_s[l] = p4.x; cy4_s[l] = p4.y; cz4_s[l] = p4.z;
        cx5_s[l] = p5.x; cy5_s[l] = p5.y; cz5_s[l] = p5.z;
        cx6_s[l] = p6.x; cy6_s[l] = p6.y; cz6_s[l] = p6.z;

        l++;
      }

    struct Point f1n[NPf], f2n[NPf], f3n[NPf], f4n[NPf], f5n[NPf], f6n[NPf];

    unsigned k = 0;

    for (unsigned i = 0; i < NPe; i++)

      for (unsigned j = 0; j < NPe; j++)
      {
        f1n[k] = T[el][0][i][j];
        f2n[k] = T[el][NX - 1][i][j];
        f3n[k] = T[el][i][0][j];
        f4n[k] = T[el][i][NY - 1][j];
        f5n[k] = T[el][i][j][0];
        f6n[k] = T[el][i][j][NZ - 1];

        k++;
      }

    long double sax, say, saz;
    long double ncx, ncy, ncz;
    long double srx, sry, srz;
    long double scx, scy, scz;

    long double p0x, p1x, p2x, p3x;
    long double p0y, p1y, p2y, p3y;
    long double p0z, p1z, p2z, p3z;

    struct Point cp = T[el][hnx][hny][hnz];
    struct Surface sx, sy, sz;
    struct Point cpe1, cpe2, cpe3, cpe4;

    p0x = T[el][0][0][0].x;
    p0y = T[el][0][0][0].y;
    p0z = T[el][0][0][0].z;
    p1x = T[el][0][NY - 1][0].x;
    p1y = T[el][0][NY - 1][0].y;
    p1z = T[el][0][NY - 1][0].z;
    p2x = T[el][0][NY - 1][NZ - 1].x;
    p2y = T[el][0][NY - 1][NZ - 1].y;
    p2z = T[el][0][NY - 1][NZ - 1].z;
    p3x = T[el][0][0][NZ - 1].x;
    p3y = T[el][0][0][NZ - 1].y;
    p3z = T[el][0][0][NZ - 1].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx1, cx1, &sx);
    ncy = fitSurface (Cy1, cy1, &sy);
    ncz = fitSurface (Cz1, cz1, &sz);

    srx = surfaceResiduals ('x', f1n, &sx);
    sry = surfaceResiduals ('y', f1n, &sy);
    srz = surfaceResiduals ('z', f1n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f1n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f1n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f1n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx1_s, cx1_s, &sx);
      ncy = fitSimpleSurface (Cy1_s, cy1_s, &sy);
      ncz = fitSimpleSurface (Cz1_s, cz1_s, &sz);

      srx = surfaceResiduals ('x', f1n, &sx);
      sry = surfaceResiduals ('y', f1n, &sy);
      srz = surfaceResiduals ('z', f1n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f1.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f1e1, &f[el].f1.edges.e1, cpf1, f[el].f1.o)) return 1;
    if (computeEdges (f1e2, &f[el].f1.edges.e2, cpf1, f[el].f1.o)) return 1;
    if (computeEdges (f1e3, &f[el].f1.edges.e3, cpf1, f[el].f1.o)) return 1;
    if (computeEdges (f1e4, &f[el].f1.edges.e4, cpf1, f[el].f1.o)) return 1;

    cpe1 = f1e1[1];
    cpe2 = f1e2[1];
    cpe3 = f1e3[1];
    cpe4 = f1e4[1];

    switch (f[el].f1.o)
    {
      case 'x':
        f[el].f1.inner = sx;
        f[el].f1.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx1_p, cx1_p, &f[el].f1.outer);
        f[el].f1.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f1.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f1.inner = sy;
        f[el].f1.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy1_p, cy1_p, &f[el].f1.outer);
        f[el].f1.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f1.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f1.inner = sz;
        f[el].f1.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz1_p, cz1_p, &f[el].f1.outer);
        f[el].f1.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f1.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }

    p0x = T[el][NX - 1][0][0].x;
    p0y = T[el][NX - 1][0][0].y;
    p0z = T[el][NX - 1][0][0].z;
    p1x = T[el][NX - 1][NY - 1][0].x;
    p1y = T[el][NX - 1][NY - 1][0].y;
    p1z = T[el][NX - 1][NY - 1][0].z;
    p2x = T[el][NX - 1][NY - 1][NZ - 1].x;
    p2y = T[el][NX - 1][NY - 1][NZ - 1].y;
    p2z = T[el][NX - 1][NY - 1][NZ - 1].z;
    p3x = T[el][NX - 1][0][NZ - 1].x;
    p3y = T[el][NX - 1][0][NZ - 1].y;
    p3z = T[el][NX - 1][0][NZ - 1].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx2, cx2, &sx);
    ncy = fitSurface (Cy2, cy2, &sy);
    ncz = fitSurface (Cz2, cz2, &sz);

    srx = surfaceResiduals ('x', f2n, &sx);
    sry = surfaceResiduals ('y', f2n, &sy);
    srz = surfaceResiduals ('z', f2n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f2n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f2n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f2n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx2_s, cx2_s, &sx);
      ncy = fitSimpleSurface (Cy2_s, cy2_s, &sy);
      ncz = fitSimpleSurface (Cz2_s, cz2_s, &sz);

      srx = surfaceResiduals ('x', f2n, &sx);
      sry = surfaceResiduals ('y', f2n, &sy);
      srz = surfaceResiduals ('z', f2n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f2.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f2e1, &f[el].f2.edges.e1, cpf2, f[el].f2.o)) return 1;
    if (computeEdges (f2e2, &f[el].f2.edges.e2, cpf2, f[el].f2.o)) return 1;
    if (computeEdges (f2e3, &f[el].f2.edges.e3, cpf2, f[el].f2.o)) return 1;
    if (computeEdges (f2e4, &f[el].f2.edges.e4, cpf2, f[el].f2.o)) return 1;

    cpe1 = f2e1[1];
    cpe2 = f2e2[1];
    cpe3 = f2e3[1];
    cpe4 = f2e4[1];

    switch (f[el].f2.o)
    {
      case 'x':
        f[el].f2.inner = sx;
        f[el].f2.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx2_p, cx2_p, &f[el].f2.outer);
        f[el].f2.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f2.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f2.inner = sy;
        f[el].f2.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy2_p, cy2_p, &f[el].f2.outer);
        f[el].f2.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f2.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f2.inner = sz;
        f[el].f2.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz2_p, cz2_p, &f[el].f2.outer);
        f[el].f2.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f2.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }

    p0x = T[el][0][0][0].x;
    p0y = T[el][0][0][0].y;
    p0z = T[el][0][0][0].z;
    p1x = T[el][NY - 1][0][0].x;
    p1y = T[el][NY - 1][0][0].y;
    p1z = T[el][NY - 1][0][0].z;
    p2x = T[el][NY - 1][0][NZ - 1].x;
    p2y = T[el][NY - 1][0][NZ - 1].y;
    p2z = T[el][NY - 1][0][NZ - 1].z;
    p3x = T[el][0][0][NZ - 1].x;
    p3y = T[el][0][0][NZ - 1].y;
    p3z = T[el][0][0][NZ - 1].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx3, cx3, &sx);
    ncy = fitSurface (Cy3, cy3, &sy);
    ncz = fitSurface (Cz3, cz3, &sz);

    srx = surfaceResiduals ('x', f3n, &sx);
    sry = surfaceResiduals ('y', f3n, &sy);
    srz = surfaceResiduals ('z', f3n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f3n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f3n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f3n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx3_s, cx3_s, &sx);
      ncy = fitSimpleSurface (Cy3_s, cy3_s, &sy);
      ncz = fitSimpleSurface (Cz3_s, cz3_s, &sz);

      srx = surfaceResiduals ('x', f3n, &sx);
      sry = surfaceResiduals ('y', f3n, &sy);
      srz = surfaceResiduals ('z', f3n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f3.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f3e1, &f[el].f3.edges.e1, cpf3, f[el].f3.o)) return 1;
    if (computeEdges (f3e2, &f[el].f3.edges.e2, cpf3, f[el].f3.o)) return 1;
    if (computeEdges (f3e3, &f[el].f3.edges.e3, cpf3, f[el].f3.o)) return 1;
    if (computeEdges (f3e4, &f[el].f3.edges.e4, cpf3, f[el].f3.o)) return 1;

    cpe1 = f3e1[1];
    cpe2 = f3e2[1];
    cpe3 = f3e3[1];
    cpe4 = f3e4[1];

    switch (f[el].f3.o)
    {
      case 'x':
        f[el].f3.inner = sx;
        f[el].f3.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx3_p, cx3_p, &f[el].f3.outer);
        f[el].f3.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f3.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f3.inner = sy;
        f[el].f3.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy3_p, cy3_p, &f[el].f3.outer);
        f[el].f3.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f3.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f3.inner = sz;
        f[el].f3.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz3_p, cz3_p, &f[el].f3.outer);
        f[el].f3.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f3.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }

    p0x = T[el][0][NY - 1][0].x;
    p0y = T[el][0][NY - 1][0].y;
    p0z = T[el][0][NY - 1][0].z;
    p1x = T[el][NY - 1][NY - 1][0].x;
    p1y = T[el][NY - 1][NY - 1][0].y;
    p1z = T[el][NY - 1][NY - 1][0].z;
    p2x = T[el][NY - 1][NY - 1][NZ - 1].x;
    p2y = T[el][NY - 1][NY - 1][NZ - 1].y;
    p2z = T[el][NY - 1][NY - 1][NZ - 1].z;
    p3x = T[el][0][NY - 1][NZ - 1].x;
    p3y = T[el][0][NY - 1][NZ - 1].y;
    p3z = T[el][0][NY - 1][NZ - 1].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx4, cx4, &sx);
    ncy = fitSurface (Cy4, cy4, &sy);
    ncz = fitSurface (Cz4, cz4, &sz);

    srx = surfaceResiduals ('x', f4n, &sx);
    sry = surfaceResiduals ('y', f4n, &sy);
    srz = surfaceResiduals ('z', f4n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f4n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f4n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f4n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx4_s, cx4_s, &sx);
      ncy = fitSimpleSurface (Cy4_s, cy4_s, &sy);
      ncz = fitSimpleSurface (Cz4_s, cz4_s, &sz);

      srx = surfaceResiduals ('x', f4n, &sx);
      sry = surfaceResiduals ('y', f4n, &sy);
      srz = surfaceResiduals ('z', f4n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f4.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f4e1, &f[el].f4.edges.e1, cpf4, f[el].f4.o)) return 1;
    if (computeEdges (f4e2, &f[el].f4.edges.e2, cpf4, f[el].f4.o)) return 1;
    if (computeEdges (f4e3, &f[el].f4.edges.e3, cpf4, f[el].f4.o)) return 1;
    if (computeEdges (f4e4, &f[el].f4.edges.e4, cpf4, f[el].f4.o)) return 1;

    cpe1 = f4e1[1];
    cpe2 = f4e2[1];
    cpe3 = f4e3[1];
    cpe4 = f4e4[1];

    switch (f[el].f4.o)
    {
      case 'x':
        f[el].f4.inner = sx;
        f[el].f4.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx4_p, cx4_p, &f[el].f4.outer);
        f[el].f4.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f4.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f4.inner = sy;
        f[el].f4.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy4_p, cy4_p, &f[el].f4.outer);
        f[el].f4.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f4.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f4.inner = sz;
        f[el].f4.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz4_p, cz4_p, &f[el].f4.outer);
        f[el].f4.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f4.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }

    p0x = T[el][0][0][0].x;
    p0y = T[el][0][0][0].y;
    p0z = T[el][0][0][0].z;
    p1x = T[el][NX - 1][0][0].x;
    p1y = T[el][NX - 1][0][0].y;
    p1z = T[el][NX - 1][0][0].z;
    p2x = T[el][NX - 1][NY - 1][0].x;
    p2y = T[el][NX - 1][NY - 1][0].y;
    p2z = T[el][NX - 1][NY - 1][0].z;
    p3x = T[el][0][NY - 1][0].x;
    p3y = T[el][0][NY - 1][0].y;
    p3z = T[el][0][NY - 1][0].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx5, cx5, &sx);
    ncy = fitSurface (Cy5, cy5, &sy);
    ncz = fitSurface (Cz5, cz5, &sz);

    srx = surfaceResiduals ('x', f5n, &sx);
    sry = surfaceResiduals ('y', f5n, &sy);
    srz = surfaceResiduals ('z', f5n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f5n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f5n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f5n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx5_s, cx5_s, &sx);
      ncy = fitSimpleSurface (Cy5_s, cy5_s, &sy);
      ncz = fitSimpleSurface (Cz5_s, cz5_s, &sz);

      srx = surfaceResiduals ('x', f5n, &sx);
      sry = surfaceResiduals ('y', f5n, &sy);
      srz = surfaceResiduals ('z', f5n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f5.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f5e1, &f[el].f5.edges.e1, cpf5, f[el].f5.o)) return 1;
    if (computeEdges (f5e2, &f[el].f5.edges.e2, cpf5, f[el].f5.o)) return 1;
    if (computeEdges (f5e3, &f[el].f5.edges.e3, cpf5, f[el].f5.o)) return 1;
    if (computeEdges (f5e4, &f[el].f5.edges.e4, cpf5, f[el].f5.o)) return 1;

    cpe1 = f5e1[1];
    cpe2 = f5e2[1];
    cpe3 = f5e3[1];
    cpe4 = f5e4[1];

    switch (f[el].f5.o)
    {
      case 'x':
        f[el].f5.inner = sx;
        f[el].f5.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx5_p, cx5_p, &f[el].f5.outer);
        f[el].f5.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f5.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f5.inner = sy;
        f[el].f5.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy5_p, cy5_p, &f[el].f5.outer);
        f[el].f5.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f5.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f5.inner = sz;
        f[el].f5.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz5_p, cz5_p, &f[el].f5.outer);
        f[el].f5.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f5.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }

    p0x = T[el][0][0][NZ - 1].x;
    p0y = T[el][0][0][NZ - 1].y;
    p0z = T[el][0][0][NZ - 1].z;
    p1x = T[el][NX - 1][0][NZ - 1].x;
    p1y = T[el][NX - 1][0][NZ - 1].y;
    p1z = T[el][NX - 1][0][NZ - 1].z;
    p2x = T[el][NX - 1][NY - 1][NZ - 1].x;
    p2y = T[el][NX - 1][NY - 1][NZ - 1].y;
    p2z = T[el][NX - 1][NY - 1][NZ - 1].z;
    p3x = T[el][0][NY - 1][NZ - 1].x;
    p3y = T[el][0][NY - 1][NZ - 1].y;
    p3z = T[el][0][NY - 1][NZ - 1].z;

    sax = fabsl (crossProduct (p1y - p0y, p1z - p0z, p3y - p0y, p3z - p0z))
        + fabsl (crossProduct (p1y - p2y, p1z - p2z, p3y - p2y, p3z - p2z));

    say = fabsl (crossProduct (p1x - p0x, p1z - p0z, p3x - p0x, p3z - p0z))
        + fabsl (crossProduct (p1x - p2x, p1z - p2z, p3x - p2x, p3z - p2z));

    saz = fabsl (crossProduct (p1x - p0x, p1y - p0y, p3x - p0x, p3y - p0y))
        + fabsl (crossProduct (p1x - p2x, p1y - p2y, p3x - p2x, p3y - p2y));

    ncx = fitSurface (Cx6, cx6, &sx);
    ncy = fitSurface (Cy6, cy6, &sy);
    ncz = fitSurface (Cz6, cz6, &sz);

    srx = surfaceResiduals ('x', f6n, &sx);
    sry = surfaceResiduals ('y', f6n, &sy);
    srz = surfaceResiduals ('z', f6n, &sz);

    scx = ncx <= MAX_NORM_1 ? maxCurvature ('x', f6n, &sx) : INFINITY;
    scy = ncy <= MAX_NORM_1 ? maxCurvature ('y', f6n, &sy) : INFINITY;
    scz = ncz <= MAX_NORM_1 ? maxCurvature ('z', f6n, &sz) : INFINITY;

    if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
        min (srx, sry, srz) > MAX_NORM_2 ||
        min (scx, scy, scz) > MAX_CURVATURE)
    {
      ncx = fitSimpleSurface (Cx6_s, cx6_s, &sx);
      ncy = fitSimpleSurface (Cy6_s, cy6_s, &sy);
      ncz = fitSimpleSurface (Cz6_s, cz6_s, &sz);

      srx = surfaceResiduals ('x', f6n, &sx);
      sry = surfaceResiduals ('y', f6n, &sy);
      srz = surfaceResiduals ('z', f6n, &sz);

      if (min (ncx, ncy, ncz) > MAX_NORM_1 ||
          min (srx, sry, srz) > MAX_NORM_2) return 1;
    }

    f[el].f6.o = orientation (sax, say, saz,
                              ncx, ncy, ncz,
                              srx, sry, srz);

    if (computeEdges (f6e1, &f[el].f6.edges.e1, cpf6, f[el].f6.o)) return 1;
    if (computeEdges (f6e2, &f[el].f6.edges.e2, cpf6, f[el].f6.o)) return 1;
    if (computeEdges (f6e3, &f[el].f6.edges.e3, cpf6, f[el].f6.o)) return 1;
    if (computeEdges (f6e4, &f[el].f6.edges.e4, cpf6, f[el].f6.o)) return 1;

    cpe1 = f6e1[1];
    cpe2 = f6e2[1];
    cpe3 = f6e3[1];
    cpe4 = f6e4[1];

    switch (f[el].f6.o)
    {
      case 'x':
        f[el].f6.inner = sx;
        f[el].f6.p = position (cp.x, face (cp.y, cp.z, &sx));
        fitPlaneSurface (Cx6_p, cx6_p, &f[el].f6.outer);
        f[el].f6.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.y, cpe2.y, cpe1.z, cpe2.z);
        f[el].f6.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.y, cpe4.y, cpe3.z, cpe4.z);
      break;

      case 'y':
        f[el].f6.inner = sy;
        f[el].f6.p = position (cp.y, face (cp.x, cp.z, &sy));
        fitPlaneSurface (Cy6_p, cy6_p, &f[el].f6.outer);
        f[el].f6.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.z, cpe2.z);
        f[el].f6.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.z, cpe4.z);
      break;

      case 'z':
        f[el].f6.inner = sz;
        f[el].f6.p = position (cp.z, face (cp.x, cp.y, &sz));
        fitPlaneSurface (Cz6_p, cz6_p, &f[el].f6.outer);
        f[el].f6.edges.ti = BOUNDARY_RATIO_3
                          * distance2D (cpe1.x, cpe2.x, cpe1.y, cpe2.y);
        f[el].f6.edges.tj = BOUNDARY_RATIO_3
                          * distance2D (cpe3.x, cpe4.x, cpe3.y, cpe4.y);
      break;
    }
  }

  return 0;
}

unsigned facesCorrections (unsigned nel, struct Point T[nel][NX][NY][NZ],
                           struct Faces f[nel])
{
  /* Computes faces corrections for all faces of all spectral elements */
  for (unsigned el = 0; el < nel; el++)
  {
    unsigned n;

    long double sgm, max_sgm, u, v, w, s_w, rsd;

    max_sgm = MAX_SGM_RATIO * f[el].ti / BOUNDARY_RATIO_2;

    n = 0; sgm = 0;

    switch (f[el].f1.o)
    {
      case 'x':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][0][j][k].y;
            v = T[el][0][j][k].z;
            w = T[el][0][j][k].x;

            s_w = face (u, v, &f[el].f1.inner);

            rsd = w - s_w;

            f[el].f1.Fn[n].u = u;
            f[el].f1.Fn[n].v = v;
            f[el].f1.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][0][j][k].x;
            v = T[el][0][j][k].z;
            w = T[el][0][j][k].y;

            s_w = face (u, v, &f[el].f1.inner);

            rsd = w - s_w;

            f[el].f1.Fn[n].u = u;
            f[el].f1.Fn[n].v = v;
            f[el].f1.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][0][j][k].x;
            v = T[el][0][j][k].y;
            w = T[el][0][j][k].z;

            s_w = face (u, v, &f[el].f1.inner);

            rsd = w - s_w;

            f[el].f1.Fn[n].u = u;
            f[el].f1.Fn[n].v = v;
            f[el].f1.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f1.sigma = sgm < max_sgm ? sgm : max_sgm;

    n = 0; sgm = 0;

    switch (f[el].f2.o)
    {
      case 'x':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][NX - 1][j][k].y;
            v = T[el][NX - 1][j][k].z;
            w = T[el][NX - 1][j][k].x;

            s_w = face (u, v, &f[el].f2.inner);

            rsd = w - s_w;

            f[el].f2.Fn[n].u = u;
            f[el].f2.Fn[n].v = v;
            f[el].f2.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][NX - 1][j][k].x;
            v = T[el][NX - 1][j][k].z;
            w = T[el][NX - 1][j][k].y;

            s_w = face (u, v, &f[el].f2.inner);

            rsd = w - s_w;

            f[el].f2.Fn[n].u = u;
            f[el].f2.Fn[n].v = v;
            f[el].f2.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][NX - 1][j][k].x;
            v = T[el][NX - 1][j][k].y;
            w = T[el][NX - 1][j][k].z;

            s_w = face (u, v, &f[el].f2.inner);

            rsd = w - s_w;

            f[el].f2.Fn[n].u = u;
            f[el].f2.Fn[n].v = v;
            f[el].f2.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f2.sigma = sgm < max_sgm ? sgm : max_sgm;

    max_sgm = MAX_SGM_RATIO * f[el].tj / BOUNDARY_RATIO_2;

    n = 0; sgm = 0;

    switch (f[el].f3.o)
    {
      case 'x':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][0][k].y;
            v = T[el][i][0][k].z;
            w = T[el][i][0][k].x;

            s_w = face (u, v, &f[el].f3.inner);

            rsd = w - s_w;

            f[el].f3.Fn[n].u = u;
            f[el].f3.Fn[n].v = v;
            f[el].f3.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][0][k].x;
            v = T[el][i][0][k].z;
            w = T[el][i][0][k].y;

            s_w = face (u, v, &f[el].f3.inner);

            rsd = w - s_w;

            f[el].f3.Fn[n].u = u;
            f[el].f3.Fn[n].v = v;
            f[el].f3.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][0][k].x;
            v = T[el][i][0][k].y;
            w = T[el][i][0][k].z;

            s_w = face (u, v, &f[el].f3.inner);

            rsd = w - s_w;

            f[el].f3.Fn[n].u = u;
            f[el].f3.Fn[n].v = v;
            f[el].f3.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f3.sigma = sgm < max_sgm ? sgm : max_sgm;

    n = 0; sgm = 0;

    switch (f[el].f4.o)
    {
      case 'x':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][NY - 1][k].y;
            v = T[el][i][NY - 1][k].z;
            w = T[el][i][NY - 1][k].x;

            s_w = face (u, v, &f[el].f4.inner);

            rsd = w - s_w;

            f[el].f4.Fn[n].u = u;
            f[el].f4.Fn[n].v = v;
            f[el].f4.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][NY - 1][k].x;
            v = T[el][i][NY - 1][k].z;
            w = T[el][i][NY - 1][k].y;

            s_w = face (u, v, &f[el].f4.inner);

            rsd = w - s_w;

            f[el].f4.Fn[n].u = u;
            f[el].f4.Fn[n].v = v;
            f[el].f4.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned k = 0; k < NZ; k++)
          {
            u = T[el][i][NY - 1][k].x;
            v = T[el][i][NY - 1][k].y;
            w = T[el][i][NY - 1][k].z;

            s_w = face (u, v, &f[el].f4.inner);

            rsd = w - s_w;

            f[el].f4.Fn[n].u = u;
            f[el].f4.Fn[n].v = v;
            f[el].f4.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
        break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f4.sigma = sgm < max_sgm ? sgm : max_sgm;

    max_sgm = MAX_SGM_RATIO * f[el].tk / BOUNDARY_RATIO_2;

    n = 0; sgm = 0;

    switch (f[el].f5.o)
    {
      case 'x':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][0].y;
            v = T[el][i][j][0].z;
            w = T[el][i][j][0].x;

            s_w = face (u, v, &f[el].f5.inner);

            rsd = w - s_w;

            f[el].f5.Fn[n].u = u;
            f[el].f5.Fn[n].v = v;
            f[el].f5.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][0].x;
            v = T[el][i][j][0].z;
            w = T[el][i][j][0].y;

            s_w = face (u, v, &f[el].f5.inner);

            rsd = w - s_w;

            f[el].f5.Fn[n].u = u;
            f[el].f5.Fn[n].v = v;
            f[el].f5.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][0].x;
            v = T[el][i][j][0].y;
            w = T[el][i][j][0].z;

            s_w = face (u, v, &f[el].f5.inner);

            rsd = w - s_w;

            f[el].f5.Fn[n].u = u;
            f[el].f5.Fn[n].v = v;
            f[el].f5.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f5.sigma = sgm < max_sgm ? sgm : max_sgm;

    n = 0; sgm = 0;

    switch (f[el].f6.o)
    {
      case 'x':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][NZ - 1].y;
            v = T[el][i][j][NZ - 1].z;
            w = T[el][i][j][NZ - 1].x;

            s_w = face (u, v, &f[el].f6.inner);

            rsd = w - s_w;

            f[el].f6.Fn[n].u = u;
            f[el].f6.Fn[n].v = v;
            f[el].f6.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'y':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][NZ - 1].x;
            v = T[el][i][j][NZ - 1].z;
            w = T[el][i][j][NZ - 1].y;

            s_w = face (u, v, &f[el].f6.inner);

            rsd = w - s_w;

            f[el].f6.Fn[n].u = u;
            f[el].f6.Fn[n].v = v;
            f[el].f6.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;

      case 'z':
        for (unsigned i = 0; i < NX; i++)

          for (unsigned j = 0; j < NY; j++)
          {
            u = T[el][i][j][NZ - 1].x;
            v = T[el][i][j][NZ - 1].y;
            w = T[el][i][j][NZ - 1].z;

            s_w = face (u, v, &f[el].f6.inner);

            rsd = w - s_w;

            f[el].f6.Fn[n].u = u;
            f[el].f6.Fn[n].v = v;
            f[el].f6.Fn[n].w = rsd;

            sgm += squarel (rsd); n++;
          }
      break;
    }

    sgm = sqrt (sgm / n); if (sgm > 3.L * max_sgm) return 1;

    f[el].f6.sigma = sgm < max_sgm ? sgm : max_sgm;
  }

  return 0;
}

void detectBoundaries (double rmin, double rmax,
                       struct SphericalPoint Tm[NEL][NX][NY][NZ],
                       bool In[NEL], double *r1, double *r2)
{
  /* Detect boundaries of the domain */
  unsigned hnx = NX / 2, hny = NY / 2, hnz = NZ / 2;

  *r1 =  INFINITY;
  *r2 = -INFINITY;

  for (unsigned el = 0; el < NEL; el++)
  {
    In[el] = false;

    double rc = Tm[el][hnx][hny][hnz].r;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          double rp = Tm[el][i][j][k].r;

          if (rc >= rmin && rc <= rmax)
          {
            In[el] = true;

            if (rp < *r1) *r1 = rp;
            if (rp > *r2) *r2 = rp;
          }
        }
  }
}

