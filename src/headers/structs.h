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

#ifndef STRUCTS_H
#define STRUCTS_H
struct CartesianPoint
{
  /* Used to store the Cartesian
     coordinates of a point */
  long double x;
  long double y;
  long double z;
};

struct SphericalPoint
{
  /* Used to store the spherical
     coordinates of a point */
  double r;
  double theta;
  double phi;
};

struct Point
{
  /* Used to store the Cartesian and
     spherical coordinates of a point */
  long double x;
  long double y;
  long double z;
  long double r;
  long double theta;
  long double phi;
};

struct Parameters
{
  /* Used to store the model parameters */
  #if VP
  long double vp;
  #endif

  #if VS
  long double vs;
  #endif

  #if RHO
  long double rho;
  #endif

  #if VPV
  long double vpv;
  #endif

  #if VPH
  long double vph;
  #endif

  #if VSV
  long double vsv;
  #endif

  #if VSH
  long double vsh;
  #endif

  #if ETA
  long double eta;
  #endif

  #if QMU
  long double qmu;
  #endif
};

struct Means
{
  /* Used to store the averages
     of the model parameters */
  struct Parameters a;
  struct Parameters g;
};

struct ElNode
{
  /* Cell to store the coordinates
     of a spectral element */
  struct Point T[M_NX][M_NY][M_NZ];
  struct ElNode *next;
};

struct PmNode
{
  /* Cell to store the model
     parameters */
  struct Parameters M[M_NX][M_NY][M_NZ];
  struct PmNode *next;
};

struct FaceNode
{
  /* Used to store the coordinates
     of a face node */
  long double u;
  long double v;
  long double w;
};

struct Boundaries
{
  /* Used to store the boundaries
     of a spectral element */
  long double xmin;
  long double xmax;
  long double ymin;
  long double ymax;
  long double zmin;
  long double zmax;
  long double rmin;
  long double rmax;
  long double tmin;
  long double tmax;
  long double pmin;
  long double pmax;
};

struct Corners
{
  /* Used to store the information
     of a pair of corners */
  char o;
  long double ti;
  long double c1;
  long double c2;
};

struct Curve
{
  /* Used to store the coefficients
     of a parabola */
  long double a;
  long double au;
  long double au2;
};

struct Edge
{
  /* Used to store the information
     of an edge */
  char o;
  char p;
  struct Curve inner;
  struct Curve outer;
  struct Corners corners;
};

struct Edges
{
  /* Used to store the information
     of the four edges of a face */
  long double ti;
  long double tj;
  struct Edge e1;
  struct Edge e2;
  struct Edge e3;
  struct Edge e4;
};

struct Surface
{
  /* Used to store the coefficients
     of a surface */
  long double a;
  long double au;
  long double au2;
  long double av;
  long double auv;
  long double au2v;
  long double av2;
  long double auv2;
  long double au2v2;
};

struct Face
{
  /* Used to store the information
     of a face */
  char o;
  char p;
  long double sigma;
  struct Surface inner;
  struct Surface outer;
  struct Edges edges;
  struct FaceNode Fn[M_NPf];
};

struct Faces
{
  /* Used to store the information
     of the six faces of a spectral
     element */
  long double ti;
  long double tj;
  long double tk;
  struct Face f1;
  struct Face f2;
  struct Face f3;
  struct Face f4;
  struct Face f5;
  struct Face f6;
};

struct Indexes
{
  /* Used to store the indexes
     of a spectral element and
     the flag to check if its
     outside the mesh */
  bool outside_mesh;
  unsigned element;
  unsigned i;
  unsigned j;
  unsigned k;
};

struct BasisIndices
{
  /* Used to store the indices
     of the basis functions */
  unsigned ri;
  unsigned ti;
  unsigned pi;
};

struct MeanModel
{
  /* Used to store the radius
     and the corresponding
     model parameters for the
     mean model */
  long double r;
  long double vam;
  long double vgm;
};
#endif
