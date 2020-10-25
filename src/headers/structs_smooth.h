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

#include "config_smooth.h"
#include "constants.h"

#ifndef STRUCTS_SMOOTH_H
#define STRUCTS_SMOOTH_H
struct Point
{
  /* Used to store the Cartesian and
     spherical coordinates of a point */
  double x;
  double y;
  double z;
  double r;
  double theta;
  double phi;
};

struct Parameters
{
  /* Used to store the model parameters */
  #if VP
  double vp;
  #endif

  #if VS
  double vs;
  #endif

  #if RHO
  double rho;
  #endif

  #if VPV
  double vpv;
  #endif

  #if VPH
  double vph;
  #endif

  #if VSV
  double vsv;
  #endif

  #if VSH
  double vsh;
  #endif

  #if ETA
  double eta;
  #endif

  #if QMU
  double qmu;
  #endif

  #if APK
  double apk;
  #endif

  #if BTK
  double btk;
  #endif

  #if RHK
  double rhk;
  #endif

  #if BCK
  double bck;
  #endif

  #if BBK
  double bbk;
  #endif

  #if BVK
  double bvk;
  #endif

  #if BHK
  double bhk;
  #endif

  #if ETK
  double etk;
  #endif

  #if HSK
  double hsk;
  #endif
};

struct Boundaries
{
  /* Used to store the boundaries
     of a spectral element */
  double rmin;
  double rmax;
  double tmin;
  double tmax;
  double pmin;
  double pmax;
};

struct Radii
{
  /* Used to store the radiuses */
  double r;
  double rv;
};

struct Profile
{
  /* Used to store the information
     of the the profile */
  double r;
  double sigma_h;
  double sigma_v;
  double search_h;
  double search_v;
};

struct Smooth
{
  /* Used to store the information
     for smoothing the parameteres */
  double sigma_h;
  double sigma_v;
  double search_h;
  double search_v;
  double sg_h;
  double sg_v;
  double sc_h;
  double sc_v;
};

struct Indexes
{
  /* Used to store the indexes
     of a spectral element */
  unsigned element;
  unsigned i;
  unsigned j;
  unsigned k;
};

struct ElNode
{
  /* Cell to store the coordinates
     of a spectral element */
  struct Point Ti[M_NX][M_NY][M_NZ];
  struct ElNode *next;
};

struct PmNode
{
  /* Cell to store the model
     parameters */
  struct Parameters M[M_NX][M_NY][M_NZ];
  struct PmNode *next;
};
#endif
