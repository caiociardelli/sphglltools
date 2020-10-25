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

#ifndef COORDINATES_H
#define COORDINATES_H
static inline double rad2Degree (double angle)
{
  /* Converts radians to degrees */
  return angle * TO_DEGREE;
}

static inline long double rad2Degreel (long double angle)
{
  /* Converts radians to degrees */
  return angle * TO_DEGREE;
}

static inline double degree2Rad (double v)
{
  /* Converts degrees to radians */
  return v * TO_RADIANS;
}

static inline long double degree2Radl (long double v)
{
  /* Converts degrees to radians */
  return v * TO_RADIANS;
}

static inline double r2Depth (double r)
{
  /* Converts Earth's normalized radius to depth in km */
  return (1.0 - r) * EARTH_R;
}

static inline long double r2Depthl (long double r)
{
  /* Converts Earth's normalized radius to depth in km */
  return (1.L - r) * EARTH_R;
}

static inline double depth2R (double depth)
{
  /* Converts depth in km to Earth's normalized radius */
    return 1.0 - (depth / EARTH_R);
}

static inline long double depth2Rl (long double depth)
{
  /* Converts depth in km to Earth's normalized radius */
    return 1.L - (depth / EARTH_R);
}

void rThetaPhi2XYZ (double r, double theta, double phi,
                    double *x, double *y, double *z);
void rThetaPhi2XYZl (long double r, long double theta, long double phi,
                     long double *x, long double *y, long double *z);

void xYZ2RThetaPhi (double x, double y, double z,
                    double *r, double *theta, double *phi);
void xYZ2RThetaPhil (long double x, long double y, long double z,
                     long double *r, long double *theta, long double *phi);

double vincenty (double t1, double p1, double t2, double p2);
long double vincentyl (long double t1, long double p1,
                       long double t2, long double p2);
#endif
