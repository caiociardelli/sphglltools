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

#include <math.h>
#include "exmath.h"
#include "constants.h"

void rThetaPhi2XYZ (double r, double theta, double phi,
                    double *x, double *y, double *z)
{
  /* Converts spherical coordinates to Cartesian */
  *x = r * sin (theta) * cos (phi);
  *y = r * sin (theta) * sin (phi);
  *z = r * cos (theta);
}

void rThetaPhi2XYZl (long double r, long double theta, long double phi,
                     long double *x, long double *y, long double *z)
{
  /* Converts spherical coordinates to Cartesian */
  *x = r * sinl (theta) * cosl (phi);
  *y = r * sinl (theta) * sinl (phi);
  *z = r * cosl (theta);
}

void xYZ2RThetaPhi (double x, double y, double z,
                    double *r, double *theta, double *phi)
{
  /* Converts Cartesian coordinates to spherical */
  *r = sqrt (square (x) + square (y) + square (z));

  if (z > -ANGLE_TOLERANCE && z <= 0.0) z = -ANGLE_TOLERANCE;
  if (z <  ANGLE_TOLERANCE && z >= 0.0) z =  ANGLE_TOLERANCE;

  *theta = atan2 (sqrt (square (x) + square (y)), z);

  if (x > -ANGLE_TOLERANCE && x <= 0.0) x = -ANGLE_TOLERANCE;
  if (x <  ANGLE_TOLERANCE && x >= 0.0) x =  ANGLE_TOLERANCE;

  *phi = atan2 (y, x);
}

void xYZ2RThetaPhil (long double x, long double y, long double z,
                     long double *r, long double *theta, long double *phi)
{
  /* Converts Cartesian coordinates to spherical */
  *r = sqrtl (squarel (x) + squarel (y) + squarel (z));

  if (z > -ANGLE_TOLERANCE && z <= 0.L) z = -ANGLE_TOLERANCE;
  if (z <  ANGLE_TOLERANCE && z >= 0.L) z =  ANGLE_TOLERANCE;

  *theta = atan2l (sqrtl (squarel (x) + squarel (y)), z);

  if (x > -ANGLE_TOLERANCE && x <= 0.L) x = -ANGLE_TOLERANCE;
  if (x <  ANGLE_TOLERANCE && x >= 0.L) x =  ANGLE_TOLERANCE;

  *phi = atan2l (y, x);
}

double vincenty (double t1, double p1, double t2, double p2)
{
  /* Computes the great circle distance using Vincenty's formula */
  double sin_t1 = sin (t1);
  double cos_t1 = cos (t1);
  double sin_t2 = sin (t2);
  double cos_t2 = cos (t2);

  double sin_dp = sin (p2 - p1);
  double cos_dp = cos (p2 - p1);

  return atan2 (sqrt (square (sin_t2 * sin_dp)
                      + square (sin_t1 * cos_t2 - cos_t1 * sin_t2 * cos_dp)),
                cos_t1 * cos_t2 + sin_t1 * sin_t2 * cos_dp);
}

long double vincentyl (long double t1, long double p1,
                       long double t2, long double p2)
{
  /* Computes the great circle distance using Vincenty's formula */
  long double sin_t1 = sinl (t1);
  long double cos_t1 = cosl (t1);
  long double sin_t2 = sinl (t2);
  long double cos_t2 = cosl (t2);

  long double sin_dp = sinl (p2 - p1);
  long double cos_dp = cosl (p2 - p1);

  return atan2l (sqrtl (squarel (sin_t2 * sin_dp)
                        + squarel (sin_t1 * cos_t2 - cos_t1 * sin_t2 * cos_dp)),
                 cos_t1 * cos_t2 + sin_t1 * sin_t2 * cos_dp);
}

