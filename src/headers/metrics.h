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

#ifndef METRICS_H
#define METRICS_H
static inline long double norm (long double x, long double y, long double z)
{
  /* Computes the L2 norm in 3D */
  return sqrtl (squarel (x) + squarel (y) + squarel (z));
}

static inline long double distance1D (long double u0, long double u1)
{
  /* Computes the distance in 1D */
  return fabsl (u1 - u0);
}

static inline long double distance2D (long double u0, long double u1,
                                      long double v0, long double v1)
{
  /* Computes the distance in 2D */
  return sqrtl (squarel (u1 - u0) + squarel (v1 - v0));
}

static inline long double distance3D (struct Point *p2, struct Point *p3)
{
  /* Computes the distance in 3D */
  return sqrtl (squarel (p2->x - p3->x) +
                squarel (p2->y - p3->y) +
                squarel (p2->z - p3->z));
}

static inline long double squaredDistance2D (long double u0, long double u1,
                                             long double v0, long double v1)
{
  /* Computes the squared distance in 2D */
  return squarel (u1 - u0) + squarel (v1 - v0);
}

static inline long double squaredDistance3D (struct Point *p2, struct Point *p3)
{
  /* Computes the squared distance in 3D */
  return squarel (p2->x - p3->x) + squarel (p2->y - p3->y) + squarel (p2->z - p3->z);
}
#endif
