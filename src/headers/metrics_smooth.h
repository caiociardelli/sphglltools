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

#include "structs_smooth.h"

#ifndef METRICS_H
#define METRICS_H
static inline double distance (struct Point *p1, struct Point *p2)
{
  return sqrt (square (p1->x - p2->x) +
               square (p1->y - p2->y) +
               square (p1->z - p2->z));
}

static inline double radialDistance (struct Point *p1, struct Point *p2)
{
  return fabs (p1->r - p2->r);
}
#endif
