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
#include "structs_smooth.h"
#include "metrics_smooth.h"

void identifyZones (unsigned n, unsigned l[n],
                    struct Point T[n][NX][NY][NZ])
{
  /* Identifies zone based on the r coordinate */
  unsigned hx = NX / 2, hy = NY / 2, hz = NZ / 2;

  for (unsigned i = 0; i < n; i++)
  {
    double r = T[i][hx][hy][hz].r;

    if (r > D410_R) l[i] = 1;

    else if (r > D650_R) l[i] = 2;

    else l[i] = 3;
  }
}

void computeElRadii (unsigned nel,
                     struct Point Ti[nel][NX][NY][NZ],
                     struct Point To[NEL][NX][NY][NZ],
                     struct Radii imr[nel],
                     struct Radii omr[NEL])
{
  /* Computes 3D distances and radial distances between
     point p and the center of all spectral elements */
  unsigned hx = NX / 2, hy = NY / 2, hz = NZ / 2;

  for (unsigned el = 0; el < nel; el++)
  {
    struct Point *c = &Ti[el][hx][hy][hz];

    imr[el].r = 0; imr[el].rv = 0;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          struct Point *p = &Ti[el][i][j][k];

          double r  = distance (c, p);
          double rv = radialDistance (c, p);

          if (r  > imr[el].r)  imr[el].r  = r;
          if (rv > imr[el].rv) imr[el].rv = rv;
        }
  }

  for (unsigned el = 0; el < NEL; el++)
  {
    struct Point *c = &To[el][hx][hy][hz];

    omr[el].r = 0; omr[el].rv = 0;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          struct Point *p = &To[el][i][j][k];

          double r  = distance (c, p);
          double rv = radialDistance (c, p);

          if (r  > omr[el].r)  omr[el].r  = r;
          if (rv > omr[el].rv) omr[el].rv = rv;
        }
  }
}

