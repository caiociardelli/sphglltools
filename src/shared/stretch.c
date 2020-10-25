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
#include "structs.h"
#include "coordinates.h"
#include "constants.h"
#include "config.h"

void computeStretchingArrays (unsigned nel, unsigned *nelic, unsigned *neloc,
                              unsigned mic[nel], unsigned moc[NEL2],
                              bool uic[nel], bool uoc[NEL2],
                              struct Point Ti[nel][NX][NY][NZ],
                              struct Point To[NEL2][NX][NY][NZ])
{
  /* Computes stretching for the spectral elements */
  unsigned elm = 0, hnx = NX / 2, hny = NY / 2, tnz = NZ - 1;

  for (unsigned el = 0; el < nel; el++)
  {
    if (Ti[el][hnx][hny][tnz].r > MIN_SURFACE_R)
    {
      mic[el] = elm++; uic[el] = true;
    }

    else uic[el] = false;
  }

  *nelic = elm; elm = 0;

  for (unsigned el = 0; el < NEL2; el++)
  {
    if (To[el][hnx][hny][tnz].r > MIN_SURFACE_R)
    {
      moc[el] = elm++; uoc[el] = true;
    }

    else uoc[el] = false;
  }

  *neloc = elm;
}

void stretchUpperCrust (unsigned nel, unsigned mic[nel], unsigned moc[NEL2],
                        unsigned nelic, unsigned imic[nelic], unsigned neloc,
                        bool uic[nel], bool uoc[NEL2], long double Z[NZ],
                        struct Point Ti[nel][NX][NY][NZ],
                        struct Point To[NEL2][NX][NY][NZ],
                        struct Point Tic[nelic][NX][NY][NZ],
                        struct Point Toc[neloc][NX][NY][NZ])
{
  /* Stretches the first layer of spectral elements to turn it into a perfectly
     spherical shell (to avoid interpolation issues caused by the topography) */
  long double m = 0.5L * (MAX_STC_UPPER_R - MIN_STC_UPPER_R);
  long double n = MIN_STC_UPPER_R - m * Z[0];

  for (unsigned elm = 0, el = 0; el < nel; el++)
  {
    if (uic[el])
    {
      imic[elm++] = el;

      for (unsigned i = 0; i < NX; i++)

        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            struct Point p = Ti[el][i][j][k];

            p.r = m * Z[k] + n;

            rThetaPhi2XYZl (p.r, p.theta, p.phi, &p.x, &p.y, &p.z);

            Tic[mic[el]][i][j][k] = p;
          }
    }
  }

  for (unsigned el = 0; el < NEL2; el++)

    if (uoc[el])

      for (unsigned i = 0; i < NX; i++)

        for (unsigned j = 0; j < NY; j++)

          for (unsigned k = 0; k < NZ; k++)
          {
            struct Point p = To[el][i][j][k];

            p.r = m * Z[k] + n;

            rThetaPhi2XYZl (p.r, p.theta, p.phi, &p.x, &p.y, &p.z);

            Toc[moc[el]][i][j][k] = p;
          }
}

