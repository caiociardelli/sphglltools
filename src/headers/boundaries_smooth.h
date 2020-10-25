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

#ifndef BOUNDARIES_SMOOTH_H
#define BOUNDARIES_SMOOTH_H
void identifyZones (unsigned n, unsigned l[n],
                    struct Point T[n][NX][NY][NZ]);

void computeElRadii (unsigned nel,
                     struct Point Ti[nel][NX][NY][NZ],
                     struct Point To[NEL][NX][NY][NZ],
                     struct Radii imr[nel],
                     struct Radii omr[NEL]);
#endif
