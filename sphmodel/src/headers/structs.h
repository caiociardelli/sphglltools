/*
 SphGLLTools/SphModel

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

#ifndef STRUCTS_H
#define STRUCTS_H
struct SphericalBoundaries
{
  /* Used to store the boundaries
     of the domain */
  double rmin;
  double rmax;
  double tmin;
  double tmax;
  double pmin;
  double pmax;
};

struct MeanModel
{
  /* Used to store the radius
     and the corresponding
     model parameter for the
     mean model */
  double r;
  double vam;
  double vgm;
};
#endif
