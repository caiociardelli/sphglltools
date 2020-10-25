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

#ifndef STRETCH_H
#define STRETCH_H
void computeStretchingArrays (unsigned nel, unsigned *nelic, unsigned *neloc,
                              unsigned mic[nel], unsigned moc[NEL2],
                              bool uic[nel], bool uoc[NEL2],
                              struct Point Ti[nel][NX][NY][NZ],
                              struct Point To[NEL2][NX][NY][NZ]);

void stretchUpperCrust (unsigned nel, unsigned mic[nel], unsigned moc[NEL2],
                        unsigned nelic, unsigned imic[nelic], unsigned neloc,
                        bool uic[nel], bool uoc[NEL2], long double Z[NZ],
                        struct Point Ti[nel][NX][NY][NZ],
                        struct Point To[NEL2][NX][NY][NZ],
                        struct Point Tic[nelic][NX][NY][NZ],
                        struct Point Toc[neloc][NX][NY][NZ]);
#endif
