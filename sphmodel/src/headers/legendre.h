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

#include "constants.h"

#ifndef LEGENDRE_H
#define LEGENDRE_H
static inline unsigned N2Nlg (unsigned N)
{
  /* Computes the total number of Associated 
     Legendre functions up to degree N */
  return N * (N + 1) / 2 + N + 1;
}

static inline unsigned mN2I (unsigned m, unsigned n)
{
  /* Maps m and n to index */
  return n * (n + 1) / 2 + m;
}

void nmlFactors (unsigned nmax, unsigned nlg, long double nF[nlg]);
void lgPmn (long double x, unsigned n, unsigned nmax, long double P[nmax + 2]);
void normalize (unsigned n, unsigned nlg, unsigned nmax,
                long double nF[nlg],
                long double P[nmax + 2], long double nP[nmax + 2]);
#endif
