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

#ifndef EXMATH_H
#define EXMATH_H
static inline double square (double v)
{
  /* Computes the square of v */
  return v * v;
}

static inline double cube (double v)
{
  /* Computes the cube of v */
  return v * v * v;
}

static inline long double squarel (long double v)
{
  /* Computes the square of v */
  return v * v;
}

static inline long double cubel (long double v)
{
  /* Computes the cube of v */
  return v * v * v;
}

static inline unsigned v2Index (double v, double v1, double dv_i)
{
  /* Maps value to index */
  return (unsigned) ((v - v1) * dv_i + 0.5);
}

static inline unsigned nnAndDg2Ns (unsigned nn, unsigned dg)
{
  /* Computes the total number of B-splines */
  return nn + dg - 1;
}

static inline unsigned nnAndDg2Nnt (unsigned nn, unsigned dg)
{
  /* Computes the total number of knots */
  return nn + 2 * dg;
}

long double powerl (long double v, unsigned n);

double roundV (double v);

long double factorial (unsigned n);
long double factorial2 (unsigned n);

static inline long double cNN (long double u, long double v,
                               unsigned i, unsigned j)
{
  /* Computes u^i * v^j */
  return powerl (u, i) * powerl (v, j);
}

static inline long double cNNN (long double u, long double v, long double w,
                                unsigned i, unsigned j, unsigned k)
{
  /* Computes u^i * v^j * w^k */
  return powerl (u, i) * powerl (v, j) * powerl (w, k);
}

unsigned gaussJordan (unsigned n, double M[n][n]);
unsigned gaussJordanl (unsigned n, long double M[n][n]);

double simpson (double x1, double x2, unsigned nx, double f[nx]);

double bSplines (double x, unsigned nnt, unsigned dg,
                 double T[nnt], unsigned i);
#endif
