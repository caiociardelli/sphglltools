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

#ifndef BOUNDARIES_H
#define BOUNDARIES_H
static inline long double edge (long double u, struct Curve *c)
{
  /* Computes the edge of a spectral element */
  return c->a + c->au * u + c->au2 * squarel (u);
}

static inline long double face (long double u, long double v, struct Surface *s)
{
  /* Computes the face of a spectral element */
  return s->a +
         s->au * u + s->av * v +
         s->auv * u * v + s->au2 * squarel (u) + s->av2 * squarel (v) +
         s->au2v * squarel (u) * v + s->auv2 * u * squarel (v) +
         s->au2v2 * squarel (u) * squarel (v);
}

static inline long double getCoordinate (struct Point *p, char o)
{
  /* Gets the coordinate according to the orientation */
  if (o == 'x') return p->x; else if (o == 'y') return p->y; else return p->z;
}

static inline bool isInEdgeShadow (long double u, struct Corners *c)
{
  /* Checks if a point is inside the edge shadow */
  return (u + c->ti < c->c1 || u - c->ti > c->c2) ? false : true;
}

void getBoundaries (unsigned nel, struct Point T[nel][NX][NY][NZ],
                    struct Boundaries b[nel]);

unsigned computeFaces (unsigned nel, struct Point T[nel][NX][NY][NZ],
                       struct Faces f[nel]);

unsigned facesCorrections (unsigned nel, struct Point T[nel][NX][NY][NZ],
                           struct Faces f[nel]);

void detectBoundaries (double rmin, double rmax,
                       struct SphericalPoint Tm[NEL][NX][NY][NZ],
                       bool In[NEL], double *r1, double *r2);
#endif
