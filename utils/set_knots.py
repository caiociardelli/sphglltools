#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
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

-----------------------------------------------------------------------------------------------

 SET_KNOTS

 USAGE
   ./utils/set_knots.py PARAMETER

 EXAMPLE
   ./utils/set_knots.py vs

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter used to fit the B-splines (vp, vs, rho, eta, vsv, etc.)

 DESCRIPTION
   Creates the files with the positions of the accurate knots, used to recreate the model from
   the spherical harmonics coefficients.

----------------------------------------------------------------------------------------------- */
"""

from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import sys
import numpy as np

def readKnots (name):
  """
  Reads knots
  """
  with open (name, 'r') as FILE:

    lines = FILE.readlines ()[1:]

    i = 0; t = np.array ([])

    while lines[i][:-1] != 'Coefficients':

      knot = float (lines[i][:-1])

      t = np.append (t, np.array (knot))

      i += 1

  return t


def writeAccurateKnots (rmin, rmax, zone):
  """
  Creates the knots files with accurate
  minimum and maximum radii.
  """
  output = str ()

  with open ('knots_Zone{}.dat'.format (zone), 'r') as FILE:

    i = 0

    for line in FILE:

      if i == 3:

        output += '{:.7f}\n'.format (rmin)

      elif i == 4:

        output += '{:.7f}\n'.format (rmax)

      else:

        output += line

      i += 1

  with open ('knots_Z{}.dat'.format (zone), 'w') as FILE:

    FILE.write (output[:-1])


def helpMenu ():
  """
  Prints help menu
  """
  help = """\n Error: wrong number of parameters on the comand line...

 SET_KNOTS

 USAGE
   ./utils/set_knots.py PARAMETER

 EXAMPLE
   ./utils/set_knots.py vs

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter used to fit the B-splines (vp, vs, rho, eta, vsv, etc.)

 DESCRIPTION
   Creates the files with the positions of the accurate knots, used to recreate the model from
   the spherical harmonics coefficients.\n"""

  print (help)


if __name__ == '__main__':

  if len (sys.argv) != 2:

    sys.exit (helpMenu ())

  prm = sys.argv[1]

  t2 = readKnots ('{}_Z2_KC.dat'.format (prm))
  t3 = readKnots ('{}_Z3_KC.dat'.format (prm))
  t4 = readKnots ('{}_Z4_KC.dat'.format (prm))

  r2min, r2max = t2[0], t2[-1]
  r3min, r3max = t3[0], t3[-1]
  r4min, r4max = t4[0], t4[-1]

  r2minf = 0.5 * (r2min + r3max)
  r2maxf = r2max + 0.5 * (r2min - r3max)
  r3minf = 0.5 * (r3min + r4max)
  r3maxf = r2minf
  r4minf = r4min - 0.5 * (r3min - r4max)
  r4maxf = r3minf

  writeAccurateKnots (r2minf, r2maxf, 2)
  writeAccurateKnots (r3minf, r3maxf, 3)
  writeAccurateKnots (r4minf, r4maxf, 4)

