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

 CREATE_VERTICAL_CROSS_SECTION

 USAGE
   ./utils/created_extra_vertical_cross_section.py PARAMETER

 EXAMPLE
   ./utils/created_extra_vertical_cross_section.py dcb

 COMMAND LINE ARGUMENTS
   PARAMETER              - derived model parameter to be created (cb, vpvs, ti, dcb, dvs, etc.)

 DESCRIPTION
   Calculates additional model parameters from the ones created by BKMNS2DD. You can calculate
   the isotropic velocities using the parameter 'iso' (or 'diso', in case you want the Vp and Vs
   perturbations). You will need the files for Vpv, Vph, Vsv, Vsh, and Eta, created at the same
   region and with the same resolution. To compute the bulk sound speed (also needs Vpv, Vph, Vsv,
   Vsh, and Eta), use parameter 'cb' (or 'dcb', in case you want the perturbations). The same
   applies to the Vp/Vs (use parameter 'vpvs'). To calculate the transverse isotropy from Vsv and
   Vsh, use the parameter 'ti'.

----------------------------------------------------------------------------------------------- */
"""

from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as floating_axes

from scipy.interpolate import interp1d
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter

EARTH_R = 6371.0


def meanModel (name):
  """
  Reads mean model
  """
  d = list ()
  v = list ()

  with open (name, 'r') as FILE:

    next (FILE)
    next (FILE)

    for line in FILE:

      data = line.split ()

      d += [float (data[0])]
      v += [float (data[2])]

  d = np.array (d)
  v = np.array (v)

  return interp1d (d, v)


def readGridFile (name):
  """
  Reads grid file
  """
  with open (name, 'r') as FILE:

    line = FILE.readline ().split ()

    nrad = int (line[2])
    ndel = int (line[3])

    gcline = line = FILE.readline ()

    radius = np.zeros (nrad * ndel)
    delta  = np.zeros (nrad * ndel)
    M      = np.zeros (nrad * ndel)

    next (FILE)

    i = 0

    for line in FILE:

      data = line.split ()

      radius[i] = float (data[0])
      delta[i]  = float (data[1])
      M[i]      = float (data[2])

      i += 1

  depth = EARTH_R - radius

  return depth, nrad, ndel, gcline, radius, delta, M


def writeOutput (name, label, depth, nrad, ndel,
                 gcline, radius, delta, M):
  """
  Writes grid file for the derived parameter
  """
  with open (name, 'w') as FILE:

    header  = '#nrad ndel: {} {}\n'.format (nrad, ndel)
    header += gcline
    header += '#radius (km)   delta (degrees)   {}\n'.format (label)

    FILE.write (header.encode ().decode ('utf-8'))

    for r, d, m in zip (radius, delta, M):

      FILE.write ('{:10.3f} {:14.3f} {:19E}\n'.format (r, d, m).encode ().decode ('utf-8'))


def helpMenu ():
  """
  Prints help menu
  """
  help = """\n Error: wrong number of parameters on the comand line...

 CREATE_VERTICAL_CROSS_SECTION

 USAGE
   ./utils/created_extra_vertical_cross_section.py PARAMETER

 EXAMPLE
   ./utils/created_extra_vertical_cross_section.py dcb

 COMMAND LINE ARGUMENTS
   PARAMETER              - derived model parameter to be created (cb, vpvs, ti, dcb, dvs, etc.)

 DESCRIPTION
   Calculates additional model parameters from the ones created by BKMNS2DD. You can calculate
   the isotropic velocities using the parameter 'iso' (or 'diso', in case you want the Vp and Vs
   perturbations). You will need the files for Vpv, Vph, Vsv, Vsh, and Eta, created at the same
   region and with the same resolution. To compute the bulk sound speed (also needs Vpv, Vph, Vsv,
   Vsh, and Eta), use parameter 'cb' (or 'dcb', in case you want the perturbations). The same
   applies to the Vp/Vs (use parameter 'vpvs'). To calculate the transverse isotropy from Vsv and
   Vsh, use the parameter 'ti'.\n"""

  print (help)


if __name__ == '__main__':

  if len (sys.argv) != 2:

    sys.exit (helpMenu ())

  prm = sys.argv[1]

  if prm == 'iso' or prm == 'cb' or prm == 'vpvs' or prm == 'diso' or prm == 'dcb':

    file1 = 'vpv_VCS.dat'
    file2 = 'vph_VCS.dat'
    file3 = 'vsv_VCS.dat'
    file4 = 'vsh_VCS.dat'
    file5 = 'eta_VCS.dat'

    depth, nrad, ndel, gcline, radius, delta, M1 = readGridFile (file1)
    depth, nrad, ndel, gcline, radius, delta, M2 = readGridFile (file2)
    depth, nrad, ndel, gcline, radius, delta, M3 = readGridFile (file3)
    depth, nrad, ndel, gcline, radius, delta, M4 = readGridFile (file4)
    depth, nrad, ndel, gcline, radius, delta, M5 = readGridFile (file5)

    M1 = np.sqrt (((8 + 4 * M5) * M2 ** 2 + 3 * M1 ** 2 +\
                   (8 - 8 * M5) * M3 ** 2) / 15)

    M2 = np.sqrt (((1 - 2 * M5) * M2 ** 2 + M1 ** 2 + 5 * M4 ** 2 +\
                   (6 + 4 * M5) * M3 ** 2) / 15)

    if prm == 'iso':

      label = ' Vp (km/s)'
      name  = 'vp_VCS.dat'

      writeOutput (name, label, depth, nrad, ndel,
                   gcline, radius, delta, M1)

      label = ' Vs (km/s)'
      name  = 'vs_VCS.dat'

      writeOutput (name, label, depth, nrad, ndel,
                   gcline, radius, delta, M2)

    elif prm == 'cb':

      M = np.sqrt (M1 ** 2 - (4.0 / 3) * M2 ** 2)

      label = ' Cbulk (km/s)'
      name  = 'cb_VCS.dat'

      writeOutput (name, label, depth, nrad, ndel,
                   gcline, radius, delta, M)

    elif prm == 'vpvs':

      M = M1 / M2

      label = '    Vp/Vs'
      name  = 'vpvs_VCS.dat'

      writeOutput (name, label, depth, nrad, ndel,
                   gcline, radius, delta, M)

    else:

      mm1 = meanModel ('1D_mean/vpv.dat')
      mm2 = meanModel ('1D_mean/vph.dat')
      mm3 = meanModel ('1D_mean/vsv.dat')
      mm4 = meanModel ('1D_mean/vsh.dat')
      mm5 = meanModel ('1D_mean/eta.dat')

      mvp = np.sqrt (((8 + 4 * mm5 (depth)) * mm2 (depth) ** 2 +\
                      3 * mm1 (depth) ** 2 +\
                      (8 - 8 * mm5 (depth)) * mm3 (depth) ** 2) / 15)

      mvs = np.sqrt (((1 - 2 * mm5 (depth)) * mm2 (depth) ** 2 +\
                      mm1 (depth) ** 2 + 5 * mm4 (depth) ** 2 +\
                      (6 + 4 * mm5 (depth)) * mm3 (depth) ** 2) / 15)

      if prm == 'diso':

        M = 100 * np.log (M1 / mvp)

        label = ' dlnVp (%)'
        name  = 'dvp_VCS.dat'

        writeOutput (name, label, depth, nrad, ndel,
                   gcline, radius, delta, M)

        M = 100 * np.log (M2 / mvs)

        label = ' dlnVs (%)'
        name  = 'dvs_VCS.dat'

        writeOutput (name, label, depth, nrad, ndel,
                     gcline, radius, delta, M)

      else:

        M   = np.sqrt (M1 ** 2 - (4.0 / 3) * M2 ** 2)
        mmv = np.sqrt (mvp ** 2 - (4.0 / 3) * mvs ** 2)

        M = 100 * np.log (M / mmv)

        label = ' dlnCbulk (%)'
        name  = 'dcb_VCS.dat'

        writeOutput (name, label, depth, nrad, ndel,
                     gcline, radius, delta, M)

  elif prm == 'ti':

    file1 = 'vsv_VCS.dat'
    file2 = 'vsh_VCS.dat'

    depth, nrad, ndel, gcline, radius, delta, M1 = readGridFile (file1)
    depth, nrad, ndel, gcline, radius, delta, M2 = readGridFile (file2)

    mm1 = meanModel ('1D_mean/vsv.dat')
    mm2 = meanModel ('1D_mean/vsh.dat')

    M = 100 * (np.log (M1 / M2) - np.log (mm1 (depth) / mm2 (depth)))

    label = '   Ti (%)'
    name  = 'ti_VCS.dat'

    writeOutput (name, label, depth, nrad, ndel,
                 gcline, radius, delta, M)

  else:

    sys.exit ('Error: {} parameter is not available!'.format (prm))

