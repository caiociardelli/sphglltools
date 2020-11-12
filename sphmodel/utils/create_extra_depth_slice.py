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

 CREATE_EXTRA_DEPTH_SLICE

 USAGE
   ./utils/created_extra_depth_slice.py PARAMETER DEPTH

 EXAMPLE
   ./utils/created_extra_depth_slice.py cb 100

 COMMAND LINE ARGUMENTS
   PARAMETER              - derived model parameter to be created (cb, vpvs, ti, dcb, dvs, etc.)
   DEPTH                  - depth in which you created the depth slice

 DESCRIPTION
Calculates additional model parameters from the ones created by BKMNS2LL. You can calculate
   the isotropic velocities using the parameter 'iso' (or 'diso', in case you want the Vp and Vs
   perturbations). You will need the files for Vpv, Vph, Vsv, Vsh, and Eta, created at the
   corresponding depth and with the same resolution. To compute the bulk sound speed (also needs
   Vpv, Vph, Vsv, Vsh, and Eta), use parameter 'cb' (or 'dcb', in case you want the perturbations).
   The same applies to the Vp/Vs (use parameter 'vpvs'). To calculate the transverse isotropy
   from Vsv and Vsh, use the parameter 'ti'.

----------------------------------------------------------------------------------------------- */
"""

from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors


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

    depth = float (line[4])
    nlat  = int (line[5])
    nlon  = int (line[6])

    phi   = np.zeros (nlat * nlon)
    theta = np.zeros (nlat * nlon)
    M     = np.zeros (nlat * nlon)

    next (FILE)

    i = 0

    for line in FILE:

      data = line.split ()

      theta[i] = float (data[0])
      phi[i]   = float (data[1])
      M[i]     = float (data[2])

      i += 1

  return depth, nlat, nlon, theta, phi, M


def writeOutput (name, label, depth, nlat, nlon, theta, phi, M):
  """
  Writes grid file for the derived parameter
  """
  with open (name, 'w') as FILE:

    header  = '#depth (km) nlat nlon: {:g} {} {}\n'.format (depth, nlat, nlon)
    header += '#latitude (degrees)   longitude (degrees)   {}\n'.format (label)

    FILE.write (header.encode ().decode ('utf-8'))

    for t, p, m in zip (theta, phi, M):

      FILE.write ('{:15.7f} {:21.7f} {:18E}\n'.format (t, p, m).encode ().decode ('utf-8'))


def helpMenu ():

  help = """\n Error: wrong number of parameters on the comand line...

 CREATE_EXTRA_DEPTH_SLICE

 USAGE
   ./utils/created_extra_depth_slice.py PARAMETER DEPTH

 EXAMPLE
   ./utils/created_extra_depth_slice.py cb 100

 COMMAND LINE ARGUMENTS
   PARAMETER              - derived model parameter to be created (cb, vpvs, ti, dcb, dvs, etc.)
   DEPTH                  - depth in which you created the depth slice

 DESCRIPTION
   Calculates additional model parameters from the ones created by BKMNS2LL. You can calculate
   the isotropic velocities using the parameter 'iso' (or 'diso', in case you want the Vp and Vs
   perturbations). You will need the files for Vpv, Vph, Vsv, Vsh, and Eta, created at the
   corresponding depth and with the same resolution. To compute the bulk sound speed (also needs
   Vpv, Vph, Vsv, Vsh, and Eta), use parameter 'cb' (or 'dcb', in case you want the perturbations).
   The same applies to the Vp/Vs (use parameter 'vpvs'). To calculate the transverse isotropy
   from Vsv and Vsh, use the parameter 'ti'.\n"""

  print (help)


if __name__ == '__main__':

  if len (sys.argv) != 3:

    sys.exit (helpMenu ())

  prm   = sys.argv[1]
  depth = sys.argv[2]

  if prm == 'iso' or prm == 'cb' or prm == 'vpvs' or prm == 'diso' or prm == 'dcb':

    file1 = 'vpv_{}_DS.dat'.format (depth)
    file2 = 'vph_{}_DS.dat'.format (depth)
    file3 = 'vsv_{}_DS.dat'.format (depth)
    file4 = 'vsh_{}_DS.dat'.format (depth)
    file5 = 'eta_{}_DS.dat'.format (depth)

    depth, nlat, nlon, theta, phi, M1 = readGridFile (file1)
    depth, nlat, nlon, theta, phi, M2 = readGridFile (file2)
    depth, nlat, nlon, theta, phi, M3 = readGridFile (file3)
    depth, nlat, nlon, theta, phi, M4 = readGridFile (file4)
    depth, nlat, nlon, theta, phi, M5 = readGridFile (file5)

    M1 = np.sqrt (((8 + 4 * M5) * M2 ** 2 + 3 * M1 ** 2 +\
                   (8 - 8 * M5) * M3 ** 2) / 15)

    M2 = np.sqrt (((1 - 2 * M5) * M2 ** 2 + M1 ** 2 + 5 * M4 ** 2 +\
                   (6 + 4 * M5) * M3 ** 2) / 15)

    if prm == 'iso':

      label = ' Vp (km/s)'
      name  = 'vp_{:g}_DS.dat'.format (depth)

      writeOutput (name, label, depth,
                   nlat, nlon, theta, phi, M1)

      label = ' Vs (km/s)'
      name  = 'vs_{:g}_DS.dat'.format (depth)

      writeOutput (name, label, depth,
                   nlat, nlon, theta, phi, M2)

    elif prm == 'cb':

      M = np.sqrt (M1 ** 2 - (4.0 / 3) * M2 ** 2)

      label = ' Cbulk (km/s)'
      name  = 'cb_{:g}_DS.dat'.format (depth)

      writeOutput (name, label, depth,
                   nlat, nlon, theta, phi, M)

    elif prm == 'vpvs':

      M = M1 / M2

      label = '    Vp/Vs'
      name  = 'vpvs_{:g}_DS.dat'.format (depth)

      writeOutput (name, label, depth,
                   nlat, nlon, theta, phi, M)

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
        name  = 'dvp_{:g}_DS.dat'.format (depth)

        writeOutput (name, label, depth,
                     nlat, nlon, theta, phi, M)

        M = 100 * np.log (M2 / mvs)

        label = ' dlnVs (%)'
        name  = 'dvs_{:g}_DS.dat'.format (depth)

        writeOutput (name, label, depth,
                     nlat, nlon, theta, phi, M)

      else:

        M   = np.sqrt (M1 ** 2 - (4.0 / 3) * M2 ** 2)
        mmv = np.sqrt (mvp ** 2 - (4.0 / 3) * mvs ** 2)

        M = 100 * np.log (M / mmv)

        label = ' dlnCbulk (%)'
        name  = 'dcb_{:g}_DS.dat'.format (depth)

        writeOutput (name, label, depth,
                     nlat, nlon, theta, phi, M)

  elif prm == 'ti':

    file1 = 'vsv_{}_DS.dat'.format (depth)
    file2 = 'vsh_{}_DS.dat'.format (depth)

    depth, nlat, nlon, theta, phi, M1 = readGridFile (file1)
    depth, nlat, nlon, theta, phi, M2 = readGridFile (file2)

    mm1 = meanModel ('1D_mean/vsv.dat')
    mm2 = meanModel ('1D_mean/vsh.dat')

    M = 100 * (np.log (M1 / M2) - np.log (mm1 (depth) / mm2 (depth)))

    label = '   Ti (%)'
    name  = 'ti_{:g}_DS.dat'.format (depth)

    writeOutput (name, label, depth,
                 nlat, nlon, theta, phi, M)

  else:

    sys.exit ('Error: {} parameter is not available!'.format (prm))

