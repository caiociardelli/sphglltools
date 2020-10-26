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

 PLOT_POWSPEC

 USAGE
   ./utils/plot_powspec.py PARAMETER

 EXAMPLE
   ./utils/plot_powspec.py vsv

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter for which the power spectrum will be plotted
                            (vph, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the 1D power spectrum created by POWSPEC.

----------------------------------------------------------------------------------------------- */
"""

from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def readPowSpec1D (name):

  with open (name, 'r') as FILE:

    line = FILE.readline ()
    data = line.split ()

    nmax = int (data[1])

    Ps = np.empty (nmax + 1)

    next (FILE)

    index = 0

    for line in FILE:

      data = line.split ()
      Ps[index] = float (data[1])

      index += 1

  Ps[0] = min (Ps[1:])
  Ps   /= max (Ps)

  return nmax, Ps


def helpMenu ():

  help = """\n Error: wrong number of parameters on the comand line...

 PLOT_POWSPEC

 USAGE
   ./utils/plot_powspec.py PARAMETER

 EXAMPLE
   ./utils/plot_powspec.py vsv

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter for which the power spectrum will be plotted
                            (vph, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the 1D power spectrum created by POWSPEC.\n"""

  print (help)


if __name__ == '__main__':

  model = 'S362ANI'

  if len (sys.argv) != 2:

    sys.exit (helpMenu ())

  prm = sys.argv[1]

  plt.figure ()

  nmax, Ps = readPowSpec1D ('{}_pwspc.dat'.format (prm))

  x  = np.linspace (1, nmax, nmax)
  Ps = np.log (np.exp (1) * Ps[1:])

  plt.plot (x, Ps, color = 'blue',
            linewidth = 3, label = model)

  plt.tick_params (axis = 'x', labelsize = 16)
  plt.tick_params (axis = 'y', labelsize = 16)
  plt.ylabel ('log power', fontsize = 18)
  plt.xlabel ('degree', fontsize = 18)
  plt.xlim (0, nmax + 1)
  plt.ylim (Ps.min () - 1, Ps.max () + 1)
  plt.legend ()
  plt.grid ()
  plt.title ('Power Spectrum', fontsize = 30, y = 1.02)

  plt.subplots_adjust (left = 0.12, bottom = 0.12, top = 0.90, right = 0.96)

  plt.savefig ('powspec_{}.pdf'.format (prm))

