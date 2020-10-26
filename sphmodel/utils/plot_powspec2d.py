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

 PLOT_POWSPEC2D

 USAGE
   ./utils/plot_powspec2d.py PARAMETER

 EXAMPLE
   ./utils/plot_powspec2d.py vsv

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter for which the power spectrum will be plotted
                            (vph, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the 2D power spectrum created by POWSPEC2D.

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


def readPowSpec2D (name):

  powspec = list ()

  with open (name, 'r') as FILE:

    line = FILE.readline ()
    data = line.split ()

    N    = int (data[2])
    nmax = int (data[3])

    next (FILE)

    for line in FILE:

      powspec += [line.split ()[1:]]

  Ps = np.zeros ((nmax + 1, N))

  for n in range (nmax + 1):

    for depth in range (N):

      Ps[n][N - depth - 1] = float (powspec[depth][n])

  Ps[0] = min (Ps[1:].ravel ())
  Ps   /= max (Ps.ravel ())

  return N, nmax, Ps.T


def helpMenu ():

  help = """\n Error: wrong number of parameters on the comand line...

 PLOT_POWSPEC2D

 USAGE
   ./utils/plot_powspec2d.py PARAMETER

 EXAMPLE
   ./utils/plot_powspec2d.py vsv

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter for which the power spectrum will be plotted
                            (vp, vs, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the 2D power spectrum created by POWSPEC2D.\n"""

  print (help)


if __name__ == '__main__':

  if len (sys.argv) != 2:

    sys.exit (helpMenu ())

  prm = sys.argv[1]

  N, nmax, Ps = readPowSpec2D ('{}_pwspc2D.dat'.format (prm))

  Ps = np.log (np.exp (1) * Ps)

  fig = plt.figure (figsize = (7, 9), dpi = 200)

  plt.imshow (Ps, origin = 'lower',
              cmap = 'Spectral_r',
              aspect = 'auto',
              interpolation = 'none',
              extent = [0, nmax + 1, N + 80, 80])

  plt.axhline (y = 410, color = 'black',
               linestyle = 'dashed',
               linewidth = 2, zorder = 1)
  plt.axhline (y = 650, color = 'black',
               linestyle = 'dashed',
               linewidth = 2, zorder = 1)

  plt.text (13, 350, '410',
            fontsize = 16,
            ha = 'center', va = 'center',
            color = 'black', zorder = 4)
  plt.text (13, 590, '650', fontsize = 16,
            ha = 'center', va = 'center',
            color = 'black', zorder = 4)

  plt.tick_params (axis = 'x', labelsize = 16)
  plt.tick_params (axis = 'y', labelsize = 16)
  plt.xlabel ('degree', fontsize = 18)
  plt.ylabel ('depth (km)', fontsize = 18)
  plt.xlim (1, nmax + 1)
  plt.ylim (N + 80, 80)
  plt.title ('2D Power Spectrum', fontsize = 30, y = 1.02)

  cax = fig.add_axes ([0.15, 0.08, 0.79, 0.04])
  cbar = plt.colorbar (cax = cax, orientation = 'horizontal')
  cax.set_xlabel ('log power', rotation = 0, fontsize = 18)
  cax.tick_params (labelsize = 16)

  plt.subplots_adjust (left = 0.18, bottom = 0.22, top = 0.92, right = 0.92)

  plt.savefig ('powspec_{}_2d.pdf'.format (prm))

