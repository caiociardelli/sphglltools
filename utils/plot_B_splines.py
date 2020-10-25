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

 PLOT_B_SPLINES

 USAGE
   ./utils/plot_B_splines.py PARAMETER

 EXAMPLE
   ./utils/plot_B_splines.py vs

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be plotted (vp, vs, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the recreation of the mean model from the B-splines against the mean model computed using
   the GLL2MEAN to check the quality of the radial expansion.

----------------------------------------------------------------------------------------------- */
"""

from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from scipy.integrate import simps

import sys
import numpy as np
import matplotlib.pyplot as plt

EARTH_R = 6371.0


def depth2R (depth):
  """
  Converts depth to normalized Earth radius
  """
  return (EARTH_R - depth) / EARTH_R


def readKnotsAndCoffs (name):
  """
  Reads knots and coefficients
  """
  with open (name, 'r') as FILE:

    lines = FILE.readlines ()[1:]

    i = 0; t = np.array ([])

    while lines[i][:-1] != 'Coefficients':

      knot = float (lines[i][:-1])

      t = np.append (t, np.array (knot))

      i += 1

    c = np.array ([])

    for line in lines[i + 1:]:

      coeff = float (line[:-1])

      c = np.append (c, np.array (coeff))

  return t, c


def readMeanModel (name):
  """
  Reads the mean model
  """
  Mm = {'r' : np.array ([]),
        'a' : np.array ([]),
        'g' : np.array ([])}

  with open (name, 'r') as FILE:

    for line in FILE:

      if line[0] == '#': continue

      data = line.split ()

      r = depth2R (float (data[0]))
      a = float (data[1])
      g = float (data[2])

      Mm['r'] = np.append (Mm['r'], np.array (r))
      Mm['a'] = np.append (Mm['a'], np.array (a))
      Mm['g'] = np.append (Mm['g'], np.array (g))

  return Mm


def aMeanModel (r, Mm):
  """
  Computes the arithmetic mean model
  """
  N = Mm['r'].size

  i = int (N * (r - Mm['r'][0]) / (Mm['r'][N - 1] - Mm['r'][0]) + 0.5)

  if i <= 0: return Mm['a'][0]
  elif i > N - 1: return Mm['a'][N - 1]

  m = (Mm['a'][i] - Mm['a'][i - 1]) / (Mm['r'][i] - Mm['r'][i - 1])

  return Mm['a'][i] + m * (r - Mm['r'][i]);


def aMM (r, Mm):
  """
  Evaluates the mean model for all r
  """
  N = r.size

  y = np.zeros (N)

  for i in range (N):

    y[i] = aMeanModel (r[i], Mm)

  return y


def bSpl (k, n, t, x):
  """
  Evaluates B_{k,n}(x)
  """
  if n == 0:

    return 1.0 if t[k] <= x and x < t[k + 1] else 0.0

  else:

    a = t[k + n] - t[k]; b = t[k + n + 1] - t[k + 1]

    c1 = (x - t[k]) / a if a > 0 else 0
    c2 = (t[k + n + 1] - x) / b if b > 0 else 0

    return c1 * bSpl(k, n - 1, t, x) + c2 * bSpl (k + 1, n - 1, t, x)


def bSpline (k, n, t, x):
  """
  Evaluates B_{k,n}
  """
  y = np.empty (x.size)

  for i in range (x.size):

    y[i] = bSpl (k, n, t, x[i])

  return y


def bSplines (n, t, x):
  """
  Evaluates B_{n}
  """
  Nk = t.size - n - 1

  B = [bSpline (k, n, t, x) for k in range (Nk)]; B[-1][-1] = 1.0

  return B


def sB (n, c, B):
  """
  Evaluates S(x)
  """
  s = np.zeros (B[0].size)

  for k in range (len (B)):

    s += c[k] * B[k]

  return s


def helpMenu ():
  """
  Prints help menu
  """
  help = """\n Error: wrong number of parameters on the comand line...

 PLOT_B_SPLINES

 USAGE
   ./utils/plot_B_splines.py PARAMETER

 EXAMPLE
   ./utils/plot_B_splines.py vs

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be plotted (vp, vs, rho, eta, vsv, etc.)

 DESCRIPTION
   Plots the recreation of the mean model from the B-splines against the mean model computed using
   the GLL2MEAN to check the quality of the radial expansion.\n"""

  print (help)


if __name__ == '__main__':

  if len (sys.argv) != 2:

    sys.exit (helpMenu ())

  N = 500
  n = 3

  prm = sys.argv[1]

  if prm[0] == 'v':

    label = '{}_{{{}}}'.format (prm[0].capitalize (),
                                prm[1:].upper ())
    unity = r'$$\,\,\left(km\,s^{{-1}}\right)'

  elif prm[0] == 'r':

    label = r'\rho'
    unity = r'$$\,\,\left(g\,{{cm}}^{{3}}\right)'

  elif prm[0] == 'e':

    label = r'\eta'
    unity = ''

  elif prm[0] == 'q':

    label = r'Q_{\mu}'
    unity = ''

  t1, c1 = readKnotsAndCoffs ('{}_Z4_KC.dat'.format (prm))
  t2, c2 = readKnotsAndCoffs ('{}_Z3_KC.dat'.format (prm))
  t3, c3 = readKnotsAndCoffs ('{}_Z2_KC.dat'.format (prm))

  t = np.append (np.append (t1, t2), t3)
  c = np.append (np.append (c1, c2), c3)

  Mm = readMeanModel ('{}.dat'.format (prm))

  r1 = np.linspace (t1[0], t1[-1], N)
  r2 = np.linspace (t2[0], t2[-1], N)
  r3 = np.linspace (t3[0], t3[-1], N)

  r = np.append (np.append (r1, r2), r3)

  B1 = bSplines (n, t1, r1)
  B2 = bSplines (n, t2, r2)
  B3 = bSplines (n, t3, r3)

  y1 = aMM (r1, Mm)
  y2 = aMM (r2, Mm)
  y3 = aMM (r3, Mm)

  y = np.append (np.append (y1, y2), y3)

  s1 = sB (n, c1, B1)
  s2 = sB (n, c2, B2)
  s3 = sB (n, c3, B3)

  s = np.append (np.append (s1, s2), s3)

  fig = plt.figure (figsize = (10, 10), dpi = 200)

  for k in range (len (B1)):

    plt.plot (3 * B1[k], r1, linewidth = 2)

  for k in range (len (B2)):

    plt.plot (3 * B2[k], r2, linewidth = 2)

  for k in range (len (B3)):

    plt.plot (3 * B3[k], r3, linewidth = 2)

  plt.plot (y, r, color = 'red',
            linewidth = 2, label = 'GLADM15')
  plt.plot (s, r, color = 'black', linestyle = 'dashed',
            linewidth = 2, label = r'$S\,(x)$')
  plt.xlabel (r'${}{}$'.format (label, unity), fontsize = 32)
  plt.ylabel (r'$r_n$', labelpad = 20, fontsize = 32, rotation = 0)
  plt.tick_params (axis = 'x', labelsize = 20)
  plt.tick_params (axis = 'y', labelsize = 20)
  plt.xlim (0, 1.1 * y.max ())
  plt.ylim (0.99 * t[0], 1.01 * t[-1])
  plt.legend (shadow = True, framealpha = 1.0,
              handlelength = 1.8, fontsize = 24)
  plt.grid (linestyle = 'dashed', linewidth = 0.1)
  plt.title ('Radial cubic B-splines', fontsize = 34, y = 1.03)

  plt.subplots_adjust (left = 0.10, right = 0.99, bottom = 0.09, top = 0.93)

  plt.savefig ('B-splines_{}.pdf'.format (prm))

