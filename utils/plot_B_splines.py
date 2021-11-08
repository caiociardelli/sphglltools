#!/usr/bin/python3
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

-----------------------------------------------------------------------------------------------
"""

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


def computeError (r, y, s):
  """
  Evaluates the relative error between S(x) and the mean model
  """
  err   = 100 * np.log (s / y)
  err_c = err.copy ()
  dt    = r[1] - r[0]

  i = 1

  while i < len (err) - 1:

    slope1 = abs ((err_c[i] - err_c[i - 1]) / dt)
    slope2 = abs ((err_c[i] - err_c[i + 1]) / dt)

    # At the seismic discontinuities, the ratio is
    # not meaningful due to limited sampling.
    # To get arround this problem, we find them by
    # measuring the slope and interpolate the values.
    if slope1 > 1000 and slope2 > 1000:

      a = (err_c[i + 2] - err_c[i - 1]) / (3 * dt)
      b = err_c[i - 1]

      err[i]     = a * dt + b
      err[i + 1] = a * 2 * dt + b

      i += 2

    i += 1

  return err


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

    label = f'{prm[0].capitalize ()}_{{{prm[1:].upper ()}}}'
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

  t1, c1 = readKnotsAndCoffs (f'{prm}_Z4_KC.dat')
  t2, c2 = readKnotsAndCoffs (f'{prm}_Z3_KC.dat')
  t3, c3 = readKnotsAndCoffs (f'{prm}_Z2_KC.dat')

  t = np.append (np.append (t1, t2), t3)
  c = np.append (np.append (c1, c2), c3)

  Mm = readMeanModel (f'{prm}.dat')

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

  fig, (axis1, axis2) = plt.subplots (1, 2, figsize = (17, 10),
                                      gridspec_kw = {'width_ratios': [3, 1]},
                                      dpi = 200)

  axis1.axhline (y = 0.987443, color = 'gray', linewidth = 2)
  axis1.axhline (y = 0.935646, color = 'gray', linewidth = 2)
  axis1.axhline (y = 0.897975, color = 'gray', linewidth = 2)
  axis1.axhline (y = 0.546225, color = 'gray', linewidth = 2)

  for k in range (len (B1)):

    axis1.plot (3 * B1[k], r1, linewidth = 2)

  for k in range (len (B2)):

    axis1.plot (3 * B2[k], r2, linewidth = 2)

  for k in range (len (B3)):

    axis1.plot (3 * B3[k], r3, linewidth = 2)

  axis1.plot (y, r, color = 'red',
              linewidth = 2, label = 'GLAD-M15')
  axis1.plot (s, r, color = 'black', linestyle = 'dashed',
              linewidth = 2, label = r'$S\,(x)$')

  axis1.set_xlabel (r'${}{}$'.format (label, unity), fontsize = 32)
  axis1.set_ylabel (r'$r_n$', labelpad = 20, fontsize = 32, rotation = 0)
  axis1.tick_params (axis = 'x', labelsize = 20)
  axis1.tick_params (axis = 'y', labelsize = 20)
  axis1.set_xlim (0, 1.1 * y.max ())
  axis1.set_ylim (0.99 * t[0], 1.01 * t[-1])
  axis1.grid (linestyle = 'dashed', linewidth = 0.1)
  axis1.set_title ('Radial cubic B-splines', fontsize = 34, y = 1.03)

  axis1.text (3.0, 0.96, 'Upper Mantle', horizontalalignment = 'center',
              color = 'darkred', fontsize = 18)
  axis1.text (3.0, 0.91, 'Transition Zone', horizontalalignment = 'center',
              color = 'darkred', fontsize = 18)
  axis1.text (3.0, 0.75, 'Lower Mantle', horizontalalignment = 'center',
              color = 'darkred', fontsize = 18)

  axis1.legend (shadow = True, framealpha = 1.0,
              handlelength = 1.8, fontsize = 22, loc = 1)

  err = computeError (r, y, s)

  xmin = -np.ceil (np.max (np.abs (err)))
  xmax =  np.ceil (np.max (np.abs (err)))

  xticks = np.linspace (xmin, xmax, 5)

  axis2.axhline (y = 0.987443, color = 'gray', linewidth = 2)
  axis2.axhline (y = 0.935646, color = 'gray', linewidth = 2)
  axis2.axhline (y = 0.897975, color = 'gray', linewidth = 2)
  axis2.axhline (y = 0.546225, color = 'gray', linewidth = 2)

  axis2.plot (np.zeros (len (r)), r, color = 'red',
              linewidth = 2, label = 'GLAD-M15')
  axis2.plot (err, r, color = 'black', linestyle = 'dashed',
              linewidth = 2, label = r'$S\,(x)$')
  axis2.set_xlabel (r'Error (%)', fontsize = 32)
  axis2.tick_params (axis = 'x', labelsize = 20)
  axis2.set_xticks (xticks)

  for tick in axis2.yaxis.get_major_ticks ():

    tick.tick1line.set_visible (False)
    tick.tick2line.set_visible (False)
    tick.label1.set_visible (False)
    tick.label2.set_visible (False)

  axis2.set_xlim (xmin, xmax)
  axis2.set_ylim (0.99 * t[0], 1.01 * t[-1])
  axis2.grid (linestyle = 'dashed', linewidth = 0.1)
  axis2.set_title ('Relative Error', fontsize = 34, y = 1.03)

  plt.subplots_adjust (left = 0.062, right = 0.985, bottom = 0.092, top = 0.925)

  plt.savefig (f'B-splines_{prm}.pdf')

