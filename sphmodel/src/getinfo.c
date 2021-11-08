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

-----------------------------------------------------------------------------------------------

 GETINFO

 USAGE
   ./bin/getinfo FILENAME

 EXAMPLE
   ./bin/getinfo vp_VCS.dat

 COMMAND-LINE ARGUMENTS
   FILENAME              - name of the grid file

 DESCRIPTION
   Reads a grid file and outputs the minimum and the maximum values, the mean, and the standard
   deviation. For all the computations, this routine ignores zeroes in the grid. The reason for
   that is because the main propose of this routine is providing information to create the color
   bars in the plotting routines. Zeroes almost always represent null values, corresponding to
   regions outside the mesh and, hence, are not meaningful.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "exmath.h"
#include "constants.h"

static void helpMenu (void)
{
  char *help_menu = "\n GETINFO"

                    "\n\n USAGE"
                    "\n     ./bin/getinfo FILENAME"

                    "\n\n EXAMPLE"
                    "\n     ./bin/getinfo vp_VCS.dat"

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    FILENAME              - name of the grid file"

                    "\n\n DESCRIPTION"
                    "\n    Reads a grid file and outputs the minimum and the maximum values, the mean, and the standard"
                    "\n    deviation. For all the computations, this routine ignores zeroes in the grid. The reason for"
                    "\n    that is because the main propose of this routine is providing information to create the color"
                    "\n    bars in the plotting routines. Zeroes almost always represent null values, corresponding to"
                    "\n    regions outside the mesh and, hence, are not meaningful.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 2)
  {
    fprintf (stderr, "Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  FILE *file = fopen (argv[1], "r");

  if (file == NULL)
  {
    fprintf (stderr, "\n Error: could not read file!\n");

    return 1;
  }

  char line[MAX_STRING_LEN];

  double min =  INFINITY;
  double max = -INFINITY;

  unsigned n = 0;
  double s = 0.0, s2 = 0.0;

  while (fgets (line, MAX_STRING_LEN, file) != NULL)
  {
    if (line[0] != '#')
    {
      char *last = strrchr (line, ' ');
      double v = atof (last + 1);

      if (v)
      {
        if (v > max) max = v;
        if (v < min) min = v;

        s += v; s2 += square (v); n++;
      }
    }
  }

  fclose (file);

  double mean = s / n;
  double sdv  = sqrt ((s2 - square (s) / n) / (n - 1));

  fprintf (stdout, "min/max/mean/sdv -> %lf %lf %lf %lf\n",
           min, max, mean, sdv);

  return 0;
}

