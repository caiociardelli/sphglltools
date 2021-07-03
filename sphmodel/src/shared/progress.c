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

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <time.h>

void progressBar (unsigned index, unsigned step, unsigned n, clock_t starttime)
{
  /* Creates a progress bar */
  if (index % step == 0 || index == n - 1)
  {
    double percentage = 100 * (double) (index + 1) / n;
    unsigned p = percentage / 2;

    double cpu_time_used = ((double) (clock () - starttime)) / (60 * CLOCKS_PER_SEC);

    fprintf (stderr, "\r [");

    for (unsigned j = 0; j < p; j++)

      fprintf (stderr, "#");

    for (unsigned j = p; j < 50; j++)

      fprintf (stderr, " ");

    fprintf (stderr, "] [%.1lf%%] [Elapsed time: %.2lf min]", percentage, cpu_time_used);
  }
}
