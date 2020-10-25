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

-----------------------------------------------------------------------------------------------

 SMOOTH

 USAGE
   mpiexec -n 24 bin/smooth INPUT_TOPO_DIRECTORY INPUT_MODEL_DIRECTORY OUTPUT_DIRECTORY

 EXAMPLE
   mpiexec -n 24 bin/smooth data/INPUT/ data/INPUT/ smooth/

 COMMAND LINE ARGUMENTS
   INPUT_TOPO_DIRECTORY   - directory containing the input topological files
   INPUT_TOPO_DIRECTORY   - directory containing the input model files
   OUTPUT_DIRECTORY       - directory where the routine will write the output files

 DESCRIPTION
   Reads the location of the topological files, the location of the model files, and the output
   directory name from the command line and smooths the desired parameters
   (defined in 'config.h').

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "metrics_smooth.h"
#include "coordinates.h"
#include "structs_smooth.h"
#include "boundaries_smooth.h"
#include "profile.h"
#include "io_smooth.h"
#include "progress.h"
#include "constants.h"
#include "config_smooth.h"

static double arcDistance (struct Point *p1, struct Point *p2)
{
  double cos_D = (p1->x * p2->x + p1->y * p2->y + p1->z * p2->z) / (p1->r * p2->r);

  if (cos_D >  1) return 0;
  if (cos_D < -1) return PI * p2->r;

  return acos (cos_D) * p2->r;
}

static unsigned r2Index (double r, double dr_i)
{
  int index = ((r - CMB_R) * dr_i + 0.5);
  int n = (int) NP;

  if (index < 0)  return 0;
  if (index >= n) return n - 1;

  return (unsigned) index;
}

static void smooth3D (int ic, double rmax, struct Smooth sm[NP],
                      unsigned nel, unsigned li[nel], unsigned lo[NEL],
                      struct Point Ti[nel][NX][NY][NZ],
                      struct Point To[NEL][NX][NY][NZ],
                      struct Parameters Mi[nel][NX][NY][NZ],
                      struct Parameters Mo[NEL][NX][NY][NZ])
{
  unsigned hx = NX / 2, hy = NY / 2, hz = NZ / 2;

  double dr_i = (NP - 1) / (MAX_SURFACE_R - CMB_R);

  struct Radii imr[nel], omr[NEL];

  computeElRadii (nel, Ti, To, imr, omr);

  bool Far[nel][NEL];

  for (unsigned elo = 0; elo < NEL; elo++)
  {
    struct Point *oc = &To[elo][hx][hy][hz];

    double search_h = 0.0;
    double search_v = 0.0;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          unsigned index = r2Index (To[elo][i][j][k].r, dr_i);

          double sm_search_h = sm[index].search_h;
          double sm_search_v = sm[index].search_v;

          if (sm_search_h > search_h) search_h = sm_search_h;
          if (sm_search_v > search_v) search_v = sm_search_v;
        }

    for (unsigned eli = 0; eli < nel; eli++)
    {
      Far[eli][elo] = true;

      if ((SEPARATE_ZONES &&
           li[eli] != lo[elo]) || oc->r > rmax) continue;

      struct Point *ic = &Ti[eli][hx][hy][hz];

      double rv = imr[eli].rv + omr[elo].rv;
      double rd = radialDistance (ic, oc);

      if (rd > search_v + rv) continue;

      double r  = imr[eli].r + omr[elo].r;
      double ad = arcDistance (ic, oc);

      if (ad > search_h + r) continue;

      double c = square (ad / (search_h + r)) +
                 square (rd / (search_v + rv));

      if (c <= 1) Far[eli][elo] = false;
    }
  }

  clock_t starttime = clock ();

  for (unsigned elo = 0; elo < NEL; elo++)
  {
    struct Point oc = To[elo][hx][hy][hz];

    if (oc.r <= rmax)
    {
      for (unsigned io = 0; io < NX; io++)

        for (unsigned jo = 0; jo < NY; jo++)

          for (unsigned ko = 0; ko < NZ; ko++)
          {
            #if VP
            Mo[elo][io][jo][ko].vp  = 0.0;
            #endif

            #if VS
            Mo[elo][io][jo][ko].vs  = 0.0;
            #endif

            #if RHO
            Mo[elo][io][jo][ko].rho = 0.0;
            #endif

            #if VPV
            Mo[elo][io][jo][ko].vpv = 0.0;
            #endif

            #if VPH
            Mo[elo][io][jo][ko].vph = 0.0;
            #endif

            #if VSV
            Mo[elo][io][jo][ko].vsv = 0.0;
            #endif

            #if VSH
            Mo[elo][io][jo][ko].vsh = 0.0;
            #endif

            #if ETA
            Mo[elo][io][jo][ko].eta = 0.0;
            #endif

            #if QMU
            Mo[elo][io][jo][ko].qmu = 0.0;
            #endif

            #if APK
            Mo[elo][io][jo][ko].apk = 0.0;
            #endif

            #if BTK
            Mo[elo][io][jo][ko].btk = 0.0;
            #endif

            #if RHK
            Mo[elo][io][jo][ko].rhk = 0.0;
            #endif

            #if BCK
            Mo[elo][io][jo][ko].bck = 0.0;
            #endif

            #if BBK
            Mo[elo][io][jo][ko].bbk = 0.0;
            #endif

            #if BVK
            Mo[elo][io][jo][ko].bvk = 0.0;
            #endif

            #if BHK
            Mo[elo][io][jo][ko].bhk = 0.0;
            #endif

            #if ETK
            Mo[elo][io][jo][ko].etk = 0.0;
            #endif

            #if HSK
            Mo[elo][io][jo][ko].hsk = 0.0;
            #endif

            double sw = 0.0;

            struct Point *op = &To[elo][io][jo][ko];

            unsigned index = r2Index (op->r, dr_i);

            double search_h = sm[index].search_h;
            double search_v = sm[index].search_v;
            double sg_h     = sm[index].sg_h;
            double sg_v     = sm[index].sg_v;
            double sc_h     = sm[index].sc_h;
            double sc_v     = sm[index].sc_v;

            for (unsigned eli = 0; eli < nel; eli++)
            {
              if (Far[eli][elo]) continue;

              for (unsigned ii = 0; ii < NX; ii++)

                for (unsigned ji = 0; ji < NY; ji++)

                  for (unsigned ki = 0; ki < NZ; ki++)
                  {
                    struct Point *ip = &Ti[eli][ii][ji][ki];

                    double rd = radialDistance (ip, op);

                    if (rd > search_v) continue;

                    double ad = arcDistance (ip, op);

                    if (ad > search_h) continue;

                    double ad_2 = square (ad);
                    double rd_2 = square (rd);

                    if (sc_h * ad_2 + sc_v * rd_2 <= 1)
                    {
                      double w = exp (-(sg_h * ad_2 + sg_v * rd_2));

                      #if VP
                      Mo[elo][io][jo][ko].vp  += w * Mi[eli][ii][ji][ki].vp;
                      #endif

                      #if VS
                      Mo[elo][io][jo][ko].vs  += w * Mi[eli][ii][ji][ki].vs;
                      #endif

                      #if RHO
                      Mo[elo][io][jo][ko].rho += w * Mi[eli][ii][ji][ki].rho;
                      #endif

                      #if VPV
                      Mo[elo][io][jo][ko].vpv += w * Mi[eli][ii][ji][ki].vpv;
                      #endif

                      #if VPH
                      Mo[elo][io][jo][ko].vph += w * Mi[eli][ii][ji][ki].vph;
                      #endif

                      #if VSV
                      Mo[elo][io][jo][ko].vsv += w * Mi[eli][ii][ji][ki].vsv;
                      #endif

                      #if VSH
                      Mo[elo][io][jo][ko].vsh += w * Mi[eli][ii][ji][ki].vsh;
                      #endif

                      #if ETA
                      Mo[elo][io][jo][ko].eta += w * Mi[eli][ii][ji][ki].eta;
                      #endif

                      #if QMU
                      Mo[elo][io][jo][ko].qmu += w * Mi[eli][ii][ji][ki].qmu;
                      #endif

                      #if APK
                      Mo[elo][io][jo][ko].apk += w * Mi[eli][ii][ji][ki].apk;
                      #endif

                      #if BTK
                      Mo[elo][io][jo][ko].btk += w * Mi[eli][ii][ji][ki].btk;
                      #endif

                      #if RHK
                      Mo[elo][io][jo][ko].rhk += w * Mi[eli][ii][ji][ki].rhk;
                      #endif

                      #if BCK
                      Mo[elo][io][jo][ko].bck += w * Mi[eli][ii][ji][ki].bck;
                      #endif

                      #if BBK
                      Mo[elo][io][jo][ko].bbk += w * Mi[eli][ii][ji][ki].bbk;
                      #endif

                      #if BVK
                      Mo[elo][io][jo][ko].bvk += w * Mi[eli][ii][ji][ki].bvk;
                      #endif

                      #if BHK
                      Mo[elo][io][jo][ko].bhk += w * Mi[eli][ii][ji][ki].bhk;
                      #endif

                      #if ETK
                      Mo[elo][io][jo][ko].etk += w * Mi[eli][ii][ji][ki].etk;
                      #endif

                      #if HSK
                      Mo[elo][io][jo][ko].hsk += w * Mi[eli][ii][ji][ki].hsk;
                      #endif

                      sw += w;
                    }
                  }
            }

            double isw = 1.0 / sw;

            #if VP
            Mo[elo][io][jo][ko].vp  *= isw;
            #endif

            #if VS
            Mo[elo][io][jo][ko].vs  *= isw;
            #endif

            #if RHO
            Mo[elo][io][jo][ko].rho *= isw;
            #endif

            #if VPV
            Mo[elo][io][jo][ko].vpv *= isw;
            #endif

            #if VPH
            Mo[elo][io][jo][ko].vph *= isw;
            #endif

            #if VSV
            Mo[elo][io][jo][ko].vsv *= isw;
            #endif

            #if VSH
            Mo[elo][io][jo][ko].vsh *= isw;
            #endif

            #if ETA
            Mo[elo][io][jo][ko].eta *= isw;
            #endif

            #if QMU
            Mo[elo][io][jo][ko].qmu *= isw;
            #endif

            #if APK
            Mo[elo][io][jo][ko].apk *= isw;
            #endif

            #if BTK
            Mo[elo][io][jo][ko].btk *= isw;
            #endif

            #if RHK
            Mo[elo][io][jo][ko].rhk *= isw;
            #endif

            #if BCK
            Mo[elo][io][jo][ko].bck *= isw;
            #endif

            #if BBK
            Mo[elo][io][jo][ko].bbk *= isw;
            #endif

            #if BVK
            Mo[elo][io][jo][ko].bvk *= isw;
            #endif

            #if BHK
            Mo[elo][io][jo][ko].bhk *= isw;
            #endif

            #if ETK
            Mo[elo][io][jo][ko].etk *= isw;
            #endif

            #if HSK
            Mo[elo][io][jo][ko].hsk *= isw;
            #endif
          }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (ic == 0) progressBar (elo, 5, NEL, starttime);
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n SMOOTH"

                    "\n\n USAGE"
                    "\n    mpiexec -n 24 bin/smooth INPUT_TOPO_DIRECTORY INPUT_MODEL_DIRECTORY OUTPUT_DIRECTORY"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 24 bin/smooth data/INPUT/ data/INPUT/ smooth/"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    INPUT_TOPO_DIRECTORY   - directory containing the input topological files"
                    "\n    INPUT_TOPO_DIRECTORY   - directory containing the input model files"
                    "\n    OUTPUT_DIRECTORY       - directory where the routine will write the output files"

                    "\n\n DESCRIPTION"
                    "\n    Reads the location of the topological files, the location of the model files, and the output"
                    "\n    directory name from the command line and smooths the desired parameters"
                    "\n    (defined in 'config.h').\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char **argv)
{
  int ic;
  int nc;

  MPI_Init (NULL, NULL);

  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  if (ic == 0 && argc != 4)
  {
    fprintf (stderr, "Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  struct Boundaries gb;

  struct Point To[NEL][NX][NY][NZ];
  struct Parameters Mo[NEL][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nReading output mesh...");

  if (checkTopoAndModelIO (readOutputMeshAndModel (ic, argv[1], argv[2],
                                                   &gb, To, Mo)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  unsigned n;

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nReading smoothing profile...");

  if (checkProfileIO (readNumberOfPoints (&n))) MPI_Abort (MPI_COMM_WORLD, 1);

  double sigma; struct Profile p[n];

  if (checkProfileIO (readSmoothingProfile (n, &sigma, p)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  struct Smooth sm[NP];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nInterpolating profile...");

  interpolate (n, p, sm);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nSmoothing profile...");

  smoothPf (sigma, sm);

  unsigned nel = 0;

  struct ElNode *lle;
  initializeElNode (&lle);

  struct PmNode *llm;
  initializePmNode (&llm);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nScanning input mesh and model files...\n\n");

  if (checkTopoAndModelIO (scanMeshAndModel (ic, argv[1], argv[2],
                                             n, p, &gb, lle, llm, &nel)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  unsigned nels; MPI_Reduce(&nel, &nels, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\nAn average of %u espectral elements were selected "
                                "in each slice.\n\n", nels / nc);

  struct Point Ti[nel][NX][NY][NZ];
  struct Parameters Mi[nel][NX][NY][NZ];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Converting mesh and model linked lists to arrays...\n");

  toArrayElAndPmNodes (nel, lle, llm, Ti, Mi);

  unsigned li[nel], lo[NEL];

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Identifying zones...\n");

  identifyZones (nel, li, Ti);
  identifyZones (NEL, lo, To);

  if (ic == 0 && SEPARATE_ZONES) fprintf (stderr, "Preserving discontinuities...\n");

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Smoothing mesh parameters...\n\n");

  smooth3D (ic, p[0].r, sm, nel, li, lo, Ti, To, Mi, Mo);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "\n\nWriting smoothed model...\n");

  if (checkModelIO (writeModel (ic, argv[3], Mo))) MPI_Abort (MPI_COMM_WORLD, 1);

  MPI_Barrier (MPI_COMM_WORLD);

  if (ic == 0) fprintf (stderr, "Done!\n");

  MPI_Finalize ();

  return 0;
}

