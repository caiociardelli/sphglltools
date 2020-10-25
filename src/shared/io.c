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

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "exmath.h"
#include "legendre.h"
#include "structs.h"
#include "coordinates.h"
#include "io.h"

void initializeElNode (struct ElNode **lle)
{
  /* Initializes a linked list to store the
     spectral elements */
  *lle = malloc (sizeof (struct ElNode));
  (*lle)->next = NULL;
}

void initializePmNode (struct PmNode **llm)
{
  /* Initializes a linked list to store the
     model parameters */
  *llm = malloc (sizeof (struct PmNode));
  (*llm)->next = NULL;
}

static void insertElNode (struct Point T[NX][NY][NZ], struct ElNode *lle,
                          unsigned *nel)
{
  /* Adds a new cell to the linked list of spectral elements */
  struct ElNode *p = lle, *new = malloc (sizeof (struct ElNode));

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)
      {
        new->T[i][j][k].x = T[i][j][k].x;
        new->T[i][j][k].y = T[i][j][k].y;
        new->T[i][j][k].z = T[i][j][k].z;

        new->T[i][j][k].r     = T[i][j][k].r;
        new->T[i][j][k].theta = T[i][j][k].theta;
        new->T[i][j][k].phi   = T[i][j][k].phi;
      }

  new->next = p->next;
  p->next = new;

  *nel += 1;
}

static void insertPmNode (struct Parameters M[NX][NY][NZ],
                          struct PmNode *llm)
{
  /* Adds a new cell to the linked list of model parameters */
  struct PmNode *p = llm, *new = malloc (sizeof (struct PmNode));

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)
      {
        #if VP
        new->M[i][j][k].vp  = M[i][j][k].vp;
        #endif

        #if VS
        new->M[i][j][k].vs  = M[i][j][k].vs;
        #endif

        #if RHO
        new->M[i][j][k].rho = M[i][j][k].rho;
        #endif

        #if VPV
        new->M[i][j][k].vpv = M[i][j][k].vpv;
        #endif

        #if VPH
        new->M[i][j][k].vph = M[i][j][k].vph;
        #endif

        #if VSV
        new->M[i][j][k].vsv = M[i][j][k].vsv;
        #endif

        #if VSH
        new->M[i][j][k].vsh = M[i][j][k].vsh;
        #endif

        #if ETA
        new->M[i][j][k].eta = M[i][j][k].eta;
        #endif

        #if QMU
        new->M[i][j][k].qmu = M[i][j][k].qmu;
        #endif
      }

  new->next = p->next;
  p->next = new;
}

static void deleteElNode (struct ElNode *lle)
{
  /* Deletes a cell of the linked list of spectral elements */
  struct ElNode *q;

  for (struct ElNode *p = lle; p != NULL; p = q)
  {
    q = p->next;
    free (p);
  }

  q = NULL;
}

static void deletePmNode (struct PmNode *llm)
{
  /* Deletes a cell of the linked list of model parameters */
  struct PmNode *q;

  for (struct PmNode *p = llm; p != NULL; p = q)
  {
    q = p->next;
    free (p);
  }

  q = NULL;
}

static unsigned readInputMeshAndModel (int ic, char *prm, unsigned nelm,
                                       struct Boundaries *gb, struct ElNode *lle,
                                       struct PmNode *llm, unsigned *nel)
{
  /* Reads the input mesh and model */
  unsigned N = nelm * NX * NY * NZ;

  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s/proc%06d_reg1_solver_data.bin", prm, ic) < 33) return 1;

  FILE *file = fopen (filename, "rb");

  if (file == NULL) return 2;

  unsigned junk;
  unsigned nspec;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nspec, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned nval;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nval, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xf[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xf, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float yf[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (yf, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float zf[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (zf, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ibool[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ibool, sizeof (unsigned), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned idoubling[nelm];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (idoubling, sizeof (unsigned), nelm, file) != nelm) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ispec_is_tiso[nelm];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ispec_is_tiso, sizeof (unsigned), nelm, file) != nelm) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xix[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xix, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiy[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiy, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiz[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etax[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etay[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etaz[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammax[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammay[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammaz[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);

  #if VP
  if (sprintf (filename, "%s/proc%06d_reg1_vp.bin", prm, ic) < 23) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vp[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VS
  if (sprintf (filename, "%s/proc%06d_reg1_vs.bin", prm, ic) < 23) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vs[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vs, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if RHO
  if (sprintf (filename, "%s/proc%06d_reg1_rho.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float rho[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (rho, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPV
  if (sprintf (filename, "%s/proc%06d_reg1_vpv.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vpv[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vpv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPH
  if (sprintf (filename, "%s/proc%06d_reg1_vph.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vph[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vph, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSV
  if (sprintf (filename, "%s/proc%06d_reg1_vsv.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsv[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSH
  if (sprintf (filename, "%s/proc%06d_reg1_vsh.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsh[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsh, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if ETA
  if (sprintf (filename, "%s/proc%06d_reg1_eta.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float eta[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (eta, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if QMU
  if (sprintf (filename, "%s/proc%06d_reg1_qmu.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float qmu[nelm][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (qmu, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  struct Point T[NX][NY][NZ];
  struct Parameters M[NX][NY][NZ];

  for (unsigned el = 0; el < nelm; el++)
  {
    long double xmin = INFINITY, xmax = -INFINITY;
    long double ymin = INFINITY, ymax = -INFINITY;
    long double zmin = INFINITY, zmax = -INFINITY;
    long double rmin = INFINITY, rmax = -INFINITY;
    long double tmin = INFINITY, tmax = -INFINITY;
    long double pmin = INFINITY, pmax = -INFINITY;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          unsigned ig = ibool[el][k][j][i] - 1;

          long double x = xf[ig];
          long double y = yf[ig];
          long double z = zf[ig];

          T[i][j][k].x = x;
          T[i][j][k].y = y;
          T[i][j][k].z = z;

          long double r, theta, phi;

          xYZ2RThetaPhil (x, y, z, &r, &theta, &phi);

          T[i][j][k].r     = r;
          T[i][j][k].theta = theta;
          T[i][j][k].phi   = phi;

          #if VP
          M[i][j][k].vp  =  vp[el][k][j][i];
          #endif

          #if VS
          M[i][j][k].vs  =  vs[el][k][j][i];
          #endif

          #if RHO
          M[i][j][k].rho = rho[el][k][j][i];
          #endif

          #if VPV
          M[i][j][k].vpv = vpv[el][k][j][i];
          #endif

          #if VPH
          M[i][j][k].vph = vph[el][k][j][i];
          #endif

          #if VSV
          M[i][j][k].vsv = vsv[el][k][j][i];
          #endif

          #if VSH
          M[i][j][k].vsh = vsh[el][k][j][i];
          #endif

          #if ETA
          M[i][j][k].eta = eta[el][k][j][i];
          #endif

          #if QMU
          M[i][j][k].qmu = qmu[el][k][j][i];
          #endif

          if (x < xmin)     xmin = x;
          if (x > xmax)     xmax = x;
          if (y < ymin)     ymin = y;
          if (y > ymax)     ymax = y;
          if (z < zmin)     zmin = z;
          if (z > zmax)     zmax = z;
          if (r < rmin)     rmin = r;
          if (r > rmax)     rmax = r;
          if (theta < tmin) tmin = theta;
          if (theta > tmax) tmax = theta;
          if (phi < pmin)   pmin = phi;
          if (phi > pmax)   pmax = phi;
        }

    if (xmin  <= gb->xmax + GB_TOLERANCE && xmax  >= gb->xmin - GB_TOLERANCE &&
        ymin  <= gb->ymax + GB_TOLERANCE && ymax  >= gb->ymin - GB_TOLERANCE &&
        zmin  <= gb->zmax + GB_TOLERANCE && zmax  >= gb->zmin - GB_TOLERANCE &&
        rmin  <= gb->rmax + GB_TOLERANCE && rmax  >= gb->rmin - GB_TOLERANCE &&
        tmin  <= gb->tmax + GB_TOLERANCE && tmax  >= gb->tmin - GB_TOLERANCE &&
       ((pmin <= gb->pmax + GB_TOLERANCE && pmax  >= gb->pmin - GB_TOLERANCE) ||
        (fabsl (ymin) <= GB_TOLERANCE || fabsl (ymax) <= GB_TOLERANCE)))
    {
      insertElNode (T, lle, nel);
      insertPmNode (M, llm);
    }
  }

  return 0;
}

unsigned scanMeshAndModel (int ic, char *prm, unsigned nelm, struct Boundaries *gb,
                           struct ElNode *lle, struct PmNode *llm, unsigned *nel)
{
  /* Reads the mesh and all model parameters */
  for (unsigned iic = 0; iic < NC; iic++)
  {
    MPI_Barrier (MPI_COMM_WORLD);

    if (ic == 0)
    {
      fprintf (stderr, "Reading topological file 'proc%06d_reg1_solver_data.bin'...\n", iic);

      #if VP
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vp.bin'...\n", iic);
      #endif

      #if VS
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vs.bin'...\n", iic);
      #endif

      #if RHO
      fprintf (stderr, "Reading model file 'proc%06d_reg1_rho.bin'...\n", iic);
      #endif

      #if VPV
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vpv.bin'...\n", iic);
      #endif

      #if VPH
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vph.bin'...\n", iic);
      #endif

      #if VSV
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vsv.bin'...\n", iic);
      #endif

      #if VSH
      fprintf (stderr, "Reading model file 'proc%06d_reg1_vsh.bin'...\n", iic);
      #endif

      #if ETA
      fprintf (stderr, "Reading model file 'proc%06d_reg1_eta.bin'...\n", iic);
      #endif

      #if QMU
      fprintf (stderr, "Reading model file 'proc%06d_reg1_qmu.bin'...\n", iic);
      #endif
    }

    unsigned r = readInputMeshAndModel (iic, prm, nelm, gb, lle, llm, nel); if (r) return r;
  }

  return 0;
}

void toArrayElAndPmNodes (unsigned nel, struct ElNode *lle, struct PmNode *llm,
                          struct Point Ti[nel][NX][NY][NZ],
                          struct Parameters M[nel][NX][NY][NZ])
{
  /* Converts the linked lists to arrays */
  struct ElNode *p = lle->next;
  struct PmNode *q = llm->next;

  for (unsigned el = 0; el < nel; el++)
  {
    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          Ti[el][i][j][k].x     = p->T[i][j][k].x;
          Ti[el][i][j][k].y     = p->T[i][j][k].y;
          Ti[el][i][j][k].z     = p->T[i][j][k].z;

          Ti[el][i][j][k].r     = p->T[i][j][k].r;
          Ti[el][i][j][k].theta = p->T[i][j][k].theta;
          Ti[el][i][j][k].phi   = p->T[i][j][k].phi;

          #if VP
          M[el][i][j][k].vp  = q->M[i][j][k].vp;
          #endif

          #if VS
          M[el][i][j][k].vs  = q->M[i][j][k].vs;
          #endif

          #if RHO
          M[el][i][j][k].rho = q->M[i][j][k].rho;
          #endif

          #if VPV
          M[el][i][j][k].vpv = q->M[i][j][k].vpv;
          #endif

          #if VPH
          M[el][i][j][k].vph = q->M[i][j][k].vph;
          #endif

          #if VSV
          M[el][i][j][k].vsv = q->M[i][j][k].vsv;
          #endif

          #if VSH
          M[el][i][j][k].vsh = q->M[i][j][k].vsh;
          #endif

          #if ETA
          M[el][i][j][k].eta = q->M[i][j][k].eta;
          #endif

          #if QMU
          M[el][i][j][k].qmu = q->M[i][j][k].qmu;
          #endif
        }

    p = p->next;
    q = q->next;
  }

  deleteElNode (lle);
  deletePmNode (llm);
}

unsigned readOutputMesh (int ic, char *prm, struct Boundaries *gb,
                         struct Point To[NEL2][NX][NY][NZ])
{
  /* Reads output file mesh */
  unsigned N = NEL2 * NX * NY * NZ;

  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s/proc%06d_reg1_solver_data.bin", prm, ic) < 33) return 1;

  FILE *file = fopen (filename, "rb");

  if (file == NULL) return 2;

  unsigned junk;
  unsigned nspec;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nspec, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned nval;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nval, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xf[NG2];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xf, sizeof (float), NG2, file) != NG2) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float yf[NG2];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (yf, sizeof (float), NG2, file) != NG2) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float zf[NG2];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (zf, sizeof (float), NG2, file) != NG2) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ibool[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ibool, sizeof (unsigned), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned idoubling[NEL2];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (idoubling, sizeof (unsigned), NEL2, file) != NEL2) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ispec_is_tiso[NEL2];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ispec_is_tiso, sizeof (unsigned), NEL2, file) != NEL2) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xix[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xix, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiy[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiy, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiz[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etax[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etay[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etaz[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammax[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammay[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammaz[NEL2][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);

  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;
  gb->rmin =  INFINITY;
  gb->rmax = -INFINITY;
  gb->tmin =  INFINITY;
  gb->tmax = -INFINITY;
  gb->pmin =  INFINITY;
  gb->pmax = -INFINITY;

  for (unsigned el = 0; el < NEL2; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          unsigned ig = ibool[el][k][j][i] - 1;

          long double x = xf[ig];
          long double y = yf[ig];
          long double z = zf[ig];

          To[el][i][j][k].x = x;
          To[el][i][j][k].y = y;
          To[el][i][j][k].z = z;

          if (x < gb->xmin) gb->xmin = x;
          if (x > gb->xmax) gb->xmax = x;
          if (y < gb->ymin) gb->ymin = y;
          if (y > gb->ymax) gb->ymax = y;
          if (z < gb->zmin) gb->zmin = z;
          if (z > gb->zmax) gb->zmax = z;

          long double r, theta, phi;

          xYZ2RThetaPhil (x, y, z, &r, &theta, &phi);

          To[el][i][j][k].r     = r;
          To[el][i][j][k].theta = theta;
          To[el][i][j][k].phi   = phi;

          if (r < gb->rmin)     gb->rmin = r;
          if (r > gb->rmax)     gb->rmax = r;
          if (theta < gb->tmin) gb->tmin = theta;
          if (theta > gb->tmax) gb->tmax = theta;
          if (phi < gb->pmin)   gb->pmin = phi;
          if (phi > gb->pmax)   gb->pmax = phi;
        }

  return 0;
}

unsigned readMeshModelAndJacobian (int ic, char *input, char *prm,
                                   struct SphericalPoint Tm[NEL][NX][NY][NZ],
                                   double M[NEL][NX][NY][NZ],
                                   double J[NEL][NX][NY][NZ])
{
  /* Reads mesh and model files and computes the Jacobian */
  unsigned N = NEL * NX * NY * NZ;

  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s/proc%06d_reg1_solver_data.bin",
               input, ic) < 25) return 1;

  FILE *file = fopen (filename, "rb");

  if (file == NULL) return 2;

  unsigned junk;
  unsigned nspec;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nspec, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned nval;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&nval, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float x[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&x, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float y[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&y, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float z[NG];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&z, sizeof (float), NG, file) != NG) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ibool[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&ibool, sizeof (unsigned), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned idoubling[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&idoubling, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ispec_is_tiso[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&ispec_is_tiso, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xix[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&xix, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiy[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&xiy, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&xiz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&etax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&etay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&etaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&gammax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&gammay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (&gammaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);

  if (sprintf (filename, "%s/proc%06d_reg1_%s.bin",
               input, ic, prm) < 15) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;

  float tM[NEL][NZ][NY][NX];

  if (fread (tM, sizeof (float), N, file) != N) return 5;

  fclose (file);

  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          M[el][i][j][k] = tM[el][k][j][i];

          unsigned ig = ibool[el][k][j][i] - 1;

          double xi = x[ig];
          double yi = y[ig];
          double zi = z[ig];

          double r, theta, phi;

          xYZ2RThetaPhi (xi, yi, zi, &r, &theta, &phi);

          Tm[el][i][j][k].r     = roundV (r);
          Tm[el][i][j][k].theta = roundV (theta);
          Tm[el][i][j][k].phi   = roundV (phi);

          double xixl = xix[el][k][j][i];
          double xiyl = xiy[el][k][j][i];
          double xizl = xiz[el][k][j][i];

          double etaxl = etax[el][k][j][i];
          double etayl = etay[el][k][j][i];
          double etazl = etaz[el][k][j][i];

          double gammaxl = gammax[el][k][j][i];
          double gammayl = gammay[el][k][j][i];
          double gammazl = gammaz[el][k][j][i];

          double k1 = xixl * (etayl * gammazl - etazl * gammayl);
          double k2 = xiyl * (etaxl * gammazl - etazl * gammaxl);
          double k3 = xizl * (etaxl * gammayl - etayl * gammaxl);

          J[el][i][j][k] = 1.0 / (k1 - k2 + k3);
        }

  return 0;
}

unsigned readBSplinesHeader (double *rmin, double *rmax,
                             unsigned *nn, unsigned *dg,
                             unsigned zone)
{
  /* Reads the header of the B-splines file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "knots_Zone%u.dat", zone) < 15) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  int f1 = fscanf (file, "%*[^\n]\n");
  int f2 = fscanf (file, "%u", nn);
  int f3 = fscanf (file, "%u", dg);
  int f4 = fscanf (file, "%lf", rmin);
  int f5 = fscanf (file, "%lf", rmax);

  fclose (file);

  if (f1)      return 3;
  if (f2 != 1) return 4;
  if (f3 != 1) return 5;
  if (f4 != 1) return 6;
  if (f5 != 1) return 7;

  return 0;
}

unsigned readKnots (double r1, double r2, unsigned dg,
                    unsigned zone, unsigned nnt, double T[nnt])
{
  /* Reads the B-splines knots */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "knots_Zone%u.dat", zone) < 5) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  for (unsigned i = 0; i <= dg; i++)
  {
    T[i] = r1;
    T[nnt - i - 1] = r2;
  }

  int f1 = fscanf (file, "%*[^\n]\n");
  int f2 = fscanf (file, "%*[^\n]\n");
  int f3 = fscanf (file, "%*[^\n]\n");
  int f4 = fscanf (file, "%*[^\n]\n");
  int f5 = fscanf (file, "%*[^\n]\n");

  if (f1 || f2 || f3 || f4 || f5)
  {
    fclose (file); return 3;
  }

  for (unsigned i = dg + 1; i < nnt - dg - 1; i++)

    if (fscanf (file, "%lf", &T[i]) != 1)
    {
      fclose (file); return 4;
    }

  fclose (file);

  return 0;
}

unsigned readMeanModelHeader (char *prm, unsigned *nl)
{
  /* Reads the header of the mean model file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s.dat", prm) < 5) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  int f1 = fscanf (file, "%*s");
  int f2 = fscanf (file, "%u", nl);

  fclose (file);

  if (f1)      return 3;
  if (f2 != 1) return 4;

  return 0;
}

unsigned readMeanModel (char *prm, unsigned nl,
                        struct MeanModel Mm[nl])
{
  /* Reads the mean model file */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s.dat", prm) < 5) return 1;

  FILE *file = fopen (filename, "r");

  if (file == NULL) return 2;

  if (fscanf (file, "%*[^\n]\n") != 0) return 3;
  if (fscanf (file, "%*[^\n]\n") != 0) return 3;

  double depth, vam, vgm;

  for (unsigned l = 0; l < nl; l++)
  {
    if (fscanf (file, "%lf %lf %lf", &depth, &vam, &vgm) != 3) return 3;

    Mm[nl - l - 1].r   = depth2R (depth);
    Mm[nl - l - 1].vam = vam;
    Mm[nl - l - 1].vgm = vgm;
  }

  fclose (file);

  return 0;
}

double aMeanModel (double r, unsigned nl, struct MeanModel Mm[nl])
{
  /* Computes the arithmetic mean model value using
     linear interpolation */
  int n = (int) nl;
  int i = (int) (n * (r - Mm[0].r) / (Mm[n - 1].r - Mm[0].r) + 0.5);

  if (i <= 0) return Mm[0].vam; else if (i > n - 1) return Mm[n - 1].vam;

  double m = (Mm[i].vam - Mm[i - 1].vam) / (Mm[i].r - Mm[i - 1].r);

  return Mm[i].vam + m * (r - Mm[i].r);
}

unsigned computeGatheringArrays (int ic, int nc,
                                 unsigned n1, unsigned n2, unsigned n3,
                                 int k[nc], int kk[nc])
{
  /* Computes arrays for MPI_Gatherv */
  int n1_l = n1 / nc;
  int rmn  = n1 % nc;

  for (int i = 0; i < nc; i++)
  {
    k[i]  = (i == 0 || i > rmn) ? n1_l * n2 * n3 : (n1_l + 1) * n2 * n3;
    kk[i] = (i == 0) ? 0 : kk[i - 1] + k[i - 1];
  }

  return (ic == 0 || ic > rmn) ? n1_l : n1_l + 1;
}

void createOutputGridBk (int ic, int nc, int k[nc], int kk[nc],
                         long double r1, long double r2,
                         unsigned np, unsigned np_l, unsigned nt, unsigned nr,
                         struct Boundaries *gb, struct Point To[np_l][nt][nr])
{
  /* Creates the file grid */
  gb->rmin =  r1 <= r2 ? r1 : r2;
  gb->rmax =  r1 >= r2 ? r1 : r2;
  gb->tmin =  0.00;
  gb->tmax =  PI;

  long double dr = (gb->rmax - gb->rmin) / (nr - 1);
  long double dp = 2 * PI / (np - 1);
  long double dt = (gb->tmax - gb->tmin) / (nt - 1);

  unsigned alpha = (ic != nc - 1) ? 0 : 1;

  gb->pmin = -PI + dp * kk[ic] / (nt * nr);
  gb->pmax = -PI + dp * ((kk[ic] + k[ic]) / (nt * nr) - alpha);

  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;

  long double r = gb->rmin;

  for (unsigned k = 0; k < nr; k++)
  {
    long double theta = gb->tmin;

    for (unsigned j = 0; j < nt; j++)
    {
      long double phi = gb->pmin;

      for (unsigned i = 0; i < np_l; i++)
      {
        To[i][j][k].r     = r;
        To[i][j][k].theta = theta;
        To[i][j][k].phi   = phi;

        long double x, y, z;

        rThetaPhi2XYZl (r, theta, phi, &x, &y, &z);

        To[i][j][k].x = x;
        To[i][j][k].y = y;
        To[i][j][k].z = z;

        if (x < gb->xmin) gb->xmin = x;
        if (x > gb->xmax) gb->xmax = x;
        if (y < gb->ymin) gb->ymin = y;
        if (y > gb->ymax) gb->ymax = y;
        if (z < gb->zmin) gb->zmin = z;
        if (z > gb->zmax) gb->zmax = z;

        phi += dp;
      }

      theta += dt;
    }

    r += dr;
  }
}

static void crossProduct2 (long double x1, long double y1, long double z1,
                           long double x2, long double y2, long double z2,
                           long double *u, long double *v, long double *w)
{
  /* Computes the cross-product of two 3D arrays */
  *u = y1 * z2 - z1 * y2;
  *v = z1 * x2 - x1 * z2;
  *w = x1 * y2 - y1 * x2;

  long double norm = sqrt (squarel (*u) + squarel (*v) + squarel (*w));

  *u = *u / norm; *v = *v / norm; *w = *w / norm;
}

static void sliceBoundaries (int ic, int nc, int k[nc], int kk[nc],
                             unsigned nr, long double dd,
                             long double r1, long double r2,
                             long double t1, long double t2,
                             long double p1, long double p2,
                             struct Boundaries *gb)
{
  /* Updates the slice boundaries according to the number
     of cores to save computational time */
  long double x1, y1, z1;
  long double x2, y2, z2;

  rThetaPhi2XYZl (1.L, t1, p1, &x1, &y1, &z1);
  rThetaPhi2XYZl (1.L, t2, p2, &x2, &y2, &z2);

  long double u, v, w;

  crossProduct2 (x1, y1, z1, x2, y2, z2, &u, &v, &w);

  long double sind, cosd;
  long double xi, yi, zi;
  long double xf, yf, zf;

  sind = sinl (dd * kk[ic] / nr);
  cosd = cosl (dd * kk[ic] / nr);

  xi = x1 * (cosd + squarel (u) * (1 - cosd))
     + y1 * (u * v * (1 - cosd) - w * sind)
     + z1 * (u * w * (1 - cosd) + v * sind);

  yi = x1 * (v * u * (1 - cosd) + w * sind)
     + y1 * (cosd + squarel (v) * (1 - cosd))
     + z1 * (v * w * (1 - cosd) - u * sind);

  zi = x1 * (w * u * (1 - cosd) - v * sind)
     + y1 * (w * v * (1 - cosd) + u * sind)
     + z1 * (cosd + squarel (w) * (1 - cosd));

  unsigned alpha = (ic != nc - 1) ? 0 : 1;

  sind = sinl (dd * ((kk[ic] + k[ic]) / nr - alpha));
  cosd = cosl (dd * ((kk[ic] + k[ic]) / nr - alpha));

  xf = x1 * (cosd + squarel (u) * (1 - cosd))
     + y1 * (u * v * (1 - cosd) - w * sind)
     + z1 * (u * w * (1 - cosd) + v * sind);

  yf = x1 * (v * u * (1 - cosd) + w * sind)
     + y1 * (cosd + squarel (v) * (1 - cosd))
     + z1 * (v * w * (1 - cosd) - u * sind);

  zf = x1 * (w * u * (1 - cosd) - v * sind)
     + y1 * (w * v * (1 - cosd) + u * sind)
     + z1 * (cosd + squarel (w) * (1 - cosd));

  long double ri, ti, pi;
  long double rf, tf, pf;

  xYZ2RThetaPhil (xi, yi, zi, &ri, &ti, &pi);
  xYZ2RThetaPhil (xf, yf, zf, &rf, &tf, &pf);

  gb->rmin =  r1 <= r2 ? r1 : r2;
  gb->rmax =  r1 >= r2 ? r1 : r2;
  gb->tmin =  ti <= tf ? ti : tf;
  gb->tmax =  ti >= tf ? ti : tf;
  gb->pmin =  pi <= pf ? pi : pf;
  gb->pmax =  pi >= pf ? pi : pf;
}

void createOutputGridDD (int ic, int nc, int k[nc], int kk[nc],
                         unsigned nd, unsigned nd_l, unsigned nr,
                         long double r1, long double r2,
                         long double t1, long double t2,
                         long double p1, long double p2,
                         long double R[nr], long double Delta_l[nd_l],
                         struct Boundaries *gb, struct Point To[nd_l][nr])
{
  /* Creates the file grid */
  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;

  long double gcarc = vincentyl (t1, p1, t2, p2);

  long double dr = (r2 - r1) / (nr - 1);
  long double dd = gcarc / (nd - 1);

  sliceBoundaries (ic, nc, k, kk, nr, dd, r1, r2, t1, t2, p1, p2, gb);

  long double x1, y1, z1;
  long double x2, y2, z2;

  rThetaPhi2XYZl (1.L, t1, p1, &x1, &y1, &z1);
  rThetaPhi2XYZl (1.L, t2, p2, &x2, &y2, &z2);

  long double u, v, w;

  crossProduct2 (x1, y1, z1, x2, y2, z2, &u, &v, &w);

  for (unsigned j = 0; j < nr; j++)
  {
    R[j] = r1 + j * dr;

    for (unsigned i = 0; i < nd_l; i++)
    {
      long double delta = dd * (i + kk[ic] / nr);

      Delta_l[i] = rad2Degreel (delta);

      long double sind = sinl (delta);
      long double cosd = cosl (delta);

      long double x = x1 * (cosd + squarel (u) * (1 - cosd))
                    + y1 * (u * v * (1 - cosd) - w * sind)
                    + z1 * (u * w * (1 - cosd) + v * sind);

      long double y = x1 * (v * u * (1 - cosd) + w * sind)
                    + y1 * (cosd + squarel (v) * (1 - cosd))
                    + z1 * (v * w * (1 - cosd) - u * sind);

      long double z = x1 * (w * u * (1 - cosd) - v * sind)
                    + y1 * (w * v * (1 - cosd) + u * sind)
                    + z1 * (cosd + squarel (w) * (1 - cosd));

      x *= R[j]; y *= R[j]; z *= R[j];

      long double r, theta, phi;

      xYZ2RThetaPhil (x, y, z, &r, &theta, &phi);

      To[i][j].x     = x;
      To[i][j].y     = y;
      To[i][j].z     = z;
      To[i][j].r     = r;
      To[i][j].theta = theta;
      To[i][j].phi   = phi;

      if (x < gb->xmin) gb->xmin = x;
      if (x > gb->xmax) gb->xmax = x;
      if (y < gb->ymin) gb->ymin = y;
      if (y > gb->ymax) gb->ymax = y;
      if (z < gb->zmin) gb->zmin = z;
      if (z > gb->zmax) gb->zmax = z;
    }
  }
}

void createOutputGridLL (int ic, int nc, int k[nc], int kk[nc],
                         long double r,
                         unsigned np, unsigned np_l, unsigned nt,
                         struct Boundaries *gb, struct Point To[np_l][nt])
{
  /* Creates the file grid */
  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;
  gb->rmin =  r - BOUNDARY_RATIO_4;
  gb->rmax =  r + BOUNDARY_RATIO_4;
  gb->tmin =  0.L;
  gb->tmax =  PI;

  long double dt = (gb->tmax - gb->tmin) / (nt - 1);
  long double dp = 2.0L * PI / (np - 1);

  unsigned alpha = (ic != nc - 1) ? 0 : 1;

  gb->pmin = -PI + dp * kk[ic] / nt;
  gb->pmax = -PI + dp * ((kk[ic] + k[ic]) / nt - alpha);

  long double theta = gb->tmin;

  for (unsigned j = 0; j < nt; j++)
  {
    long double phi = gb->pmin;

    for (unsigned i = 0; i < np_l; i++)
    {
      To[i][j].r     = r;
      To[i][j].theta = theta;
      To[i][j].phi   = phi;

      long double x, y, z;

      rThetaPhi2XYZl (r, theta, phi, &x, &y, &z);

      To[i][j].x = x;
      To[i][j].y = y;
      To[i][j].z = z;

      if (x < gb->xmin) gb->xmin = x;
      if (x > gb->xmax) gb->xmax = x;
      if (y < gb->ymin) gb->ymin = y;
      if (y > gb->ymax) gb->ymax = y;
      if (z < gb->zmin) gb->zmin = z;
      if (z > gb->zmax) gb->zmax = z;

      phi += dp;
    }

    theta += dt;
  }
}

void createOutputGrid1D (int ic, int nc, int k[nc], int kk[nc],
                         long double r1, long double r2,
                         unsigned np, unsigned np_l, unsigned nt, unsigned nr,
                         struct Boundaries *gb,
                         long double R[nr],
                         long double Theta[nt],
                         long double Phi[np_l])
{
  /* Creates the file grid */
  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;
  gb->rmin =  r1 <= r2 ? r1 : r2;
  gb->rmax =  r1 >= r2 ? r1 : r2;
  gb->tmin =  0.00;
  gb->tmax =  PI;

  long double dr = (gb->rmax - gb->rmin) / (nr - 1);
  long double dp = 2 * PI / (np - 1);
  long double dt = (gb->tmax - gb->tmin) / (nt - 1);

  unsigned alpha = (ic != nc - 1) ? 0 : 1;

  gb->pmin = -PI + dp * kk[ic] / (nt * nr);
  gb->pmax = -PI + dp * ((kk[ic] + k[ic]) / (nt * nr) - alpha);

  long double r = gb->rmin;

  for (unsigned k = 0; k < nr; k++)
  {
    R[k] = r;

    long double theta = gb->tmin;

    for (unsigned j = 0; j < nt; j++)
    {
      Theta[j] = theta;

      long double phi = gb->pmin;

      for (unsigned i = 0; i < np_l; i++)
      {
        Phi[i] = phi;

        long double x, y, z;

        rThetaPhi2XYZl (r, theta, phi, &x, &y, &z);

        if (x < gb->xmin) gb->xmin = x;
        if (x > gb->xmax) gb->xmax = x;
        if (y < gb->ymin) gb->ymin = y;
        if (y > gb->ymax) gb->ymax = y;
        if (z < gb->zmin) gb->zmin = z;
        if (z > gb->zmax) gb->zmax = z;

        phi += dp;
      }

      theta += dt;
    }

    r += dr;
  }
}

void createOutputGridPf (long double r1, long double r2,
                         long double t, long double p,
                         unsigned nr,
                         long double R[nr],
                         struct Boundaries *gb, struct Point To[nr])
{
  /* Creates the file grid */
  gb->xmin =  INFINITY;
  gb->xmax = -INFINITY;
  gb->ymin =  INFINITY;
  gb->ymax = -INFINITY;
  gb->zmin =  INFINITY;
  gb->zmax = -INFINITY;
  gb->rmin =  r1 <= r2 ? r1 : r2;
  gb->rmax =  r1 >= r2 ? r1 : r2;
  gb->tmin =  t - GB_TOLERANCE;
  gb->tmax =  t + GB_TOLERANCE;
  gb->pmin =  p - GB_TOLERANCE;
  gb->pmax =  p + GB_TOLERANCE;

  long double xr, yr, zr;

  rThetaPhi2XYZl (1.L, t, p, &xr, &yr, &zr);

  long double dr = (r2 - r1) / (nr - 1);

  for (unsigned i = 0; i < nr; i++)
  {
    R[i] = r1 + i * dr;

    long double x = xr * R[i], y = yr * R[i], z = zr * R[i];

    long double r, theta, phi;

    xYZ2RThetaPhil (x, y, z, &r, &theta, &phi);

    To[i].x     = x;
    To[i].y     = y;
    To[i].z     = z;
    To[i].r     = r;
    To[i].theta = theta;
    To[i].phi   = phi;

    if (x < gb->xmin) gb->xmin = x;
    if (x > gb->xmax) gb->xmax = x;
    if (y < gb->ymin) gb->ymin = y;
    if (y > gb->ymax) gb->ymax = y;
    if (z < gb->zmin) gb->zmin = z;
    if (z > gb->zmax) gb->zmax = z;
  }
}

static inline char machineEndianess (void)
{
  /* Returns the machine endianness */
  int c = 1; return (*(char*) &c == 1) ? '<' : '>';
}

unsigned writeModelBk (char *prm,
                       unsigned np, unsigned nt, unsigned nr,
                       struct Parameters Mo[np][nt][nr])
{
  /* Writes the crustal block model */
  unsigned n = np * nt * nr;

  char e = machineEndianess ();
  char name[MAX_STRING_LEN];

  #if VP
  float  vp[np][nt][nr];
  #endif

  #if VS
  float  vs[np][nt][nr];
  #endif

  #if RHO
  float rho[np][nt][nr];
  #endif

  #if VPV
  float vpv[np][nt][nr];
  #endif

  #if VPH
  float vph[np][nt][nr];
  #endif

  #if VSV
  float vsv[np][nt][nr];
  #endif

  #if VSH
  float vsh[np][nt][nr];
  #endif

  #if ETA
  float eta[np][nt][nr];
  #endif

  #if QMU
  float qmu[np][nt][nr];
  #endif

  for (unsigned i = 0; i < np; i++)

    for (unsigned j = 0; j < nt; j++)

      for (unsigned k = 0; k < nr; k++)
      {
        #if VP
        vp[i][j][k]  = (float) Mo[i][j][k].vp;
        if (vp[i][j][k] < 0)  vp[i][j][k] = 0;
        #endif

        #if VS
        vs[i][j][k]  = (float) Mo[i][j][k].vs;
        if (vs[i][j][k] < 0)  vs[i][j][k] = 0;
        #endif

        #if RHO
        rho[i][j][k] = (float) Mo[i][j][k].rho;
        if (rho[i][j][k] < 0) rho[i][j][k] = 0;
        #endif

        #if VPV
        vpv[i][j][k] = (float) Mo[i][j][k].vpv;
        if (vpv[i][j][k] < 0) vpv[i][j][k] = 0;
        #endif

        #if VPH
        vph[i][j][k] = (float) Mo[i][j][k].vph;
        if (vph[i][j][k] < 0) vph[i][j][k] = 0;
        #endif

        #if VSV
        vsv[i][j][k] = (float) Mo[i][j][k].vsv;
        if (vsv[i][j][k] < 0) vsv[i][j][k] = 0;
        #endif

        #if VSH
        vsh[i][j][k] = (float) Mo[i][j][k].vsh;
        if (vsh[i][j][k] < 0) vsh[i][j][k] = 0;
        #endif

        #if ETA
        eta[i][j][k] = (float) Mo[i][j][k].eta;
        if (eta[i][j][k] < 0) eta[i][j][k] = 0;
        #endif

        #if QMU
        qmu[i][j][k] = (float) Mo[i][j][k].qmu;
        if (qmu[i][j][k] < 0) qmu[i][j][k] = 0;
        #endif
      }

  FILE *file;

  #if VP
  if (sprintf (name, "%s/vp.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vp, sizeof (float), n, file);
  #endif

  #if VS
  if (sprintf (name, "%s/vs.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vs, sizeof (float), n, file);

  fclose (file);
  #endif

  #if RHO
  if (sprintf (name, "%s/rho.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (rho, sizeof (float), n, file);

  fclose (file);
  #endif

  #if VPV
  if (sprintf (name, "%s/vpv.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vpv, sizeof (float), n, file);

  fclose (file);
  #endif

  #if VPH
  if (sprintf (name, "%s/vph.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vph, sizeof (float), n, file);

  fclose (file);
  #endif

  #if VSV
  if (sprintf (name, "%s/vsv.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vsv, sizeof (float), n, file);

  fclose (file);
  #endif

  #if VSH
  if (sprintf (name, "%s/vsh.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (vsh, sizeof (float), n, file);

  fclose (file);
  #endif

  #if ETA
  if (sprintf (name, "%s/eta.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (eta, sizeof (float), n, file);

  fclose (file);
  #endif

  #if QMU
  if (sprintf (name, "%s/qmu.bin", prm) < 8) return 1;

  file = fopen (name, "wb");

  if (file == NULL) return 2;

  fwrite (&e, sizeof (char), 1, file);
  fwrite (&np, sizeof (unsigned), 1, file);
  fwrite (&nt, sizeof (unsigned), 1, file);
  fwrite (&nr, sizeof (unsigned), 1, file);
  fwrite (qmu, sizeof (float), n, file);

  fclose (file);
  #endif

  return 0;
}

unsigned writeModelDD (char *argv[], unsigned nd, unsigned nr,
                       long double R[nr], long double Delta[nd],
                       struct Parameters Mo[nd][nr])
{
  /* Writes vertical cross section */
  char name[MAX_STRING_LEN];

  FILE *file;

  #if VP
  if (sprintf (name, "%s/vp_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vp (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vp);

  fclose (file);
  #endif

  #if VS
  if (sprintf (name, "%s/vs_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vs (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vs);

  fclose (file);
  #endif

  #if RHO
  if (sprintf (name, "%s/rho_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Rho (g/cm^3)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].rho);

  fclose (file);
  #endif

  #if VPV
  if (sprintf (name, "%s/vpv_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vpv (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vpv);

  fclose (file);
  #endif

  #if VPH
  if (sprintf (name, "%s/vph_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vph (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vph);

  fclose (file);
  #endif

  #if VSV
  if (sprintf (name, "%s/vsv_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vsv (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vsv);

  fclose (file);
  #endif

  #if VSH
  if (sprintf (name, "%s/vsh_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "#radius (km)   delta (degrees)    Vsh (km/s)\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].vsh);

  fclose (file);
  #endif

  #if ETA
  if (sprintf (name, "%s/eta_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "##radius (km)   delta (degrees)       Eta\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].eta);

  fclose (file);
  #endif

  #if QMU
  if (sprintf (name, "%s/qmu_VCS.dat", argv[10]) < 12) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#nrad ndel: %u %u\n", nr, nd);
  fprintf (file, "#lat1 lon1 lat2 lon2: %s %s %s %s\n",
           argv[3], argv[4], argv[5], argv[6]);
  fprintf (file, "##radius (km)   delta (degrees)       Qmu\n");

  for (unsigned i = 0; i < nr; i++)

    for (unsigned j = 0; j < nd; j++)

      fprintf (file, "%10.3Lf %14.3Lf %19LE\n", R[i] * EARTH_R,
                                                Delta[j],
                                                Mo[j][i].qmu);

  fclose (file);
  #endif

  return 0;
}

unsigned writeModelGLL (int ic, char *prm,
                        struct Parameters M[NEL2][NX][NY][NZ])
{
  /* Writes the new model */
  unsigned N = NEL2 * NX * NY * NZ;
  unsigned size = 4 * N;

  char filename[MAX_STRING_LEN];

  #if VP
  float  vp[NEL2][NZ][NY][NX];
  #endif

  #if VS
  float  vs[NEL2][NZ][NY][NX];
  #endif

  #if RHO
  float rho[NEL2][NZ][NY][NX];
  #endif

  #if VPV
  float vpv[NEL2][NZ][NY][NX];
  #endif

  #if VPH
  float vph[NEL2][NZ][NY][NX];
  #endif

  #if VSV
  float vsv[NEL2][NZ][NY][NX];
  #endif

  #if VSH
  float vsh[NEL2][NZ][NY][NX];
  #endif

  #if ETA
  float eta[NEL2][NZ][NY][NX];
  #endif

  #if QMU
  float qmu[NEL2][NZ][NY][NX];
  #endif

  for (unsigned el = 0; el < NEL2; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          #if VP
          vp[el][k][j][i]  = (float) M[el][i][j][k].vp;
          #endif

          #if VS
          vs[el][k][j][i]  = (float) M[el][i][j][k].vs;
          #endif

          #if RHO
          rho[el][k][j][i] = (float) M[el][i][j][k].rho;
          #endif

          #if VPV
          vpv[el][k][j][i] = (float) M[el][i][j][k].vpv;
          #endif

          #if VPH
          vph[el][k][j][i] = (float) M[el][i][j][k].vph;
          #endif

          #if VSV
          vsv[el][k][j][i] = (float) M[el][i][j][k].vsv;
          #endif

          #if VSH
          vsh[el][k][j][i] = (float) M[el][i][j][k].vsh;
          #endif

          #if ETA
          eta[el][k][j][i] = (float) M[el][i][j][k].eta;
          #endif

          #if QMU
          qmu[el][k][j][i] = (float) M[el][i][j][k].qmu;
          #endif
        }

  FILE *file;

  #if VP
  if (sprintf (filename, "%s/proc%06d_reg1_vp.bin", prm, ic) < 23) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vp, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if VS
  if (sprintf (filename, "%s/proc%06d_reg1_vs.bin", prm, ic) < 23) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vs, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if RHO
  if (sprintf (filename, "%s/proc%06d_reg1_rho.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (rho, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if VPV
  if (sprintf (filename, "%s/proc%06d_reg1_vpv.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vpv, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if VPH
  if (sprintf (filename, "%s/proc%06d_reg1_vph.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vph, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if VSV
  if (sprintf (filename, "%s/proc%06d_reg1_vsv.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vsv, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if VSH
  if (sprintf (filename, "%s/proc%06d_reg1_vsh.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (vsh, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if ETA
  if (sprintf (filename, "%s/proc%06d_reg1_eta.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (eta, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if QMU
  if (sprintf (filename, "%s/proc%06d_reg1_qmu.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (qmu, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  return 0;
}

unsigned writeModelLL (char *prm, struct Boundaries *gb,
                       long double depth, unsigned np, unsigned nt,
                       struct Parameters Mo[np][nt])
{
  /* Writes depth slice */
  char name[MAX_STRING_LEN];

  double Theta[nt], Phi[np];

  double t = gb->tmin;
  double p = gb->pmin;

  double dt = (gb->tmax - gb->tmin) / (nt - 1);
  double dp = 2 * PI / (np - 1);

  for (unsigned i = 0; i < nt; i++)
  {
    Theta[nt - i - 1] = 90 - rad2Degreel (t); t += dt;
  }

  for (unsigned j = 0; j < np; j++)
  {
    Phi[j] = rad2Degreel (p); p += dp;
  }

  FILE *file;

  #if VP
  if (sprintf (name, "%s/vp_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vp (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vp);

  fclose (file);
  #endif

  #if VS
  if (sprintf (name, "%s/vs_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vs (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vs);

  fclose (file);
  #endif

  #if RHO
  if (sprintf (name, "%s/rho_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Rho (g/cm^3)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].rho);

  fclose (file);
  #endif

  #if VPV
  if (sprintf (name, "%s/vpv_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vpv (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vpv);

  fclose (file);
  #endif

  #if VPH
  if (sprintf (name, "%s/vph_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vph (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vph);

  fclose (file);
  #endif

  #if VSV
  if (sprintf (name, "%s/vsv_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vsv (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vsv);

  fclose (file);
  #endif

  #if VSH
  if (sprintf (name, "%s/vsh_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)    Vsh (km/s)\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].vsh);

  fclose (file);
  #endif

  #if ETA
  if (sprintf (name, "%s/eta_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)       Eta\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].eta);

  fclose (file);
  #endif

  #if QMU
  if (sprintf (name, "%s/qmu_%Lg_DS.dat", prm, depth) < 10) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#depth (km) nlat nlon: %Lg %u %u\n", depth, nt, np);
  fprintf (file, "#latitude (degrees)   longitude (degrees)       Qmu\n");

  for (unsigned i = 0; i < nt; i++)

    for (unsigned j = 0; j < np; j++)

      fprintf (file, "%15.7lf %21.7lf %18LE\n", Theta[i], Phi[j],
                                                Mo[j][nt - i - 1].qmu);

  fclose (file);
  #endif

  return 0;
}

unsigned writeModelPf (char *argv[], unsigned nr,
                       long double R[nr], struct Parameters Mo[nr])
{
  /* Writes vertical cross section */
  char name[MAX_STRING_LEN];

  FILE *file;

  #if VP
  if (sprintf (name, "%s/vp_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vp (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vp);

  fclose (file);
  #endif

  #if VS
  if (sprintf (name, "%s/vs_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vs (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vs);

  fclose (file);
  #endif

  #if RHO
  if (sprintf (name, "%s/rho_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Rho (g/cm^3)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].rho);

  fclose (file);
  #endif

  #if VPV
  if (sprintf (name, "%s/vpv_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vpv (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vpv);

  fclose (file);
  #endif

  #if VPH
  if (sprintf (name, "%s/vph_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vph (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vph);

  fclose (file);
  #endif

  #if VSV
  if (sprintf (name, "%s/vsv_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vsv (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vsv);

  fclose (file);
  #endif

  #if VSH
  if (sprintf (name, "%s/vsh_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)     Vsh (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].vsh);

  fclose (file);
  #endif

  #if ETA
  if (sprintf (name, "%s/eta_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)        Eta\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].eta);

  fclose (file);
  #endif

  #if QMU
  if (sprintf (name, "%s/qmu_PF.dat", argv[7]) < 11) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#lat lon: %s %s\n", argv[3], argv[4]);
  fprintf (file, "#depth (km)        Qmu\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%9.3Lf %17LE\n", r2Depthl (R[i]), Mo[i].qmu);

  fclose (file);
  #endif

  return 0;
}

unsigned writeModel1D (char *prm, unsigned nr,
                       long double R[nr], struct Means Mo[nr])
{
  /* Writes vertical cross section */
  char name[MAX_STRING_LEN];

  FILE *file;

  #if VP
  if (sprintf (name, "%s/vp.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)    Vp_arithmetic (km/s)    Vp_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vp, Mo[i].g.vp);

  fclose (file);
  #endif

  #if VS
  if (sprintf (name, "%s/vs.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)    Vs_arithmetic (km/s)    Vs_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vs, Mo[i].g.vs);

  fclose (file);
  #endif

  #if RHO
  if (sprintf (name, "%s/rho.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)   Rho_arithmetic (g/cm^3)   Rho_geometric (g/cm^3)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.rho, Mo[i].g.rho);

  fclose (file);
  #endif

  #if VPV
  if (sprintf (name, "%s/vpv.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)   Vpv_arithmetic (km/s)   Vpv_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vpv, Mo[i].g.vpv);

  fclose (file);
  #endif

  #if VPH
  if (sprintf (name, "%s/vph.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)   Vph_arithmetic (km/s)   Vph_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vph, Mo[i].g.vph);

  fclose (file);
  #endif

  #if VSV
  if (sprintf (name, "%s/vsv.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)   Vsv_arithmetic (km/s)   Vsv_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vsv, Mo[i].g.vsv);

  fclose (file);
  #endif

  #if VSH
  if (sprintf (name, "%s/vsh.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)   Vsh_arithmetic (km/s)   Vsh_geometric (km/s)\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.vsh, Mo[i].g.vsh);

  fclose (file);
  #endif

  #if ETA
  if (sprintf (name, "%s/eta.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)      Eta_arithmetic          Eta_geometric\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.eta, Mo[i].g.eta);

  fclose (file);
  #endif

  #if QMU
  if (sprintf (name, "%s/qmu.dat", prm) < 8) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#ndep: %u\n", nr);
  fprintf (file, "#depth (km)      Qmu_arithmetic          Qmu_geometric\n");

  for (int i = nr - 1; i >= 0; i--)

    fprintf (file, "%8.1Lf %20.4LE %23.4LE\n", r2Depthl (R[i]),
                                               Mo[i].a.qmu, Mo[i].g.qmu);

  fclose (file);
  #endif

  return 0;
}

unsigned writeCoefficients (char *output, char *prm, unsigned zone,
                            unsigned nmax, unsigned ns, unsigned nlg,
                            double A[ns][nlg], double B[ns][nlg])
{
  /* Writes coefficients */
  char name[MAX_STRING_LEN];

  if (sprintf (name, "%s/mns_Z%u_%s.dat", output, zone, prm) < 15) return 1;

  FILE *file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "#radial and spherical harmonics coefficients"
                 "  (A^m_n)_s  (B^m_n)_s  (m <= n <= Nmax)"
                 "  (s < NS)  Nmax = %u  NS = %u\n", nmax, ns);

  for (unsigned s = 0; s < ns; s++)

    for (unsigned n = 0; n <= nmax; n++)

      for (unsigned m = 0; m <= n; m++)
      {
        unsigned mni = mN2I (m, n);

        sprintf (name, "(A^%u_%u)_%u", m, n, s);
        fprintf (file, " %-15s = %13E\n", name, A[s][mni]);
        sprintf (name, "(B^%u_%u)_%u", m, n, s);
        fprintf (file, " %-15s = %13E\n", name, B[s][mni]);
      }

  fclose (file);

  return 0;
}

unsigned writeKnotsAndCoeffs (char *output, char *prm, unsigned zone,
                              unsigned nnt, double T[nnt],
                              unsigned ns, double C[ns])
{
  /* Writes B-splines knots and coefficients */
  char name[MAX_STRING_LEN];

  FILE *file;

  if (sprintf (name, "%s/%s_Z%u_KC.dat", output, prm, zone) < 14) return 1;

  file = fopen (name, "w");

  if (file == NULL) return 2;

  fprintf (file, "Knots\n");

  for (unsigned i = 0; i < nnt; i++)

    fprintf (file, "%.8lf\n", T[i]);

  fprintf (file, "Coefficients\n");

  for (unsigned i = 0; i < ns; i++)

    fprintf (file, "%.6lf\n", C[i]);

  fclose (file);

  return 0;
}

unsigned checkTopoIO (unsigned rvalue)
{
  /* Checks IO of the topological files */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Unable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "Error opening topological file...\n");
    break;

    case 3:
      fprintf (stderr, "Error reading topological file...\n");
    break;
  }

  return rvalue;
}

unsigned checkTopoAndModelIO (unsigned rvalue)
{
  /* Checks IO of the topological and model files */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Unable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "Error opening topological file...\n");
    break;

    case 3:
      fprintf (stderr, "Error reading topological file...\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not open model file!\n");
    break;

    case 5:
      fprintf (stderr, "Error: could not read model file!\n");
    break;
  }

  return rvalue;
}

unsigned checkModelIO (unsigned rvalue)
{
  /* Checks IO of the model files */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Unable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "Error opening file...\n");
    break;

    case 3:
      fprintf (stderr, "Error writing model file...\n");
    break;
  }

  return rvalue;
}

unsigned checkBSplinesHeaderIO (unsigned rvalue)
{
  /* Checks IO of the B-splines file header */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file header!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read number of knots!\n");
    break;

    case 5:
      fprintf (stderr, "Error: could not read b-splines degree!\n");
    break;

    case 6:
      fprintf (stderr, "Error: could not read minimum radius!\n");
    break;

    case 7:
      fprintf (stderr, "Error: could not read maximum radius!\n");
    break;
  }

  return rvalue;
}

unsigned checkKnotsIO (unsigned rvalue)
{
  /* Checks IO of the B-splines file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read header information!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read knots positions!\n");
    break;

    case 5:
      fprintf (stderr, "Error: could not read minimum radius!\n");
    break;

    case 6:
      fprintf (stderr, "Error: could not read maximum radius!\n");
    break;
  }

  return rvalue;
}

unsigned checkMeanModelHeaderIO (unsigned rvalue)
{
  /* Checks IO of the mean model file header */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file header!\n");
    break;

    case 4:
      fprintf (stderr, "Error: could not read number of points!\n");
    break;
  }

  return rvalue;
}

unsigned checkMeanModelIO (unsigned rvalue)
{
  /* Checks IO of the mean model file */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not open file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file!\n");
    break;
  }

  return rvalue;
}

unsigned checkCoefficientsIO (unsigned rvalue)
{
  /* Checks IO of the model files */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Unable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "Error opening file...\n");
    break;
  }

  return rvalue;
}

