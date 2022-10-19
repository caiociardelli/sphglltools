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
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "structs_smooth.h"
#include "exmath.h"
#include "metrics_smooth.h"
#include "coordinates.h"

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

static void insertElNode (struct Point Ti[NX][NY][NZ], struct ElNode *lle,
                          unsigned *nel)
{
  /* Adds a new cell to the linked list of spectral elements */
  struct ElNode *p = lle, *new = malloc (sizeof (struct ElNode));

  for (unsigned i = 0; i < NX; i++)

    for (unsigned j = 0; j < NY; j++)

      for (unsigned k = 0; k < NZ; k++)
      {
        new->Ti[i][j][k].x = Ti[i][j][k].x;
        new->Ti[i][j][k].y = Ti[i][j][k].y;
        new->Ti[i][j][k].z = Ti[i][j][k].z;

        new->Ti[i][j][k].r     = Ti[i][j][k].r;
        new->Ti[i][j][k].theta = Ti[i][j][k].theta;
        new->Ti[i][j][k].phi   = Ti[i][j][k].phi;
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

        #if GCP
        new->M[i][j][k].gcp = M[i][j][k].gcp;
        #endif

        #if GSP
        new->M[i][j][k].gsp = M[i][j][k].gsp;
        #endif

        #if MU0
        new->M[i][j][k].mu0 = M[i][j][k].mu0;
        #endif

        #if APK
        new->M[i][j][k].apk = M[i][j][k].apk;
        #endif

        #if BTK
        new->M[i][j][k].btk = M[i][j][k].btk;
        #endif

        #if RHK
        new->M[i][j][k].rhk = M[i][j][k].rhk;
        #endif

        #if BCK
        new->M[i][j][k].bck = M[i][j][k].bck;
        #endif

        #if BBK
        new->M[i][j][k].bbk = M[i][j][k].bbk;
        #endif

        #if BVK
        new->M[i][j][k].bvk = M[i][j][k].bvk;
        #endif

        #if BHK
        new->M[i][j][k].bhk = M[i][j][k].bhk;
        #endif

        #if ETK
        new->M[i][j][k].etk = M[i][j][k].etk;
        #endif

        #if GCK
        new->M[i][j][k].gck = M[i][j][k].gck;
        #endif

        #if GSK
        new->M[i][j][k].gsk = M[i][j][k].gsk;
        #endif

        #if HSK
        new->M[i][j][k].hsk = M[i][j][k].hsk;
        #endif
      }

  new->next = p->next;
  p->next = new;
}

static void deleteElNodes (struct ElNode *lle)
{
  /* Deletes all cells of the linked list of spectral elements */
  struct ElNode *q;

  for (struct ElNode *p = lle; p != NULL; p = q)
  {
    q = p->next;
    free (p);
  }

  q = NULL;
}

static void deletePmNodes (struct PmNode *llm)
{
  /* Deletes all cells of the linked list of model parameters */
  struct PmNode *q;

  for (struct PmNode *p = llm; p != NULL; p = q)
  {
    q = p->next;
    free (p);
  }

  q = NULL;
}

unsigned readInputMeshAndModel (int ic, char *prm1, char *prm2,
                                unsigned n, struct Profile p[n],
                                struct Boundaries *gb, struct ElNode *lle,
                                struct PmNode *llm, unsigned *nel)
{
  /* Reads the input mesh and model */
  unsigned N = NEL * NX * NY * NZ;

  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s/proc%06d_reg1_solver_data.bin", prm1, ic) < 33) return 1;

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

  unsigned ibool[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ibool, sizeof (unsigned), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned idoubling[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (idoubling, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ispec_is_tiso[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ispec_is_tiso, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xix[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xix, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiy[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiy, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);

  #if VP
  if (sprintf (filename, "%s/proc%06d_reg1_vp.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VS
  if (sprintf (filename, "%s/proc%06d_reg1_vs.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vs[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vs, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if RHO
  if (sprintf (filename, "%s/proc%06d_reg1_rho.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float rho[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (rho, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPV
  if (sprintf (filename, "%s/proc%06d_reg1_vpv.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vpv[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vpv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPH
  if (sprintf (filename, "%s/proc%06d_reg1_vph.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vph[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vph, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSV
  if (sprintf (filename, "%s/proc%06d_reg1_vsv.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsv[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSH
  if (sprintf (filename, "%s/proc%06d_reg1_vsh.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsh[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsh, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if ETA
  if (sprintf (filename, "%s/proc%06d_reg1_eta.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float eta[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (eta, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if QMU
  if (sprintf (filename, "%s/proc%06d_reg1_qmu.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float qmu[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (qmu, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCP
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gcp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gcp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCP
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gsp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gsp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if MU0
  if (sprintf (filename, "%s/proc%06d_reg1_mu0.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float mu0[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (mu0, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if APK
  if (sprintf (filename, "%s/proc%06d_reg1_alpha_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float apk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (apk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BTK
  if (sprintf (filename, "%s/proc%06d_reg1_beta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float btk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (btk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if RHK
  if (sprintf (filename, "%s/proc%06d_reg1_rho_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float rhk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (rhk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BCK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_c_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bck[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bck, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BBK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_beta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bbk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bbk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BVK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betav_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bvk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bvk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BHK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betah_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bhk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bhk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if ETK
  if (sprintf (filename, "%s/proc%06d_reg1_eta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float etk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (etk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCK
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gck[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gck, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCK
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gsk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gsk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if HSK
  if (sprintf (filename, "%s/proc%06d_reg1_hess_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float hsk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (hsk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  struct Point Ti[NX][NY][NZ];
  struct Parameters M[NX][NY][NZ];

  double tlc_h = 0, tlc_v = 0;

  for (unsigned i = 0; i < n; i++)
  {
    if (tlc_h < p[i].search_h) tlc_h = p[i].search_h;
    if (tlc_v < p[i].search_v) tlc_v = p[i].search_v;
  }

  tlc_h += GB_TOLERANCE;
  tlc_v += GB_TOLERANCE;

  for (unsigned el = 0; el < NEL; el++)
  {
    bool elisin = false;

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          unsigned ig = ibool[el][k][j][i] - 1;

          double x = xf[ig];
          double y = yf[ig];
          double z = zf[ig];

          Ti[i][j][k].x = x;
          Ti[i][j][k].y = y;
          Ti[i][j][k].z = z;

          double r, theta, phi;

          xYZ2RThetaPhi (x, y, z, &r, &theta, &phi);

          Ti[i][j][k].r     = r;
          Ti[i][j][k].theta = theta;
          Ti[i][j][k].phi   = phi;

          #if VP
          M[i][j][k].vp  = vp[el][k][j][i];
          #endif

          #if VS
          M[i][j][k].vs  = vs[el][k][j][i];
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

          #if APK
          M[i][j][k].apk = apk[el][k][j][i];
          #endif

          #if BTK
          M[i][j][k].btk = btk[el][k][j][i];
          #endif

          #if RHK
          M[i][j][k].rhk = rhk[el][k][j][i];
          #endif

          #if BCK
          M[i][j][k].bck = bck[el][k][j][i];
          #endif

          #if BBK
          M[i][j][k].bbk = bbk[el][k][j][i];
          #endif

          #if BVK
          M[i][j][k].bvk = bvk[el][k][j][i];
          #endif

          #if BHK
          M[i][j][k].bhk = bhk[el][k][j][i];
          #endif

          #if ETK
          M[i][j][k].etk = etk[el][k][j][i];
          #endif

          #if HSK
          M[i][j][k].hsk = hsk[el][k][j][i];
          #endif

          if (r      >= gb->rmin - tlc_v && r     <= gb->rmax + tlc_v &&
              theta  >= gb->tmin - tlc_h && theta <= gb->tmax + tlc_h &&
              ((phi  >= gb->pmin - tlc_h && phi   <= gb->pmax + tlc_h) ||
              (theta <= tlc_h || theta >= PI - tlc_h)))

              elisin = true;
        }

    if (elisin)
    {
      insertElNode (Ti, lle, nel);
      insertPmNode (M, llm);
    }
  }

  return 0;
}

unsigned readOutputMeshAndModel (int ic, char *prm1, char *prm2,
                                 struct Boundaries *gb,
                                 struct Point To[NEL][NX][NY][NZ],
                                 struct Parameters Mo[NEL][NX][NY][NZ])
{
  /* Reads the output mesh and model */
  unsigned N = NEL * NX * NY * NZ;

  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "%s/proc%06d_reg1_solver_data.bin", prm1, ic) < 33) return 1;

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

  unsigned ibool[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ibool, sizeof (unsigned), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned idoubling[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (idoubling, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  unsigned ispec_is_tiso[NEL];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (ispec_is_tiso, sizeof (unsigned), NEL, file) != NEL) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xix[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xix, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiy[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiy, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float xiz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (xiz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float etaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (etaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammax[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammax, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammay[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammay, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  float gammaz[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;
  if (fread (gammaz, sizeof (float), N, file) != N) return 3;
  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);

  #if VP
  if (sprintf (filename, "%s/proc%06d_reg1_vp.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VS
  if (sprintf (filename, "%s/proc%06d_reg1_vs.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vs[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vs, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if RHO
  if (sprintf (filename, "%s/proc%06d_reg1_rho.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float rho[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (rho, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPV
  if (sprintf (filename, "%s/proc%06d_reg1_vpv.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vpv[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vpv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VPH
  if (sprintf (filename, "%s/proc%06d_reg1_vph.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vph[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vph, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSV
  if (sprintf (filename, "%s/proc%06d_reg1_vsv.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsv[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsv, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if VSH
  if (sprintf (filename, "%s/proc%06d_reg1_vsh.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float vsh[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (vsh, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if ETA
  if (sprintf (filename, "%s/proc%06d_reg1_eta.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float eta[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (eta, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if QMU
  if (sprintf (filename, "%s/proc%06d_reg1_qmu.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float qmu[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (qmu, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCP
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gcp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gcp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GSP
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gsp[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gsp, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if MU0
  if (sprintf (filename, "%s/proc%06d_reg1_mu0.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float mu0[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (mu0, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if APK
  if (sprintf (filename, "%s/proc%06d_reg1_alpha_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float apk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (apk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BTK
  if (sprintf (filename, "%s/proc%06d_reg1_beta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float btk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (btk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if RHK
  if (sprintf (filename, "%s/proc%06d_reg1_rho_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float rhk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (rhk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BCK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_c_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bck[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bck, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BBK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_beta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bbk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bbk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BVK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betav_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bvk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bvk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if BHK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betah_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float bhk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (bhk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if ETK
  if (sprintf (filename, "%s/proc%06d_reg1_eta_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float etk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (etk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GCK
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gck[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gck, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if GSK
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float gsk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (gsk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  #if HSK
  if (sprintf (filename, "%s/proc%06d_reg1_hess_kernel.bin", prm2, ic) < 24) return 1;

  file = fopen (filename, "rb");

  if (file == NULL) return 4;

  float hsk[NEL][NZ][NY][NX];

  if (fread (&junk, sizeof (unsigned), 1, file) != 1) return 5;
  if (fread (hsk, sizeof (float), N, file) != N) return 5;

  fclose (file);
  #endif

  gb->rmin =  INFINITY;
  gb->rmax = -INFINITY;
  gb->tmin =  INFINITY;
  gb->tmax = -INFINITY;
  gb->pmin =  INFINITY;
  gb->pmax = -INFINITY;

  for (unsigned el = 0; el < NEL; el++)

    for (unsigned i = 0; i < NX; i++)

      for (unsigned j = 0; j < NY; j++)

        for (unsigned k = 0; k < NZ; k++)
        {
          unsigned ig = ibool[el][k][j][i] - 1;

          double x = (double) xf[ig];
          double y = (double) yf[ig];
          double z = (double) zf[ig];

          To[el][i][j][k].x = x;
          To[el][i][j][k].y = y;
          To[el][i][j][k].z = z;

          double r, theta, phi;

          xYZ2RThetaPhi (x, y, z, &r, &theta, &phi);

          To[el][i][j][k].r     = r;
          To[el][i][j][k].theta = theta;
          To[el][i][j][k].phi   = phi;

          if (r < gb->rmin)     gb->rmin = r;
          if (r > gb->rmax)     gb->rmax = r;
          if (theta < gb->tmin) gb->tmin = theta;
          if (theta > gb->tmax) gb->tmax = theta;
          if (phi < gb->pmin)   gb->pmin = phi;
          if (phi > gb->pmax)   gb->pmax = phi;

          #if VP
          Mo[el][i][j][k].vp  = vp[el][k][j][i];
          #endif

          #if VS
          Mo[el][i][j][k].vs  = vs[el][k][j][i];
          #endif

          #if RHO
          Mo[el][i][j][k].rho = rho[el][k][j][i];
          #endif

          #if VPV
          Mo[el][i][j][k].vpv = vpv[el][k][j][i];
          #endif

          #if VPH
          Mo[el][i][j][k].vph = vph[el][k][j][i];
          #endif

          #if VSV
          Mo[el][i][j][k].vsv = vsv[el][k][j][i];
          #endif

          #if VSH
          Mo[el][i][j][k].vsh = vsh[el][k][j][i];
          #endif

          #if ETA
          Mo[el][i][j][k].eta = eta[el][k][j][i];
          #endif

          #if QMU
          Mo[el][i][j][k].qmu = qmu[el][k][j][i];
          #endif

          #if GCP
          Mo[el][i][j][k].gcp = gcp[el][k][j][i];
          #endif

          #if GSP
          Mo[el][i][j][k].gsp = gsp[el][k][j][i];
          #endif

          #if MU0
          Mo[el][i][j][k].mu0 = mu0[el][k][j][i];
          #endif

          #if APK
          Mo[el][i][j][k].apk = apk[el][k][j][i];
          #endif

          #if BTK
          Mo[el][i][j][k].btk = btk[el][k][j][i];
          #endif

          #if RHK
          Mo[el][i][j][k].rhk = rhk[el][k][j][i];
          #endif

          #if BCK
          Mo[el][i][j][k].bck = bck[el][k][j][i];
          #endif

          #if BBK
          Mo[el][i][j][k].bbk = bbk[el][k][j][i];
          #endif

          #if BVK
          Mo[el][i][j][k].bvk = bvk[el][k][j][i];
          #endif

          #if BHK
          Mo[el][i][j][k].bhk = bhk[el][k][j][i];
          #endif

          #if ETK
          Mo[el][i][j][k].etk = etk[el][k][j][i];
          #endif

          #if GCK
          Mo[el][i][j][k].gck = gck[el][k][j][i];
          #endif

          #if GSK
          Mo[el][i][j][k].gsk = gsk[el][k][j][i];
          #endif

          #if HSK
          Mo[el][i][j][k].hsk = hsk[el][k][j][i];
          #endif
        }

  return 0;
}

unsigned scanMeshAndModel (int ic, char *prm1, char *prm2,
                           unsigned n, struct Profile p[n],
                           struct Boundaries *gb,
                           struct ElNode *lle, struct PmNode *llm,
                           unsigned *nel)
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

      #if GCP
      fprintf (stderr, "Reading model file 'proc%06d_reg1_Gc_prime.bin'...\n", iic);
      #endif

      #if GSP
      fprintf (stderr, "Reading model file 'proc%06d_reg1_Gs_prime.bin'...\n", iic);
      #endif

      #if MU0
      fprintf (stderr, "Reading model file 'proc%06d_reg1_mu0.bin'...\n", iic);
      #endif

      #if APK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_alpha_kernel.bin'...\n", iic);
      #endif

      #if BTK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_beta_kernel.bin'...\n", iic);
      #endif

      #if RHK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_rho_kernel.bin'...\n", iic);
      #endif

      #if BCK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_bulk_c_kernel.bin'...\n", iic);
      #endif

      #if BBK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_bulk_beta_kernel.bin'...\n", iic);
      #endif

      #if BVK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_bulkv_beta_kernel.bin'...\n", iic);
      #endif

      #if BHK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_bulkh_beta_kernel.bin'...\n", iic);
      #endif

      #if ETK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_eta_kernel.bin'...\n", iic);
      #endif

      #if GCK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_Gc_prime_kernel.bin'...\n", iic);
      #endif

      #if GSK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_Gs_prime_kernel.bin'...\n", iic);
      #endif

      #if HSK
      fprintf (stderr, "Reading model file 'proc%06d_reg1_hess_kernel.bin'...\n", iic);
      #endif
    }

    unsigned r = readInputMeshAndModel (iic, prm1, prm2, n, p,
                                        gb, lle, llm, nel); if (r) return r;
  }

  return 0;
}

void toArrayElAndPmNodes (unsigned nel, struct ElNode *lle,
                          struct PmNode *llm,
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
          Ti[el][i][j][k].x     = p->Ti[i][j][k].x;
          Ti[el][i][j][k].y     = p->Ti[i][j][k].y;
          Ti[el][i][j][k].z     = p->Ti[i][j][k].z;

          Ti[el][i][j][k].r     = p->Ti[i][j][k].r;
          Ti[el][i][j][k].theta = p->Ti[i][j][k].theta;
          Ti[el][i][j][k].phi   = p->Ti[i][j][k].phi;

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

          #if GCP
          M[el][i][j][k].gcp = q->M[i][j][k].gcp;
          #endif

          #if GSP
          M[el][i][j][k].gsp = q->M[i][j][k].gsp;
          #endif

          #if MU0
          M[el][i][j][k].mu0 = q->M[i][j][k].mu0;
          #endif

          #if APK
          M[el][i][j][k].apk = q->M[i][j][k].apk;
          #endif

          #if BTK
          M[el][i][j][k].btk = q->M[i][j][k].btk;
          #endif

          #if RHK
          M[el][i][j][k].rhk = q->M[i][j][k].rhk;
          #endif

          #if BCK
          M[el][i][j][k].bck = q->M[i][j][k].bck;
          #endif

          #if BBK
          M[el][i][j][k].bbk = q->M[i][j][k].bbk;
          #endif

          #if BVK
          M[el][i][j][k].bvk = q->M[i][j][k].bvk;
          #endif

          #if BHK
          M[el][i][j][k].bhk = q->M[i][j][k].bhk;
          #endif

          #if ETK
          M[el][i][j][k].etk = q->M[i][j][k].etk;
          #endif

          #if GCK
          M[el][i][j][k].gck = q->M[i][j][k].gck;
          #endif

          #if GSK
          M[el][i][j][k].gsk = q->M[i][j][k].gsk;
          #endif

          #if HSK
          M[el][i][j][k].hsk = q->M[i][j][k].hsk;
          #endif
        }

    p = p->next;
    q = q->next;
  }

  deleteElNodes (lle);
  deletePmNodes (llm);
}

unsigned readNumberOfPoints (unsigned *n)
{
  /* Reads number of points in the
     smoothing profile */
  FILE *file = fopen ("setup/parameters.cfg", "r");

  if (file == NULL) return 1;

  char line[MAX_STRING_LEN];

  if (!fgets (line, sizeof (line), file)) return 2;
  if (!fgets (line, sizeof (line), file)) return 2;
  if (!fgets (line, sizeof (line), file)) return 2;

  if (fscanf (file, "%u", n) != 1) return 2;

  return 0;
}

unsigned readSmoothingProfile (unsigned n, double *sigma,
                               struct Profile p[n])
{
  /* Reads smoothing profile */
  FILE *file = fopen ("setup/parameters.cfg", "r");

  if (file == NULL) return 1;

  char line[MAX_STRING_LEN];

  if (!fgets (line, sizeof (line), file)) return 2;
  if (!fgets (line, sizeof (line), file)) return 2;
  if (!fgets (line, sizeof (line), file)) return 2;
  if (!fgets (line, sizeof (line), file)) return 2;

  if (fscanf (file, "%lf", sigma) != 1) return 3;

  unsigned l = 0;

  double depth, sigma_h, sigma_v;

  for (unsigned i = 0; i < n; i++)
  {
    l += fscanf (file, "%lf %lf %lf", &depth, &sigma_h, &sigma_v);

    p[i].r = depth2R (depth);

    if (i > 0 && p[i].r >= p[i - 1].r) return 3;

    p[i].sigma_h  = (sigma_h + WATER_LEVEL) / EARTH_R;
    p[i].sigma_v  = (sigma_v + WATER_LEVEL) / EARTH_R;
    p[i].search_h = TRUNCATE * p[i].sigma_h;
    p[i].search_v = TRUNCATE * p[i].sigma_v;
  }

  fclose (file);

  if (l != 3 * n) return 2;

  return 0;
}

unsigned writeModel (int ic, char *prm,
                     struct Parameters M[NEL][NX][NY][NZ])
{
  /* Writes smoothed model */
  unsigned N = NEL * NX * NY * NZ;
  unsigned size = 4 * N;

  char filename[MAX_STRING_LEN];

  #if VP
  float vp[NEL][NZ][NY][NX];
  #endif

  #if VS
  float vs[NEL][NZ][NY][NX];
  #endif

  #if RHO
  float rho[NEL][NZ][NY][NX];
  #endif

  #if VPV
  float vpv[NEL][NZ][NY][NX];
  #endif

  #if VPH
  float vph[NEL][NZ][NY][NX];
  #endif

  #if VSV
  float vsv[NEL][NZ][NY][NX];
  #endif

  #if VSH
  float vsh[NEL][NZ][NY][NX];
  #endif

  #if ETA
  float eta[NEL][NZ][NY][NX];
  #endif

  #if QMU
  float qmu[NEL][NZ][NY][NX];
  #endif

  #if GCP
  float gcp[NEL][NZ][NY][NX];
  #endif

  #if GSP
  float gsp[NEL][NZ][NY][NX];
  #endif

  #if MU0
  float mu0[NEL][NZ][NY][NX];
  #endif

  #if APK
  float apk[NEL][NZ][NY][NX];
  #endif

  #if BTK
  float btk[NEL][NZ][NY][NX];
  #endif

  #if RHK
  float rhk[NEL][NZ][NY][NX];
  #endif

  #if BCK
  float bck[NEL][NZ][NY][NX];
  #endif

  #if BBK
  float bbk[NEL][NZ][NY][NX];
  #endif

  #if BVK
  float bvk[NEL][NZ][NY][NX];
  #endif

  #if BHK
  float bhk[NEL][NZ][NY][NX];
  #endif

  #if ETK
  float etk[NEL][NZ][NY][NX];
  #endif

  #if GCK
  float gck[NEL][NZ][NY][NX];
  #endif

  #if GSK
  float gsk[NEL][NZ][NY][NX];
  #endif

  #if HSK
  float hsk[NEL][NZ][NY][NX];
  #endif

  for (unsigned el = 0; el < NEL; el++)

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

          #if GCP
          gcp[el][k][j][i] = (float) M[el][i][j][k].gcp;
          #endif

          #if GSP
          gsp[el][k][j][i] = (float) M[el][i][j][k].gsp;
          #endif

          #if MU0
          mu0[el][k][j][i] = (float) M[el][i][j][k].mu0;
          #endif

          #if APK
          apk[el][k][j][i] = (float) M[el][i][j][k].apk;
          #endif

          #if BTK
          btk[el][k][j][i] = (float) M[el][i][j][k].btk;
          #endif

          #if RHK
          rhk[el][k][j][i] = (float) M[el][i][j][k].rhk;
          #endif

          #if BCK
          bck[el][k][j][i] = (float) M[el][i][j][k].bck;
          #endif

          #if BBK
          bbk[el][k][j][i] = (float) M[el][i][j][k].bbk;
          #endif

          #if BVK
          bvk[el][k][j][i] = (float) M[el][i][j][k].bvk;
          #endif

          #if BHK
          bhk[el][k][j][i] = (float) M[el][i][j][k].bhk;
          #endif

          #if ETK
          etk[el][k][j][i] = (float) M[el][i][j][k].etk;
          #endif

          #if GCK
          gck[el][k][j][i] = (float) M[el][i][j][k].gck;
          #endif

          #if GSK
          gsk[el][k][j][i] = (float) M[el][i][j][k].gsk;
          #endif

          #if HSK
          hsk[el][k][j][i] = (float) M[el][i][j][k].hsk;
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

  #if GCP
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (gcp, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if GSP
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (gsp, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if MU0
  if (sprintf (filename, "%s/proc%06d_reg1_mu0.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (mu0, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if APK
  if (sprintf (filename, "%s/proc%06d_reg1_alpha_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (apk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if BTK
  if (sprintf (filename, "%s/proc%06d_reg1_beta_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (btk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if RHK
  if (sprintf (filename, "%s/proc%06d_reg1_rho_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (rhk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if BCK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_c_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (bck, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if BBK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_beta_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (bbk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if BVK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betav_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (bvk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if BHK
  if (sprintf (filename, "%s/proc%06d_reg1_bulk_betah_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (bhk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if ETK
  if (sprintf (filename, "%s/proc%06d_reg1_eta_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (etk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if GCK
  if (sprintf (filename, "%s/proc%06d_reg1_Gc_prime_kernel.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (gck, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if GSK
  if (sprintf (filename, "%s/proc%06d_reg1_Gs_prime_kernel.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (gsk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  #if HSK
  if (sprintf (filename, "%s/proc%06d_reg1_hess_kernel_smooth.bin", prm, ic) < 24) return 1;

  file = fopen (filename, "wb");

  if (file == NULL) return 2;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;
  if (fwrite (hsk, sizeof (float), N, file) != N) return 3;
  if (fwrite (&size, sizeof (unsigned), 1, file) != 1) return 3;

  fclose (file);
  #endif

  return 0;
}

unsigned checkProfileIO (unsigned rvalue)
{
  /* Checks IO of the smoothing profile */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "\nError opening the smoothing profile...\n");
    break;

    case 2:
      fprintf (stderr, "\nError reading the smoothing profile...\n");
    break;

    case 3:
      fprintf (stderr, "\nError in the smoothing profile format: "
                       "each new depth must be larger than the previous one!\n");
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
      fprintf (stderr, "\nUnable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "\nError opening topological file...\n");
    break;

    case 3:
      fprintf (stderr, "\nError reading topological file...\n");
    break;

    case 4:
      fprintf (stderr, "\nError: could not open model file!\n");
    break;

    case 5:
      fprintf (stderr, "\nError: could not read model file!\n");
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
      fprintf (stderr, "\nUnable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "\nError opening topological file...\n");
    break;

    case 3:
      fprintf (stderr, "\nError writing model file...\n");
    break;
  }

  return rvalue;
}

