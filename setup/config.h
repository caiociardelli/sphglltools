#include <stdbool.h>

#ifndef SETUP_H
#define SETUP_H
static const unsigned NC = 24; /* Total number of slices of the input mesh */

static const unsigned NEL = 8896;   /* Total number of spectral elements on each slice of the input GLL model */
static const unsigned NG  = 592913; /* Total number of points on each slice of the input GLL model */

static const unsigned NEL1 = 8896;   /* Total number of spectral elements on each slice of the input GLL model */
static const unsigned NG1  = 592913; /* Total number of points on each slice of the input GLL model */
static const unsigned NEL2 = 8896;   /* Total number of spectral elements on each slice of the output GLL model */
static const unsigned NG2  = 592913; /* Total number of points on each slice of the output GLL model */

static const bool SEPARATE_ZONES = false; /* Smooth upper mantle, tranzition zone and lower mantle separately
                                             (preserve discontinuities) */

#define VP    0 /* Set to 1 to enable the Vp parameter */
#define VS    0 /* Set to 1 to enable the Vs parameter */
#define RHO   0 /* Set to 1 to enable the Rho parameter */
#define VPV   0 /* Set to 1 to enable the Vpv parameter */
#define VPH   0 /* Set to 1 to enable the Vph parameter */
#define VSV   1 /* Set to 1 to enable the Vsv parameter */
#define VSH   0 /* Set to 1 to enable the Vsh parameter */
#define ETA   0 /* Set to 1 to enable the Eta parameter */
#define QMU   0 /* Set to 1 to enable the Qmu parameter */
#define GCP   0 /* Set to 1 to enable the Gc_prime parameter */
#define GSP   0 /* Set to 1 to enable the Gs_prime parameter */
#define MU0   0 /* Set to 1 to enable the mu0 parameter */
#define APK   0 /* Set to 1 to enable the Alpha Kernel parameter */
#define BTK   0 /* Set to 1 to enable the Beta Kernel parameter */
#define RHK   0 /* Set to 1 to enable the Rho Kernel parameter */
#define BCK   0 /* Set to 1 to enable the Bulk_c Kernel parameter */
#define BBK   0 /* Set to 1 to enable the Bulk_beta Kernel parameter */
#define BVK   0 /* Set to 1 to enable the Bulk_betav Kernel parameter */
#define BHK   0 /* Set to 1 to enable the Bulk_betah Kernel parameter */
#define ETK   0 /* Set to 1 to enable the Eta Kernel parameter */
#define GCK   0 /* Set to 1 to enable the Gc_prime Kernel parameter */
#define GSK   0 /* Set to 1 to enable the Gs_prime Kernel parameter */
#define HSK   0 /* Set to 1 to enable the Hessian Kernel parameter */
#endif
