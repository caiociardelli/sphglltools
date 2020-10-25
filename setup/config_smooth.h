#include <stdbool.h>

#ifndef SETUP_SMOOTH_H
#define SETUP_SMOOTH_H
static const unsigned NC = 24; /* Total number of slices of the input mesh */

static const unsigned NEL = 8896;   /* Total number of spectral elements on each slice of the mesh*/
static const unsigned NG  = 592913; /* Total number of points on each slice of the mesh */

static const bool SEPARATE_ZONES = false; /* Smooth upper mantle, tranzition zone and lower mantle separately
                                             (preserve discontinuities) */

#define VP    1 /* Set to 1 to enable smoothing of the Vp parameter */
#define VS    0 /* Set to 1 to enable smoothing of the Vs parameter */
#define RHO   0 /* Set to 1 to enable smoothing of the Rho parameter */
#define VPV   0 /* Set to 1 to enable smoothing of the Vpv parameter */
#define VPH   0 /* Set to 1 to enable smoothing of the Vph parameter */
#define VSV   0 /* Set to 1 to enable smoothing of the Vsv parameter */
#define VSH   0 /* Set to 1 to enable smoothing of the Vsh parameter */
#define ETA   0 /* Set to 1 to enable smoothing of the Eta parameter */
#define QMU   0 /* Set to 1 to enable smoothing of the Qmu parameter */
#define APK   0 /* Set to 1 to enable smoothing of the Alpha Kernel parameter */
#define BTK   0 /* Set to 1 to enable smoothing of the Beta Kernel parameter */
#define RHK   0 /* Set to 1 to enable smoothing of the Rho Kernel parameter */
#define BCK   0 /* Set to 1 to enable smoothing of the Bulk_c Kernel parameter */
#define BBK   0 /* Set to 1 to enable smoothing of the Bulk_beta Kernel parameter */
#define BVK   0 /* Set to 1 to enable smoothing of the Bulk_betav Kernel parameter */
#define BHK   0 /* Set to 1 to enable smoothing of the Bulk_betah Kernel parameter */
#define ETK   0 /* Set to 1 to enable smoothing of the Eta Kernel parameter */
#define HSK   0 /* Set to 1 to enable smoothing of the Hessian Kernel parameter */
#endif
