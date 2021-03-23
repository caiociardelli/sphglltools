#ifndef SETUP_H
#define SETUP_H
static const unsigned NC = 24; /* Total number of slices of the input mesh */

static const unsigned NEL = 8896;   /* Total number of spectral elements on each slice of the input GLL model */
static const unsigned NG  = 592913; /* Total number of points on each slice of the input GLL model */

static const unsigned NEL1 = 8896;   /* Total number of spectral elements on each slice of the input GLL model */
static const unsigned NG1  = 592913; /* Total number of points on each slice of the input GLL model */
static const unsigned NEL2 = 8896;   /* Total number of spectral elements on each slice of the output GLL model */
static const unsigned NG2  = 592913; /* Total number of points on each slice of the output GLL model */

#define VP   1 /* Set to 1 to enable the interpolation of the Vp parameter */
#define VS   1 /* Set to 1 to enable the interpolation of the Vs parameter */
#define RHO  1 /* Set to 1 to enable the interpolation of the Rho parameter */
#define VPV  1 /* Set to 1 to enable the interpolation of the Vpv parameter */
#define VPH  1 /* Set to 1 to enable the interpolation of the Vph parameter */
#define VSV  1 /* Set to 1 to enable the interpolation of the Vsv parameter */
#define VSH  1 /* Set to 1 to enable the interpolation of the Vsh parameter */
#define ETA  1 /* Set to 1 to enable the interpolation of the Eta parameter */
#define QMU  0 /* Set to 1 to enable the interpolation of the Qmu parameter */
#endif
