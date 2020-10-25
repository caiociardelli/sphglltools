#include <stdbool.h>

/* Do not change anything in this file unless you know very well what you are doing! */
#ifndef CONSTANTS_H
#define CONSTANTS_H
static const unsigned MAX_STRING_LEN = 200; /* Buffer size for file names */

static const unsigned NX = 5; /* Number of spectral elements in X direction */
static const unsigned NY = 5; /* Number of spectral elements in Y direction */
static const unsigned NZ = 5; /* Number of spectral elements in Z direction */

static const unsigned NN    = 27; /* Number of nodes on each spectral element */
static const unsigned NPe   = 5;  /* Number of points on each edge of a spectral element */
static const unsigned NNe   = 3;  /* Number of nodes on each edge of a spectral element */
static const unsigned NPf   = 25; /* Number of points on each face of a spectral element */
static const unsigned NNf   = 9;  /* Number of nodes on each face of a spectral element */
static const unsigned NCF_s = 6;  /* Number of coefficients in the simple surface */
static const unsigned NCF_p = 3;  /* Number of coefficients in the plane surface */
static const unsigned SSM   = 11; /* Size of matrix for simple surface */

static const long double ANGLE_TOLERANCE = 1E-10; /* Water level to convert from Cartesian to spherical coordinates */
static const long double TOLERANCE       = 1E-15; /* Tolerance to detect a ill-conditioned matrix */
static const long double TRUNCATE        = 3.L;   /* Number of sigmas to be used as search radius */
static const long double WATER_LEVEL     = 1e-3;  /* Water level to prevent divisions by zero in the smoothing routine */
static const long double WATER_LEVEL_1   = 1E-15; /* Water level to prevent divisions by zero */
static const long double WATER_LEVEL_2   = 1E-8;  /* Water level to determine the best face orientation */
static const long double MAX_NORM_1      = 1E8;   /* Maximum allowed norm for the coefficients */
static const long double MAX_NORM_2      = 5E-3;  /* Maximum allowed norm for the residuals */
static const long double MAX_CURVATURE   = 1E2;   /* Maximum allowed mean surface curvature */
static const long double MAX_SGM_RATIO   = 1E-1;  /* Maximum uncertainty/spectral element size ratio */

static const long double MIN_XI    = -1.05L; /* Minimum value for the xi coordinate */
static const long double MAX_XI    =  1.05L; /* Maximum value for the xi coordinate */
static const long double MIN_ETA   = -1.05L; /* Minimum value for the eta coordinate */
static const long double MAX_ETA   =  1.05L; /* Maximum value for the eta coordinate */
static const long double MIN_GAMMA = -1.05L; /* Minimum value for the gamma coordinate */
static const long double MAX_GAMMA =  1.05L; /* Maximum value for the gamma coordinate */

static const unsigned MAX_HITR  = 10; /* Max number of iterations to find the GLL nodes and weights */
static const unsigned MAX_NITER = 10; /* Max number of iterations to find the interpolation point */

static const long double EPSILON = 1E-12; /* Tolerance to convergence of the Newton's method */

static const long double GB_TOLERANCE     = 5E-3; /* Tolerance at the boundaries to read the input mesh */
static const long double BOUNDARY_RATIO_1 = 1E-1; /* Rough tolerance ratio at boundaries of the spectral elements */
static const long double BOUNDARY_RATIO_2 = 5E-3; /* Tolerance ratio at boundaries of the spectral elements */
static const long double BOUNDARY_RATIO_3 = 1E-1; /* Tolerance ratio at boundaries of the faces and edges */
static const long double BOUNDARY_RATIO_4 = 1E-3; /* Tolerance ratio around the radius of the depth slice */

static const long double D650_R          = 0.897975L; /* Normalized radius for 650 km depth */
static const long double D410_R          = 0.935646L; /* Normalized radius for 410 km depth */
static const long double CMB_R           = 0.546225L; /* Normalized radius of the core-mantle boundary */
static const long double MOHO_R          = 0.987443L; /* Normalized minimum radius of the crust-mantle boundary */
static const long double MIN_SURFACE_R   = 0.998744L; /* Normalized radius corresponding to the maximum surface depth */
static const long double MAX_SURFACE_R   = 1.000942L; /* Normalized radius corresponding to the maximum surface altitude */
static const long double MIN_STC_UPPER_R = 0.993722L; /* Normalized radius corresponding to the bottom of the stretched upper
                                                         crust */
static const long double MAX_STC_UPPER_R = 1.000000L; /* Normalized radius corresponding to the top of the stretched upper
                                                         crust */

static const long double MAX_PERTURBATION = 2E-1; /* Maximum allowed perturbation between the Lagrange and the trilinear
                                                     interpolations for the upper crust */

static const unsigned NP  = 10000; /* Number of interpolation points in the smoothing profile */
static const unsigned NPT = 100001; /* Number of interpolation points in for the basis arrays (must be odd!) */

static const unsigned ROUND = 7; /* Round mesh coordinates to this number of decimal places */

static const long double PI = 3.14159265358979323846L; /* Value of Pi */

static const long double TO_DEGREE  = 180.L / 3.14159265358979323846L; /* Constant to convert from radians to degrees */
static const long double TO_RADIANS = 3.14159265358979323846L / 180.L; /* Constant to convert from degrees to radians */

static const long double EARTH_R = 6371.L; /* Earth radius */

                                      /* K1 and K2 are constants used to prevent underflow
                                         in the Associated Legendre funtions */
static const long double K1 = 1E300;
static const long double K2 = 1E-150;

#define M_NX  5  /* Macro for NX (to be used inside structures) */
#define M_NY  5  /* Macro for NY (to be used inside structures) */
#define M_NZ  5  /* Macro for NZ (to be used inside structures) */

#define M_NPf 25 /* Macro for NPf (to be used inside structures) */
#endif
