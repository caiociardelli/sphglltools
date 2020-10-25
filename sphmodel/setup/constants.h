/* Do not change anything in this file unless you know very well what you are doing! */
#ifndef CONSTANTS_H
#define CONSTANTS_H
static const unsigned MAX_STRING_LEN = 200; /* Buffer size for file names */

static const unsigned NPT = 100001; /* Number of interpolation points in the radial basis (must be odd!) */

static const double ANGLE_TOLERANCE = 1e-10; /* Water level to convert from Cartesian to spherical coordinates */
static const double WATER_LEVEL     = 1e-15; /* Water level to prevent divisions by zero */

static const double CMB_R  = 0.54622508; /* Normalized core-mantle boundary radius */
static const double MOHO_R = 0.98744310; /* Normalized Moho radius */
static const double TOP_R  = 1.00078481; /* Normalized radius corresponding to the maximum surface altitude */

static const double PI = 3.14159265358979323; /* Value of Pi */

static const double TO_DEGREE  = 180.0 / 3.14159265358979323; /* Constant to convert from radians to degrees */
static const double TO_RADIANS = 3.14159265358979323 / 180.0; /* Constant to convert from degrees to radians */

static const double EARTH_R = 6371.0; /* Earth radius */

                                      /* K1 and K2 are constants used to prevent underflow
                                         in the Associated Legendre funtions */
static const long double K1 = 1E300;
static const long double K2 = 1E-150;
#endif
