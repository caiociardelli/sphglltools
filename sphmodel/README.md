# SphModel
-----------------------

`SphModel` is a subpackage of `SphGLLTools`, a toolbox for visualization, processing and spherical harmonics analysis of tomographic models defined on GLL meshes.

Author: Caio Ciardelli

If you use SphModel, please, cite the following paper:

Ciardelli, C., Bozdağ, E. and Peter, D., 2020. SphGLLTools: A set of routines for visualization, processing, sharing, and spherical harmonics analysis of tomographic models defined on GLL meshes. Computer & Geosciences, submitted.

## Installation
-----------------------

Instructions on how to install and use `SphModel` are available in the PDF manual of `SphGLLTools`: User_Manual.pdf

## Development
-----------------------

Development is hosted on GitHub in the [caio.ciardelli/sphglltools repository](https://github.com/caio.ciardelli/sphglltools).

# YOUR_MODEL_NAME
-----------------------

Your publication

# Plate boundaries
-----------------------

The plate boundaries were taken from http://geoscience.wisc.edu/~chuck/MORVEL/citation.html (access on June 20, 2020)

DeMets, C., Gordon, R. G., and Argus, D. F., 2010. Geologically current plate motions, Geophysical Journal International, v. 181, no. 1, p. 1-80, doi: 10.1111/j.1365-246X.2009.04491.x

# Model parametrization
-----------------------

This parametrization includes PARAMETERS (UNITIES) (e.g. Vp (km/s), Vs (km/s) and Rho (g/cm^3)). The perturbations for the above mentioned parameters as well as some derived parameters such as bulk sound speed, Vp/Vs ratios and transverse isotropy can also be computed using these routines.

The model is separated in four zones:

Zone 1: Represents the crust using a block model. The horizontal resolution is ? degrees both in latitude 
and longitude. The vertical resolution is ? km. This zone stretches from ? km above sea level to ? km depth.

Zone 2: Represents the upper mantle using a combination of spherical harmonics up to degree ? and ? cubic 
B-splines. This zone stretches from ? km to ? km depth.

Zone 3: Represents the transition zone using a combination of spherical harmonics up to degree ? and ? cubic 
B-splines. This zone stretches from ? km to ? km depth.

Zone 4: Represents the lower mantle using a combination of spherical harmonics up to degree ? and ? cubic 
b-splines. This zone stretches from ? km to ? km depth.

# Usage
-----------------------

All routines have a help menu that shows up wherever you run them with no or with a wrong number of command-line parameters. The
same menu is also at the beginning of each source code. For more details, please, refer to the manual of `SphGLLTools`.

#Contact
-----------------------

Is you have any questions, suggestions and bug reports you can email *caio.ciardelli@gmail.com*

