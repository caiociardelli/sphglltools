#!/usr/bin/env bash

# SphGLLTools

# Author: Caio Ciardelli, University of São Paulo, October 2020

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#-----------------------------------------------------------------------------------------------

# PLOT_PROFILE

# USAGE
#   ./utils/plot_profile.bash PARAMETER

# EXAMPLE
#   ./utils/plot_profile.bash vpv

# COMMAND LINE ARGUMENTS
#   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)

# DESCRIPTION
#   Plots the profile for the grid files created by GLL2PF or by CREATE_EXTRA_PROFILE.

#-----------------------------------------------------------------------------------------------

title="S362ANI"

if [ $1 = "vp" ]; then
  label="V@-P@- (km s@+-1@+)"
elif [ $1 = "vpv" ]; then
  label="V@-PV@- (km s@+-1@+)"
elif [ $1 = "vph" ]; then
  label="V@-PH@- (km s@+-1@+)"
elif [ $1 = "vs" ]; then
  label="V@-S@- (km s@+-1@+)"
elif [ $1 = "vsv" ]; then
  label="V@-SV@- (km s@+-1@+)"
elif [ $1 = "vsh" ]; then
  label="V@-SH@- (km s@+-1@+)"
elif [ $1 = "eta" ]; then
  label="@~\150@~"
elif [ $1 = "rho" ]; then
  label="@~\162@~ (g cm@+-3@+)"
elif [ $1 = "qmu" ]; then
  label="Q@-@~\155@~@-"
elif [ $1 = "cb" ]; then
  label="C@-Bulk@- (km s@+-1@+)"
elif [ $1 = "ti" ]; then
  label="T@-i@- (%)"
elif [ $1 = "vpvs" ]; then
  label="V@-P@-/V@-S@-"
elif [ $1 = "dvp" ]; then
  label="dlnV@-P@- (%)"
elif [ $1 = "dvpv" ]; then
  label="dlnV@-PV@- (%)"
elif [ $1 = "dvph" ]; then
  label="dlnV@-PH@- (%)"
elif [ $1 = "dvs" ]; then
  label="dlnV@-S@- (%)"
elif [ $1 = "dvsv" ]; then
  label="dlnV@-SV@- (%)"
elif [ $1 = "dvsh" ]; then
  label="dlnV@-SH@- (%)"
elif [ $1 = "deta" ]; then
  label="dln@~\150@~ (%)"
elif [ $1 = "drho" ]; then
  label="dln@~\162@~ (%)"
elif [ $1 = "dqmu" ]; then
  label="dQ@-@~\155@~@- (%)"
elif [ $1 = "dcb" ]; then
  label="dlnC@-Bulk@- (%)"
else
  echo "Error: $1 parameter is not available!"
  exit 1
fi

help()
{
  echo " PLOT_PROFILE"
  echo ""
  echo " USAGE"
  echo "   ./utils/plot_profile.bash PARAMETER"
  echo ""
  echo " EXAMPLE"
  echo "   ./utils/plot_profile.bash vp"
  echo ""
  echo " COMMAND LINE ARGUMENTS"
  echo "   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)"
  echo ""
  echo " DESCRIPTION"
  echo "   Plots the profile for the grid files created by GLL2PF or by CREATE_EXTRA_PROFILE."
  echo ""
}

if [ "$#" -ne 1 ]; then
  help
  exit 1
fi

filename="$1"_PF.dat
grdname="$1"_PF.grd
output="$1"_PF

ndep=$(head -n1 $filename | cut -f2 -d' ')

line=$(sed '2q;d' $filename)

lat=$(echo "$line" | awk '{print $3}')
lon=$(echo "$line" | awk '{print $4}')

string=$(gmt info $filename)

dstring=$(echo "$string" | cut -f2 | awk -F '[</>]' '{print $2" "$3}')
vstring=$(echo "$string" | cut -f3 | awk -F '[</>]' '{print $2" "$3}')

dmin=$(echo "$dstring" | awk '{print $1}')
dmax=$(echo "$dstring" | awk '{print $2}')
vmin=$(echo "$vstring" | awk '{print $1}')
vmax=$(echo "$vstring" | awk '{print $2}')

awk -v min=$vmin -v max=$vmax 'BEGIN {printf "Min = %E Max = %E\n", min, max}'
echo 'Creating figure...'

vrange=$(echo "$vstring" | awk '{print $2 - $1}')

vmin=$(echo "$vmin - 0.1 * $vrange" | bc -l)
vmax=$(echo "$vmax + 0.1 * $vrange" | bc -l)

gmt begin $output pdf
  gmt set FONT_TITLE 20p,100

  gmt basemap -R$vmin/$vmax/$dmin/$dmax -JX5c/-7.5c -BWSne+t$title -Bxafg+l"$label" -Bya500fg500+l"depth (km)"

  gmt plot "$1"_PF.dat -: -W1.5p,blue -h2

  gmt inset begin -DJTR+w5c -M0.1c
    gmt coast -JG$lon/$lat/2.4c -Rg -Bafg -Dc -A10000 -Glightgray -Wthinnest

    gmt plot extra/plate_boundaries.dat -JG$lon/$lat/2.4c -: -W0.6p,brown

    echo "$lon $lat" | gmt plot -JG$lon/$lat/2.4c -Sc0.2c -Ggreen -W0.03c,black
  gmt inset end
gmt end

echo 'Figure saved!'
