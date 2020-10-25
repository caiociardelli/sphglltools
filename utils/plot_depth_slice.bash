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

# PLOT_DEPTH_SLICE

# USAGE
#   ./utils/plot_depth_slice.bash PARAMETER DEPTH CBMIN CBMAX

# EXAMPLE
#   ./utils/plot_depth_slice.bash vsv 100

# COMMAND LINE ARGUMENTS
#   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)
#   DEPTH                  - depth in which you created the depth slice
#   CBMIN (optional)       - minimum of the color bar
#   CBMAX (optional)       - maximum of the color bar

# DESCRIPTION
#   Plots the depth slices for the grid files created by GLL2LL or by CREATE_EXTRA_DEPTH_SLICE.

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
  echo " PLOT_DEPTH_SLICE"
  echo ""
  echo " USAGE"
  echo "   ./utils/plot_depth_slice.bash PARAMETER DEPTH CBMIN CBMAX"
  echo ""
  echo " EXAMPLE"
  echo "   ./utils/plot_depth_slice.bash vsv 100"
  echo ""
  echo " COMMAND LINE ARGUMENTS"
  echo "   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)"
  echo "   DEPTH                  - depth in which you created the depth slice"
  echo "   CBMIN (optional)       - minimum of the color bar"
  echo "   CBMAX (optional)       - maximum of the color bar"
  echo " DESCRIPTION"
  echo "   Plots the depth slices for the grid files created by GLL2LL or by CREATE_EXTRA_DEPTH_SLICE."
  echo ""
}

if [ "$#" -ne 2 ] && [ "$#" -ne 4 ]; then
  help
  exit 1
fi

filename="$1"_"$2"_DS.dat
grdname="$1"_"$2"_DS.grd
output="$1"_"$2"_DS

np=$(head -n1 $filename | cut -f7 -d' ')
nt=$(head -n1 $filename | cut -f6 -d' ')

dp=$(echo "360 / ($np - 1)" | bc -l)
dt=$(echo "180 / ($nt - 1)" | bc -l)

info=$(./bin/getinfo $filename)

mean=$(echo $info | awk '{print $5}')
stdv=$(echo $info | awk '{print $6}')

cbmin=$(echo "$mean - 3 * $stdv" | bc -l)
cbmax=$(echo "$mean + 3 * $stdv" | bc -l)

if [ "$#" -eq 4 ]; then
  cbmin=$3
  cbmax=$4
fi

string=$(gmt info $filename | cut -f4 | awk -F '[</>]' '{print $2" "$3}')

min=$(echo $string | awk '{print $1}')
max=$(echo $string | awk '{print $2}')

awk -v min=$min -v max=$max 'BEGIN {printf "Min = %E Max = %E\n", min, max}'
echo 'Creating figure...'

gmt xyz2grd $filename -G$grdname -Rg -I$dp/$dt -: -h3
gmt makecpt -Cextra/tomo.cpt -T$cbmin/$cbmax > tomo_rescaled.cpt

gmt begin $output pdf
gmt set FONT_TITLE 20p,100

  gmt grdimage $grdname -JR0/12c -Rg -Bxa90fg90 -Bya30fg30 -BWeSn+t"$title ($2 km depth)" -Ctomo_rescaled.cpt

  gmt coast -W0.2,50

  gmt plot extra/plate_boundaries.dat -: -W0.7p
  gmt plot extra/hotspots.dat -Sc0.15c -Ggreen -W0.03c,black -h1 -l"Hotspots"

  gmt colorbar -Ctomo_rescaled.cpt -Baf -DJBC+e -B+l"$label"
  gmt legend -Dn0.85/0.95+w2.3c -C0.2c/0.2c -F+p+gazure1+r+s
gmt end

echo 'Figure saved!'
