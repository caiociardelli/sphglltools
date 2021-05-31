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

# PLOT_VERTICAL_CROSS_SECTION

# USAGE
#   ./utils/plot_vertical_cross_section.bash PARAMETER CBMIN CBMAX

# EXAMPLE
#   ./utils/plot_vertical_cross_section.bash vpvs

# COMMAND LINE ARGUMENTS
#   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)
#   CBMIN (optional)       - minimum of the color bar
#   CBMAX (optional)       - maximum of the color bar

# DESCRIPTION
#   Plots the vertical cross-section for the grid files created by GLL2DD or by
#   CREATE_EXTRA_VERTICAL_CROSS_SECTION.

#-----------------------------------------------------------------------------------------------

title="GLAD-M15"

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
  echo " PLOT_VERTICAL_CROSS_SECTION"
  echo ""
  echo " USAGE"
  echo "   ./utils/plot_vertical_cross_section.bash PARAMETER"
  echo ""
  echo " EXAMPLE"
  echo "   ./utils/plot_vertical_cross_section.bash vp"
  echo ""
  echo " COMMAND LINE ARGUMENTS"
  echo "   PARAMETER              - model parameter to be plotted (cb, vpvs, ti, dcb, dvsv, etc.)"
  echo "   CBMIN (optional)       - minimum of the color bar"
  echo "   CBMAX (optional)       - maximum of the color bar"
  echo ""
  echo " DESCRIPTION"
  echo "   Plots the vertical cross-section for the grid files created by GLL2DD or by"
  echo "   CREATE_EXTRA_VERTICAL_CROSS_SECTION."
  echo ""
}

if [ "$#" -ne 1 ] && [ "$#" -ne 3 ]; then
  help
  exit 1
fi

REARTH=6371
R410=5961
R650=5721
KM2DEGREE=111.2
TOPOSCALE=0.08

filename="$1"_VCS.dat
grdname="$1"_VCS.grd
output="$1"_VCS

nd=$(head -n1 $filename | cut -f4 -d' ')
nr=$(head -n1 $filename | cut -f3 -d' ')

line=$(sed '2q;d' $filename)

lat1=$(echo $line | awk '{printf "%lf", $5}')
lon1=$(echo $line | awk '{printf "%lf", $6}')
lat2=$(echo $line | awk '{printf "%lf", $7}')
lon2=$(echo $line | awk '{printf "%lf", $8}')

lat0=$(echo "($lat1 + $lat2) / 2" | bc -l)
lon0=$(echo "($lon1 + $lon2) / 2" | bc -l)

info=$(./bin/getinfo $filename)

mean=$(echo $info | awk '{printf "%lf", $5}')
stdv=$(echo $info | awk '{printf "%lf", $6}')

if [ "$#" -eq 3 ]; then
  cbmin=$2
  cbmax=$3
else
  cbmin=$(echo "$mean - 3 * $stdv" | bc -l)
  cbmax=$(echo "$mean + 3 * $stdv" | bc -l)
fi

string=$(gmt info $filename)

dstring=$(echo "$string" | cut -f3 | awk -F '[</>]' '{print $2" "$3}')
rstring=$(echo "$string" | cut -f2 | awk -F '[</>]' '{print $2" "$3}')
vstring=$(echo "$string" | cut -f4 | awk -F '[</>]' '{print $2" "$3}')

drange=$(echo $dstring | awk '{printf "%lf", $2 - $1}')
rrange=$(echo $rstring | awk '{printf "%lf", $2 - $1}')

if [ $drange > 180 ]; then
  lon0=$(echo "$lon0 + 180" | bc -l)
fi

dd=$(echo "$drange / ($nd - 1)" | bc -l)
dr=$(echo "$rrange / ($nr - 1)" | bc -l)

dmin=$(echo "$dstring" | awk '{printf "%lf", $1}')
dmax=$(echo "$dstring" | awk '{printf "%lf", $2}')
rmin=$(echo "$rstring" | awk '{printf "%lf", $1}')
rmax=$(echo "$rstring" | awk '{printf "%lf", $2}')
vmin=$(echo "$vstring" | awk '{printf "%lf", $1}')
vmax=$(echo "$vstring" | awk '{printf "%lf", $2}')

ermax=$(echo "1.01 * $rmax" | bc -l)

dmean=$(echo "($dmin + $dmax) / 2" | bc -l)
rmean=$(echo "($rmin + $rmax) / 2" | bc -l)

dlabeld=$(echo "$dmean" | bc -l)
dlabelr=$(echo "0.77 * $rmin" | bc -l)
rlabeld=$(echo "$dmin - 0.25 * $drange" | bc -l)
rlabelr=$(echo "1.2 * $rmean" | bc -l)

gmt project -C$lon1/$lat1 -E$lon2/$lat2 -G10 -Q > great_circle_points.xyp
gmt grdtrack -Gextra/earth_relief_15m.grd great_circle_points.xyp | awk '{print $3" "$4}' > profile.dat

string=$(gmt info profile.dat | cut -f3 | awk -F '[</>]' '{print $2" "$3}')

emin=$(echo $string | awk '{printf "%lf", $1}')
emax=$(echo $string | awk '{printf "%lf", $2}')
erange=$(echo $string | awk '{printf "%lf", $2 - $1}')
scale=$(echo "$erange / ($TOPOSCALE * $rrange)" | bc -l)
shift=$(echo "((0.2 * $erange) - $emin) / $scale" | bc -l)
Rmax=$(echo "1.015 * ($rmax + $shift)" | bc -l)

awk -v km2dg=$KM2DEGREE -v rmax=$rmax -v shift=$shift -v scale=$scale\
    '{print $1 / km2dg" "rmax + shift + $2 / scale}' profile.dat > profile_rescaled.dat
awk -v dmin=$dmin -v dmax=$dmax -v rmax=$rmax\
    'BEGIN {for (delta = dmax; delta >= 0; delta -= 0.1) printf "%.3f %.1f\n", delta, rmax}' >> profile_rescaled.dat

awk -v min=$vmin -v max=$vmax 'BEGIN {printf "Min = %E Max = %E\n", min, max}'
echo 'Creating figure...'

gmt xyz2grd $filename -G$grdname -R$dmin/$dmax/$rmin/$rmax -I$dd/$dr -: -h4
gmt makecpt -Cextra/tomo.cpt -T$cbmin/$cbmax > tomo_rescaled.cpt

gmt begin $output pdf
  gmt set FONT_TITLE 20p,100

  gmt grdimage $grdname -JPa10c/$dmean -R$dmin/$dmax/$rmin/$Rmax -BWeS+t$title -Baf -Ctomo_rescaled.cpt

  gmt plot profile_rescaled.dat -W0.03c -L -Glightgray

  gmt text -N <<< "$dlabeld $dlabelr @~\104@~"
  gmt text -F+a$dmean -N <<< "$rlabeld $rlabelr r (km)"

  echo "$dmin $Rmax" | gmt plot -Skrflag/1.2c -Ggreen -W0.03c,black -N
  echo "$dmax $Rmax" | gmt plot -Skrflag/1.2c -Gmagenta -W0.03c,black -N

  gmt plot -W0.02c,0,-- << DISCONTINUITIES
$(for delta in $(seq 0 360); do echo $delta $R410; done)
$(for delta in $(seq 0 360); do echo $delta $R650; done)
DISCONTINUITIES

  gmt colorbar -Ctomo_rescaled.cpt -Baf -DJBC+e -B+l"$label"

  gmt inset begin -DJTR+w10c -M0.1c
    gmt coast -JG$lon0/$lat0/2.5c -Rg -Bafg -Dc -A10000 -Glightgray -Wthinnest

    gmt plot extra/plate_boundaries.dat -JG$lon0/$lat0/2.5c -: -W0.6p,brown
    gmt plot great_circle_points.xyp -JG$lon0/$lat0/2.5c -W1.0p

    echo "$lon1 $lat1" | gmt plot -JG$lon0/$lat0/2.5c -Skrflag/0.8c -Ggreen -W0.03c,black -N
    echo "$lon2 $lat2" | gmt plot -JG$lon0/$lat0/2.5c -Skrflag/0.8c -Gmagenta -W0.03c,black -N
  gmt inset end
gmt end

echo 'Figure saved!'
