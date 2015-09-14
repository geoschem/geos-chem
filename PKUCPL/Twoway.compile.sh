#!/bin/bash

#Record the beginning time for compiling
date
declare -i account=0

#Set the two-way coupled run (used in Makefile_header.mk)
COUPLE="COUPLE=yes"
#COUPLECH="COUPLECH=4x5ch"
#COUPLENA="COUPLENA=4x5na"
#COUPLEEU="COUPLEEU=4x5eu"
COUPLECH="COUPLECH=2x25ch"
COUPLENA="COUPLENA=2x25na"
COUPLEEU="COUPLEEU=2x25eu"

#Set the nesting option
NESTEU="NEST=eu"
NESTNA="NEST=na"
NESTCH="NEST=ch"

#Set the MET input
MET="MET=geos5"

#Set the GRID input
GRID="GRID=2x25"
#GRID="GRID=4x5"
GRIDNEST="GRID=05x0666"

#Set the name for geos
geos="geos"
geosna="geosna"
geoseu="geoseu"
geosch="geosch"

#Option for making
option="-j4"

#Set the code dir and run dir
dir_code="/home/yanyy/geoschem/run/src/Code.v9-02.two-way/"
dir_coupler="PKUCPL/"
dir_geos="geos/"
dir_geosna="geosna/"
dir_geoseu="geoseu/"
dir_geosch="geosch/"
dir_bin="bin/"
dir_compile='Twoway.compile.sh'

#Compile the code for nested model in Aisa
if [ ${#COUPLECH} != 0 ]; then
rm -f $dir_code$dir_coupler$dir_geosch$geos
cd $dir_code
make realclean
make $option $MET $GRIDNEST $NESTCH $COUPLE $COUPLECH
cp $dir_code$dir_bin$geos $dir_code$dir_coupler$dir_geosch
fi

#Compile the code for nested model in North America
if [ ${#COUPLENA} != 0 ]; then
rm -f $dir_code$dir_coupler$dir_geosna$geos
cd $dir_code
make realclean
make $option $MET $GRIDNEST $NESTNA $COUPLE $COUPLENA
cp $dir_code$dir_bin$geos $dir_code$dir_coupler$dir_geosna
fi

#Compile the code for nested model in Europe
if [ ${#COUPLEEU} != 0 ]; then
rm -f $dir_code$dir_coupler$dir_geoseu$geos
cd $dir_code
make realclean
make $option $MET $GRIDNEST $NESTEU $COUPLE $COUPLEEU
cp $dir_code$dir_bin$geos $dir_code$dir_coupler$dir_geoseu
fi

#Compile the code for global model
rm -f $dir_code$dir_coupler$dir_geos$geos
cd $dir_code
make realclean
make $option $MET $GRID $COUPLE $COUPLENA $COUPLECH $COUPLEEU
cp $dir_code$dir_bin$geos $dir_code$dir_coupler$dir_geos

#Compile the PKUCPL
cd $dir_code$dir_coupler
rm -f PKUCPL.o run
#ifort -O2 -o run PKUCPL.F90
ifort -g -check bounds -traceback -o run PKUCPL.F90

