#! /usr/bin/env bash

FC=/usr/local/apps/gcc/4.8.1/LP64/bin/gfortran
#FC=/usr/local/apps/intel/parallel_studio_xe_2013/bin/ifort

GRIB_API_DIR=/home/rd/nawd/local
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRIB_API_DIR/lib

SRC="\
  src/common/common_module.f90\
  src/mesh/elements_module.f90\
  src/mesh/lagrangep0_module.f90\
  src/mesh/lagrangep1_module.f90\
  src/mesh/grid_module.f90\
  src/solver/model_module.f90\
  src/solver/mpdata_module.f90\
  src/solver/shallow_water_module.f90\
  src/io/read_joanna_module.f90\
  src/io/gmsh_module.f90\
  src/io/grib_module.f90\
  src/shallow_water.f90\
"

GRIB="\
  -I$GRIB_API_DIR/include\
  -L$GRIB_API_DIR/lib -lgrib_api_f90 -lgrib_api -ljasper\
"
$FC -O3 $SRC $GRIB -o shallow_water &&\
  ./shallow_water &&\
  rm *.mod shallow_water
