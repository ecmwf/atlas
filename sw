#! /usr/bin/env bash

# This script compiles the sources, runs the shallow_water executable
# and removes the *.mod files

# Environment variables that must be defined: 
# - FC
# - GRIB_API_DIR

# NOTE: make sure that grib_api is compiled with same compiler as $FC

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
GRIB_API_DIR=/lus/scratch/ecmwf/esm/grib_api/1.11.0/cray/82/

GRIB="\
  -I$GRIB_API_DIR/include\
  -L$GRIB_API_DIR/lib -lgrib_api_f90 -lgrib_api \
"

ftn -O0 -emf $SRC $GRIB -o shallow_water 
  rm *.mod
