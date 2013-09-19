#! /usr/bin/env bash

if [[ $(uname -n) == "coolcat.local" ]]; then
  FC=ifort
else
  FC=/usr/local/apps/intel/parallel_studio_xe_2013/bin/ifort
fi
FC=/usr/local/apps/gcc/4.8.1/LP64/bin/gfortran
$FC -g -O0 \
  src/mesh/elements_module.f90\
  src/mesh/lagrangep0_module.f90\
  src/mesh/lagrangep1_module.f90\
  src/mesh/grid_module.f90\
  src/solver/model_module.f90\
  src/solver/mpdata_module.f90\
  src/solver/shallow_water_module.f90\
  src/io/read_joanna_module.f90\
  src/io/gmsh_module.f90\
  src/shallow_water.f90\
  -o shallow_water &&\
  ./shallow_water &&\
  rm *.mod shallow_water
