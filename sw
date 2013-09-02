#! /usr/bin/env bash

if [ $(uname -n)=="coolcat.local" ]; then
  FC=gfortran
  echo mac
else
  FC=/usr/local/apps/intel/parallel_studio_xe_2013/bin/ifort
fi

$FC -g -O0 \
  mesh/elements_module.f90\
  mesh/lagrangep0_module.f90\
  mesh/lagrangep1_module.f90\
  mesh/grid_module.f90\
  solver/model_module.f90\
  solver/mpdata_module.f90\
  solver/shallow_water_module.f90\
  io/read_joanna_module.f90\
  io/gmsh_module.f90\
  shallow_water.f90\
  -o shallow_water &&\
  ./shallow_water &&\
  rm *.mod shallow_water