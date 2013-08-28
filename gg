#! /usr/bin/env bash
/usr/local/apps/intel/parallel_studio_xe_2013/bin/ifort -g -O0 \
  elements_module.f90\
  lagrangep0_module.f90\
  lagrangep1_module.f90\
  grid_module.f90\
  state_module.f90\
  read_joanna_module.f90\
  gmsh_module.f90\
  mpdata_module.f90\
  shallow_water_module.f90\
  shallow_water.f90\
  -o shallow_water\
  && ./shallow_water && rm *.mod && rm shallow_water