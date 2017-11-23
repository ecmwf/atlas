#!/bin/bash

# CC=mpicc
# CXX=mpicxx
# FC=mpif90
ulimit -s unlimited

# Initialise module environment if it is not
if [[ ! $(command -v module > /dev/null 2>&1) ]]; then
  . /usr/local/apps/module/init/bash
fi

# Minimum required cmake (3.1) when using fckit only
module unload cmake
module load cmake/3.7.1
