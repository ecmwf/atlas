#! /usr/bin/env bash

# This makefile compiles the sources, for the shallow_water executable

# Environment variables that must be defined: 
# - FC (default=ftn)
# - GRIB_API_DIR can be commented out

# NOTE: make sure that grib_api is compiled with same compiler as $FC


# On Cray, you can find GRIB_API_DIR here
#GRIB_API_DIR=/lus/scratch/ecmwf/esm/grib_api/1.11.0/cray/82

ifdef GRIB_API_DIR
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRIB_API_DIR/lib
  GRIB= -DENABLE_GRIB -L$(GRIB_API_DIR)/lib -lgrib_api_f90 -lgrib_api -ljasper -I$(GRIB_API_DIR)/include
endif

#FC=ftn
FCFLAGS= -O3
# -hfp3 -rad
#-fstack-arrays
# -rad -Ktrap=fp -hfp3 

# Possible usefule flags for crayftn 8.2
# -hfp3  (float operations)
# -R b   (bounds checking)
# -rad   (optimisation listing, plus decompiled)
# -Ktrap=fp (for aborting at NaN)
#
# Possible extra flags to optimize gfortran 4.8
# -fstack-arrays (allocate automatic arrays on stack)

KERNEL_INC= -I./src/common -I./src/mesh -I./src/io  -I./ 
%.o : %.F90
	echo $@
	$(FC) $(GRIB) $(KERNEL_INC) $(FCFLAGS) -c -o $@ $<


KERNEL_SRC= \
	src/common/common_module.F90 \
	src/mesh/elements_module.F90 \
	src/mesh/lagrangep0_module.F90 \
	src/mesh/lagrangep1_module.F90 \
	src/mesh/grid_module.F90 \
	src/mesh/datastruct_module.F90 \
    src/io/read_joanna_module.F90 \
  	src/io/gmsh_module.F90 \
	src/io/grib_module.F90

KERNEL_OBJ=$(KERNEL_SRC:.F90=.o)

SHALLOW_WATER_SRC= \
	src/shallow_water.F90

SHALLOW_WATER_OBJ=$(SHALLOW_WATER_SRC:.F90=.o)

all: shallow_water

kernel: $(KERNEL_OBJ)

shallow_water: kernel $(SHALLOW_WATER_OBJ)
	$(FC) $(GRIB) $(KERNEL_INC) $(FCFLAGS) -o shallow_water $(KERNEL_OBJ) $(SHALLOW_WATER_OBJ)

clean:
	rm -f *.o             *.mod
	rm -f src/*.o         src/*.mod
	rm -f src/common/*.o  src/common/*.mod
	rm -f src/mesh/*.o    src/mesh/*.mod
	rm -f src/io/*.o      src/io/*.mod

