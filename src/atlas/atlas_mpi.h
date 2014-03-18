#ifndef atlas_mpi_h
#define atlas_mpi_h

#include "atlas/atlas_config.h"
#include "atlas/atlas_defines.h"

#ifdef HAVE_MPI
#include <mpi.h>
#else
#include "atlas/mpistubs/mpi.h"
#endif

#endif
