/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_omp_h
#define atlas_omp_h

#include "atlas/atlas_config.h"

#ifdef HAVE_OMP
#include <omp.h>
#else
  #include "eckit/log/Log.h"
  inline void omp_set_num_threads(int) { eckit::Log::warning() << "\nWARNING: OpenMP not available!\n" << std::endl; }
  inline int omp_get_max_threads() { return 0; }
#endif

#endif
