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

#ifdef ATLAS_HAVE_OMP
#include <omp.h>
#else
/// Minimal set of stubs (see atlas_omp.cc)
inline void omp_set_num_threads(int num_threads)
{
}
inline int omp_get_num_threads(void)
{
 return 1;
}
inline int omp_get_max_threads(void)
{
 return 1;
}
inline int omp_get_thread_num(void)
{
 return 0;
}
inline int omp_get_num_procs(void)
{
 return 1;
}
inline int omp_in_parallel(void)
{
 return 0;
}
inline void omp_set_dynamic(int dynamic_threads)
{
}
inline int omp_get_dynamic(void)
{
 return 0;
}
inline void omp_set_nested(int nested)
{
}
inline int omp_get_nested(void)
{
 return 0;
}

#endif

#endif
