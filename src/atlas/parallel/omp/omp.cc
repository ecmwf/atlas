/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/parallel/omp/omp.h"

#ifdef _OPENMP
#include <omp.h>
#else
extern void omp_set_num_threads(int num_threads);
extern int omp_get_num_threads(void);
extern int omp_get_max_threads(void);
extern int omp_get_thread_num(void);
extern int omp_get_num_procs(void);
extern int omp_in_parallel(void);
extern int omp_set_dynamic(int dynamic_threads);
extern int omp_get_dynamic(void);
extern void omp_set_nested(int nested);
extern int omp_get_nested(void);
#endif

#pragma weak omp_set_num_threads
#pragma weak omp_get_num_threads
#pragma weak omp_get_max_threads
#pragma weak omp_get_thread_num
#pragma weak omp_get_num_procs
#pragma weak omp_in_parallel
#pragma weak omp_set_dynamic
#pragma weak omp_get_dynamic
#pragma weak omp_set_nested
#pragma weak omp_get_nested

void atlas_omp_set_num_threads(int num_threads)
{
  if( omp_set_num_threads ) omp_set_num_threads(num_threads);
}

int atlas_omp_get_num_threads(void)
{
  return omp_get_num_threads ? omp_get_num_threads() : 1;
}

int atlas_omp_get_max_threads(void)
{
  return omp_get_max_threads ? omp_get_max_threads() : 1;
}

int atlas_omp_get_thread_num(void)
{
  return omp_get_thread_num ? omp_get_thread_num() : 0;
}

int atlas_omp_get_num_procs(void)
{
  return omp_get_num_procs ? omp_get_num_procs() : 1;
}

int atlas_omp_in_parallel(void)
{
  return omp_in_parallel ? omp_in_parallel() : 0;
}

void atlas_omp_set_dynamic(int dynamic_threads)
{
  if( omp_set_dynamic ) omp_set_dynamic(dynamic_threads);
}

int atlas_omp_get_dynamic(void)
{
  return omp_get_dynamic ? omp_get_dynamic() : 0;
}

void atlas_omp_set_nested(int nested)
{
  if( omp_set_nested ) omp_set_nested(nested);
}

int atlas_omp_get_nested(void)
{
  return omp_get_nested ? omp_get_nested() : 0;
}

