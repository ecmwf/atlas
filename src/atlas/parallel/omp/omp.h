/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_omp_h
#define atlas_omp_h

#include "atlas/library/config.h"

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

#define ATLAS_STR(x) #x
#define ATLAS_STRINGIFY(x) ATLAS_STR(x)
#define ATLAS_CONCATENATE(X,Y) X Y

#ifdef ATLAS_HAVE_OMP
#define atlas_omp_pragma(x) \
  _Pragma( ATLAS_STRINGIFY( x ) )
#else
#define atlas_omp_pragma(x)
#endif

#define atlas_omp_parallel_for atlas_omp_pragma(omp parallel for) for
#define atlas_omp_for atlas_omp_pragma(omp for) for
#define atlas_omp_parallel atlas_omp_pragma(omp parallel)
#define atlas_omp_critical atlas_omp_pragma(omp critical)

template <typename T>
class atlas_scoped_helper
{
public:
  atlas_scoped_helper(T p): value(p), once_(true) {}
  bool once() const { return once_; }
  void done() { once_=false; }
  T value;
private:
  bool once_;
};
#define atlas_scoped(T,VAR,VAL) \
  for( atlas_scoped_helper<T> VAR(VAL); VAR.once(); VAR.done() )

#define atlas_omp_critical_ordered \
    atlas_scoped(const size_t, _nthreads, omp_get_num_threads()) \
    atlas_omp_pragma( omp for ordered schedule(static,1) )\
    for( size_t _thread=0; _thread<_nthreads.value; ++_thread )\
      atlas_omp_pragma( omp ordered )

#endif
