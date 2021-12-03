/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/config.h"

void atlas_omp_set_num_threads(int num_threads);
int atlas_omp_get_num_threads(void);
int atlas_omp_get_max_threads(void);
int atlas_omp_get_thread_num(void);
int atlas_omp_get_num_procs(void);
int atlas_omp_in_parallel(void);
void atlas_omp_set_dynamic(int dynamic_threads);
int atlas_omp_get_dynamic(void);
void atlas_omp_set_nested(int nested);
int atlas_omp_get_nested(void);

#if ATLAS_HAVE_OMP
#define ATLAS_OMP_STR(x) #x
#define ATLAS_OMP_STRINGIFY(x) ATLAS_OMP_STR(x)
#define atlas_omp_pragma(x) _Pragma(ATLAS_OMP_STRINGIFY(x))
#else
#define atlas_omp_pragma(x)
#endif

#define atlas_omp_parallel_for atlas_omp_pragma(omp parallel for schedule(guided) ) for
#define atlas_omp_for atlas_omp_pragma(omp for schedule(guided)) for
#define atlas_omp_parallel atlas_omp_pragma(omp parallel)
#define atlas_omp_critical atlas_omp_pragma(omp critical)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <typename T>
class atlas_omp_scoped_helper {
public:
    atlas_omp_scoped_helper(T p): value(p), once_(true) {}
    bool once() const { return once_; }
    void done() { once_ = false; }
    T value;

private:
    bool once_;
};

#define _atlas_omp_scoped(T, VAR, VAL) for (atlas_omp_scoped_helper<T> VAR(VAL); VAR.once(); VAR.done())
#endif

#define atlas_omp_critical_ordered                                          \
    _atlas_omp_scoped(const size_t, _nthreads, atlas_omp_get_num_threads()) \
    atlas_omp_pragma( omp for ordered schedule(static,1) )\
    for( size_t _thread=0; _thread<_nthreads.value; ++_thread )\
      atlas_omp_pragma( omp ordered )

#undef _atlas_omp_scoped
