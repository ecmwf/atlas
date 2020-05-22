/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/os/BackTrace.h"

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

#if ATLAS_HAVE_FORTRAN
#include "fckit/Log.h"
#endif

namespace atlas {

std::string backtrace() {
    return eckit::BackTrace::dump();
}

namespace detail {

void debug_parallel_here( const eckit::CodeLocation& here ) {
    const auto& comm = mpi::comm();
    comm.barrier();
    Log::info() << "DEBUG_PARALLEL() @ " << here << std::endl;
}

void debug_parallel_what( const eckit::CodeLocation& here, const std::string& what ) {
    const auto& comm = mpi::comm();
    comm.barrier();
    Log::info() << "DEBUG_PARALLEL(" << what << ") @ " << here << std::endl;
}

}  // namespace detail

Log::Channel& Log::info() {
    return atlas::Library::instance().infoChannel();
}

Log::Channel& Log::warning() {
    return atlas::Library::instance().warningChannel();
}

Log::Channel& Log::trace() {
    return atlas::Library::instance().traceChannel();
}

Log::Channel& Log::debug() {
    return atlas::Library::instance().debugChannel();
}

void Log::addFortranUnit( int unit, Style style, const char* prefix ) {
#if ATLAS_HAVE_FORTRAN
    fckit::Log::addFortranUnit( unit, fckit::Log::Style( style ), prefix );
#else
/*NOTIMP*/
#endif
}

void Log::setFortranUnit( int unit, Style style, const char* prefix ) {
#if ATLAS_HAVE_FORTRAN
    fckit::Log::setFortranUnit( unit, fckit::Log::Style( style ), prefix );
#else
/*NOTIMP*/
#endif
}

}  // namespace atlas
