/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "FFT.h"

#include <cstdlib>

#include "atlas/library/defines.h"


#if ATLAS_HAVE_POCKETFFT

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Trace.h"

ATLAS_SUPPRESS_WARNINGS_PUSH
ATLAS_SUPPRESS_WARNINGS_BOTH_INLINE_NOINLINE
#define POCKETFFT_CACHE_SIZE 4000
// POCKETFFT_CACHE_SIZE should keep O8000 grid fft's in cache
#include "pocketfft_hdronly.h"
ATLAS_SUPPRESS_WARNINGS_PUSH

namespace atlas::linalg {

    void pocketfft::do_plan_inverse_c2r(size_t /*size_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
    void pocketfft::do_plan_inverse_c2r_many(size_t /*howmany*/, size_t /*size_in*/, size_t /*size_out*/, size_t /*dist_in*/, size_t /*dist_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
    void pocketfft::do_inverse_c2r(const size_t size_out, std::complex<double>* in, double* out) const {
        ATLAS_TRACE("pocketfft::inverse_c2r");
        static const ::pocketfft::stride_t stride_in{sizeof(std::complex<double>)}; // in bytes
        static const ::pocketfft::stride_t stride_out{sizeof(double)}; // in bytes
        constexpr size_t axis = 0;
        constexpr double fct = 1.0;
        constexpr size_t nthreads = 1;
        constexpr bool direction = ::pocketfft::BACKWARD;
        ::pocketfft::shape_t shape_out{size_out};
        ::pocketfft::c2r(shape_out, stride_in, stride_out, axis, direction, in, out, fct, nthreads);
    }
    void pocketfft::do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
        ATLAS_TRACE("pocketfft::inverse_c2r");
        static const ::pocketfft::stride_t stride_in{sizeof(std::complex<double>)}; // in bytes
        static const ::pocketfft::stride_t stride_out{sizeof(double)}; // in bytes
        const size_t idist = dist_in;
        const size_t odist = dist_out;
        constexpr size_t axis = 0;
        constexpr double fct = 1.0;
        constexpr size_t nthreads = 1;
        constexpr bool direction = ::pocketfft::BACKWARD;
        const ::pocketfft::shape_t shape_out{size_out};

        atlas_omp_parallel_for(size_t j=0; j<howmany; ++j) {
            ::pocketfft::c2r(shape_out, stride_in, stride_out, axis, direction, in+j*idist, out+j*odist, fct, nthreads);
        }
    }

}  // namespace atlas::linalg

#else

#include "atlas/runtime/Exception.h"
#define THROW_ATLAS_POCKETFFT_NOT_SUPPORTED() throw_Exception("Atlas was not compiled with pocketfft support", Here())

namespace atlas::linalg {

    void pocketfft::do_plan_inverse_c2r(size_t /*size_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
    void pocketfft::do_plan_inverse_c2r_many(size_t /*howmany*/, size_t /*size_in*/, size_t /*size_out*/, size_t /*dist_in*/, size_t /*dist_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
    void pocketfft::do_inverse_c2r(const size_t /*size_out*/, std::complex<double>* /*in*/, double* /*out*/) const {
        THROW_ATLAS_POCKETFFT_NOT_SUPPORTED();
    }
    void pocketfft::do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
        THROW_ATLAS_POCKETFFT_NOT_SUPPORTED();
    }

}  // namespace atlas::linalg

#endif
