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
#include <mutex>
#include <sstream>

#include "eckit/log/Bytes.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "atlas/library/defines.h"

#include "atlas/parallel/omp/omp.h"

#if ATLAS_HAVE_FFTW
#include "fftw3.h"
#endif

ATLAS_SUPPRESS_WARNINGS_PUSH
ATLAS_SUPPRESS_WARNINGS_BOTH_INLINE_NOINLINE
#define POCKETFFT_CACHE_SIZE 4000
// POCKETFFT_CACHE_SIZE should keep O8000 grid fft's in cache
#include "pocketfft_hdronly.h"
ATLAS_SUPPRESS_WARNINGS_PUSH

namespace atlas::linalg {

    void pocketfft::do_prepare_inverse_c2r(size_t /*size_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
    void pocketfft::do_prepare_inverse_c2r_many(size_t /*howmany*/, size_t /*size_in*/, size_t /*size_out*/, size_t /*dist_in*/, size_t /*dist_out*/, std::complex<double>* /*in*/, double* /*out*/) const {}
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

#define THROW_ATLAS_FFTW_NOT_SUPPORTED() throw_Exception("Atlas was not compiled with FFTW support", Here())

struct FFTW_Plan {
#if ATLAS_HAVE_FFTW
    ~FFTW_Plan() {
        destroy();
    }
    void create_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) {
        ATLAS_TRACE("FFTW::create_plan_inverse_c2r");
        plan = fftw_plan_dft_c2r_1d(size_out, reinterpret_cast<fftw_complex*>(in), out, fftw_plan_flags);
        constructed_ = true;
    }
    void create_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) {
        ATLAS_TRACE("FFTW::create_plan_inverse_c2r");
        constexpr int rank = 1; /* 1d transforms */
        const int n[] = {static_cast<int>(size_out)}; /* each transform of length size_out */
        constexpr int istride = 1;
        constexpr int ostride = 1;
        const int idist = dist_in;
        const int odist = dist_out;
        const int* inembed = nullptr;
        const int* onembed = nullptr;
        plan = fftw_plan_many_dft_c2r(rank, n, howmany, reinterpret_cast<fftw_complex*>(in),
            inembed, istride, idist, out, onembed, ostride, odist, fftw_plan_flags);
        constructed_ = true;
    }

    void inverse_c2r(std::complex<double>* in, double* out) {
        ATLAS_TRACE("FFTW::inverse_c2r");
        ATLAS_ASSERT(constructed_);
        fftw_execute_dft_c2r(plan, reinterpret_cast<fftw_complex*>(in), out);
    }

    void destroy() {
        if (constructed_) {
            fftw_destroy_plan(plan);
        }
        constructed_ = false;
    }
    fftw_plan plan;
    bool constructed_{false};
    int fftw_plan_flags{FFTW_ESTIMATE};
#else
    void inverse_c2r(std::complex<double>* in, double* out) {
        THROW_ATLAS_FFTW_NOT_SUPPORTED();
    }
    void create_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) {
        THROW_ATLAS_FFTW_NOT_SUPPORTED();
    }
    void create_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) {
        THROW_ATLAS_FFTW_NOT_SUPPORTED();
    }
#endif
};

struct FFTW::FFTW_Plans {
    mutable std::map<std::int64_t, std::unique_ptr<FFTW_Plan>> fftw_plans_;
    mutable std::mutex fftw_plans_mutex;  // protects fftw_plans_

    FFTW_Plan& get_plan(size_t size_out, std::complex<double>* in, double* out) const {
        std::lock_guard lock(fftw_plans_mutex);

        bool inverse = true;
        bool inplace = (static_cast<void*>(in) == static_cast<void*>(out));
        bool aligned = ((reinterpret_cast<size_t>(in) & 15) | (reinterpret_cast<size_t>(out) & 15)) == 0;
        int n1 = size_out;
        std::int64_t key = ((n1 << 3) | (inverse << 2) | (inplace << 1) | aligned) << 1;
        auto it = fftw_plans_.find(key);
        if (it == fftw_plans_.end()) {
            it = fftw_plans_.emplace(key, new FFTW_Plan()).first;
            it->second->create_inverse_c2r(size_out, in, out);
        }
        return *(it->second);
    }

    FFTW_Plan& get_plan_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
        std::lock_guard lock(fftw_plans_mutex);

        bool inverse = true;
        bool inplace = (static_cast<void*>(in) == static_cast<void*>(out));
        bool aligned = ((reinterpret_cast<size_t>(in) & 15) | (reinterpret_cast<size_t>(out) & 15)) == 0;
        int n0 = howmany;
        int n1 = size_out;
        std::int64_t key = (((((std::int64_t)n0) << 30) | (n1 << 3) | (inverse << 2) | (inplace << 1) | aligned) << 1) + 1;
        auto it = fftw_plans_.find(key);
        if (it == fftw_plans_.end()) {
            it = fftw_plans_.emplace(key, new FFTW_Plan()).first;
            it->second->create_inverse_c2r_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
        }
        return *(it->second);
    }
};

FFTW::FFTW() {
    plans_ = new FFTW_Plans();
}

FFTW::~FFTW() {
    delete plans_;
}

void FFTW::do_prepare_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const {
    plans_->get_plan(size_out, in, out);
}
void FFTW::do_prepare_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
    plans_->get_plan_many(howmany, size_in, size_out, dist_in, dist_out, in, out);
}
void FFTW::do_inverse_c2r(size_t size_out, std::complex<double>* in, double* out) const {
    plans_->get_plan(size_out, in, out).inverse_c2r(in, out);
}
void FFTW::do_inverse_c2r_many(size_t howmany, size_t size_in, size_t size_out, size_t dist_in, size_t dist_out, std::complex<double>* in, double* out) const {
    plans_->get_plan_many(howmany, size_in, size_out, dist_in, dist_out, in, out).inverse_c2r(in, out);
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace atlas::linalg
