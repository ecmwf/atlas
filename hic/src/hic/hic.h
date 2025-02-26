/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include "hic/hic_config.h"
#include "hic/hic_runtime.h"

#include <sstream>
#include <stdexcept>

#if HIC_BACKEND_DUMMY
#define HIC_CALL(val)
#else
#define HIC_CALL(val) hic_assert((val), #val, __FILE__, __LINE__)
#endif

inline void hic_assert(hicError_t err, const char* const func, const char* const file, const int line) {
    // Ignore errors when HIP/CUDA runtime is unloaded or deinitialized.
    // This happens when calling HIP/CUDA after main has ended, e.g. in teardown of static variables calling `hicFree`
    //   --> ignore hicErrorDeinitialized (a.k.a. cudaErrorCudartUnloading / hipErrorDeinitialized)
    if (err != hicSuccess && err != hicErrorDeinitialized) {
        std::ostringstream msg;
        msg << "HIC Runtime Error [code=" << err << "] at: " << file << " + " << line << " : " << func << "\n";
        msg << "  Reason: " << hicGetErrorString(err);
        throw std::runtime_error(msg.str());
    }
}

#define HIC_CHECK_KERNEL_LAUNCH() HIC_CALL(hicPeekAtLastError())
#define HIC_CHECK_LAST_ERROR() HIC_CALL(hicGetLastError())
