/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/parallel/Locate.h"
#include "atlas/field/Field.h"

namespace atlas::functionspace {

class Locator : public ::atlas::parallel::Locator {
public:
    Locator(const FunctionSpace& fs);

    using ::atlas::parallel::Locator::locate;
    void locate(
        // input
        span<const gidx_t> global_index, const gidx_t global_index_base,
        // output
        span<int> partition, const int partition_base,
        span<idx_t> remote_index, const idx_t remote_index_base) const override;

private:
    vector<int> distribution_array_;
    std::function<int(size_t)> distribution_function_;
    size_t distribution_size_;
    const Field fs_global_index_;
    const Field fs_ghost_;
    std::string mpi_comm_;
};

} // namespace atlas::functionspace
