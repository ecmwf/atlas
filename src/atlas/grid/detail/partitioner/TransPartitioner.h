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

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

/// @class TransPartitioner
/// @brief Equal regions partitioning algorithm computed by the IFS trans
/// library
///
/// The user is advised to use the "EqualRegionsPartitioner" class instead. This
/// implementation is only here to to guarantee the exact same distribution
/// as IFS is using. The difference with "EqualRegionsPartitioner" is minimal.
/// (a few points may be assigned to different partitions).
class TransPartitioner : public Partitioner {
public:
    /// @brief Constructor
    TransPartitioner();

    TransPartitioner(const eckit::Parametrisation&);

    TransPartitioner(const idx_t nb_partitions, const eckit::Parametrisation& = util::NoConfig());

    virtual ~TransPartitioner();

    using Partitioner::partition;
    /// Warning: this function temporariliy allocates a new Trans, but without the
    /// computations
    /// of the spectral coefficients (LDGRIDONLY=TRUE)
    virtual void partition(const Grid&, int part[]) const;

    int nb_bands() const;

    int nb_regions(int b) const;

    virtual std::string type() const { return "ectrans"; }

private:
    size_t nbands_;
    std::vector<size_t> nregions_;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
