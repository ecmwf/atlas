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

#include <array>
#include <functional>
#include <type_traits>
#include <memory>

#include "atlas/array/DataType.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/functionspace/detail/StructuredColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/Vertical.h"
#include "atlas/library/config.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"
#include "atlas/util/vector.h"

namespace atlas {
namespace parallel {
class GatherScatter;
//class HaloExchange;
//class Checksum;
}  // namespace parallel
}  // namespace atlas

namespace atlas {
class Field;
//class FieldSet;
class Grid;
class StructuredGrid;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
class Partitioner;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {

// -------------------------------------------------------------------

class BlockStructuredColumns : public FunctionSpaceImpl {
public:
    BlockStructuredColumns(const Grid&, const eckit::Configuration& = util::NoConfig());

    BlockStructuredColumns(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());

    BlockStructuredColumns(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());

    BlockStructuredColumns(const Grid&, const grid::Distribution&, const Vertical&,
                      const eckit::Configuration& = util::NoConfig());

    BlockStructuredColumns(const Grid&, const Vertical&, const eckit::Configuration& = util::NoConfig());

    BlockStructuredColumns(const Grid&, const Vertical&, const grid::Partitioner&,
                      const eckit::Configuration& = util::NoConfig());

    static std::string static_type() { return "BlockStructuredColumns"; }
    std::string type() const override { return static_type(); }
    std::string distribution() const override { return structuredcolumns_->distribution(); }

    Field createField(const eckit::Configuration&) const override;
    Field createField(const Field&, const eckit::Configuration&) const override;

    using FunctionSpaceImpl::scatter;
    void scatter(const FieldSet&, FieldSet&) const override;
    void scatter(const Field&, Field&) const override;
    using FunctionSpaceImpl::gather;
    void gather(const FieldSet&, FieldSet&) const override;
    void gather(const Field&, Field&) const override;

    idx_t size() const override { return structuredcolumns_->size(); }
    idx_t index(idx_t jblk, idx_t jrof) const {
        return jblk * nproma_ + jrof; // local index;
    }
    idx_t nproma() const { return nproma_; }
    idx_t nblks() const { return nblks_; }

    const Vertical& vertical() const { return structuredcolumns_->vertical(); }
    const StructuredGrid& grid() const override { return structuredcolumns_->grid(); }

    idx_t levels() const { return structuredcolumns_->levels(); }
    Field lonlat() const override { return structuredcolumns_->lonlat(); }
    Field xy() const { return structuredcolumns_->xy(); }
    Field z() const { return structuredcolumns_->z(); }
    Field partition() const override { return structuredcolumns_->partition(); }
    Field global_index() const override { return structuredcolumns_->global_index(); }
    Field remote_index() const override { return structuredcolumns_->remote_index(); }
    Field index_i() const { return structuredcolumns_->index_i(); }
    Field index_j() const { return structuredcolumns_->index_j(); }
    Field ghost() const override { return structuredcolumns_->ghost(); }
    size_t footprint() const override { return structuredcolumns_->footprint(); }
    idx_t part() const override { return structuredcolumns_->part(); }
    idx_t nb_parts() const override { return structuredcolumns_->nb_parts(); }
    const StructuredColumns& structuredcolumns() const { return *structuredcolumns_; }
    idx_t block_begin(idx_t jblk) const { return jblk * nproma_; }
    idx_t block_size(idx_t jblk) const { return (jblk != nblks_-1 ? nproma_ : endblk_size_); }
    idx_t k_begin() const { return vertical().k_begin(); }
    idx_t k_end() const { return vertical().k_end(); }

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;

private:  // methods
    array::ArrayShape config_shape(const eckit::Configuration&) const;
    array::ArrayAlignment config_alignment(const eckit::Configuration&) const;
    array::ArraySpec config_spec(const eckit::Configuration&) const;

private:  // data
    idx_t nproma_;
    idx_t endblk_size_;
    idx_t nblks_;

    detail::StructuredColumns* structuredcolumns_;
    functionspace::StructuredColumns structuredcolumns_handle_;

    void setup(const eckit::Configuration& config);
};

// -------------------------------------------------------------------

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
