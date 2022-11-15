/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/BlockStructuredColumns.h"

#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>

#include "eckit/utils/MD5.h"

#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/domain.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/StructuredPartitionPolygon.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/fill.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/detail/Cache.h"

namespace {
using namespace atlas;

template <class ValueType>
void block_copy(const Field sloc, Field loc, const functionspace::detail::BlockStructuredColumns& fs) {
    if (sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                        loc_v(jblk, jvar, jlev, jrof) = sloc_v(blk_begin+jrof, jlev, jvar);
                    }
                }
            }
        }
    }
    else if (not sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    loc_v(jblk, jlev, jrof) = sloc_v(blk_begin+jrof, jlev);
                }
            }
        }
    }
    else if (sloc.variables() and not sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    loc_v(jblk, jvar, jrof) = sloc_v(blk_begin+jrof, jvar);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                loc_v(jblk, jrof) = sloc_v(blk_begin+jrof);
            }
        }
    }
}

template <class ValueType>
void rev_block_copy(const Field loc, Field sloc, const functionspace::detail::BlockStructuredColumns& fs) {
    if (loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                        sloc_v(blk_begin+jrof, jlev, jvar) = loc_v(jblk, jvar, jlev, jrof);
                    }
                }
            }
        }
    }
    else if (not loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    sloc_v(blk_begin+jrof, jlev) = loc_v(jblk, jlev, jrof);
                }
            }
        }
    }
    else if (loc.variables() and not loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    sloc_v(blk_begin+jrof, jvar) = loc_v(jblk, jvar, jrof);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                sloc_v(blk_begin+jrof) = loc_v(jblk, jrof);
            }
        }
    }
}

}// namespace

namespace atlas {
namespace functionspace {
namespace detail {


array::ArrayShape BlockStructuredColumns::config_shape(const eckit::Configuration& config) const {
    array::ArrayShape shape;

    bool global = false;
    config.get("global", global);
    if (global) {
        return structuredcolumns_->config_shape(config);
    }
    else {
        shape.emplace_back(nblks_);
        idx_t variables(0);
        config.get("variables", variables);
        if (variables > 0) {
            shape.emplace_back(variables);
        }
        idx_t levels(structuredcolumns_->levels());
        config.get("levels", levels);
        if (levels > 0) {
            shape.emplace_back(levels);
        }
        shape.emplace_back(nproma_);
    }
    return shape;
}

array::ArraySpec BlockStructuredColumns::config_spec(const eckit::Configuration& config) const {
    return array::ArraySpec(config_shape(config), structuredcolumns_->config_alignment(config));
}

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const eckit::Configuration& config):
    BlockStructuredColumns::BlockStructuredColumns(grid, grid::Partitioner(), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Partitioner& p, const eckit::Configuration& config):
    BlockStructuredColumns(grid, Vertical(config), p, config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution,
                                     const eckit::Configuration& config):
    BlockStructuredColumns(grid, distribution, Vertical(config), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution, const Vertical& vertical,
                                     const eckit::Configuration& config):
    structuredcolumns_(new StructuredColumns(grid, distribution, vertical, config)),
    structuredcolumns_handle_(structuredcolumns_){
    setup(config);
}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const eckit::Configuration& config):
    BlockStructuredColumns(grid, vertical, grid::Partitioner(), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const grid::Partitioner& p,
                                     const eckit::Configuration& config):
    structuredcolumns_(new StructuredColumns(grid, vertical, p, config)),
    structuredcolumns_handle_(structuredcolumns_){
    setup(config);
}

// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
Field BlockStructuredColumns::createField(const eckit::Configuration& options) const {
    Field field(structuredcolumns_->config_name(options), structuredcolumns_->config_datatype(options), config_spec(options));
    structuredcolumns_->set_field_metadata(options, field);
    field.set_functionspace(this);
    return field;
}

Field BlockStructuredColumns::createField(const Field& other, const eckit::Configuration& config) const {
    return createField(option::datatype(other.datatype()) | option::levels(other.levels()) |
                       option::variables(other.variables()) |
                       option::type(other.metadata().getString("type", "scalar")) | config);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void BlockStructuredColumns::scatter(const FieldSet& global_fieldset, FieldSet& local_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());
    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& glb      = global_fieldset[f];
        auto config = option::datatype(glb.datatype()) | option::levels(glb.levels()) | option::variables(glb.variables());
        auto sloc             = structuredcolumns_->createField(config);
        const idx_t nb_fields = 1;
        idx_t root(0);
        glb.metadata().get("owner", root);

        Field& loc            = local_fieldset[f];
        glb.metadata().broadcast(loc.metadata(), root);
        loc.metadata().set("global", false);

        if (sloc.datatype().kind() == array::DataType::kind<int>()) {
            parallel::Field<int const> glb_field(structuredcolumns_->make_leveled_view<const int>(glb));
            parallel::Field<int> sloc_field(structuredcolumns_->make_leveled_view<int>(sloc));
            structuredcolumns_->scatter().scatter(&glb_field, &sloc_field, nb_fields, root);
            block_copy<int>(sloc, loc, *this);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<long>()) {
            parallel::Field<long const> glb_field(structuredcolumns_->make_leveled_view<const long>(glb));
            parallel::Field<long> sloc_field(structuredcolumns_->make_leveled_view<long>(sloc));
            structuredcolumns_->scatter().scatter(&glb_field, &sloc_field, nb_fields, root);
            block_copy<long>(sloc, loc, *this);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<float>()) {
            parallel::Field<float const> glb_field(structuredcolumns_->make_leveled_view<const float>(glb));
            parallel::Field<float> sloc_field(structuredcolumns_->make_leveled_view<float>(sloc));
            structuredcolumns_->scatter().scatter(&glb_field, &sloc_field, nb_fields, root);
            block_copy<float>(sloc, loc, *this);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<double>()) {
            parallel::Field<double const> glb_field(structuredcolumns_->make_leveled_view<const double>(glb));
            parallel::Field<double> sloc_field(structuredcolumns_->make_leveled_view<double>(sloc));
            structuredcolumns_->scatter().scatter(&glb_field, &sloc_field, nb_fields, root);
            block_copy<double>(sloc, loc, *this);
        }
        else {
            throw_Exception("datatype not supported", Here());
        }
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void BlockStructuredColumns::scatter(const Field& global, Field& local) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void BlockStructuredColumns::gather(const FieldSet& local_fieldset, FieldSet& global_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());
    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& loc      = local_fieldset[f];
        auto config = option::datatype(loc.datatype()) | option::levels(loc.levels()) | option::variables(loc.variables());
        auto sloc             = structuredcolumns_->createField(config);


        idx_t root(0);
        sloc.metadata().set("global", false);
        Field& glb            = global_fieldset[f];
        glb.metadata().broadcast(sloc.metadata(), root);
        const idx_t nb_fields = 1;
        glb.metadata().get("owner", root);

        if (sloc.datatype().kind() == array::DataType::kind<int>()) {
            rev_block_copy<int>(loc, sloc, *this);
            parallel::Field<int const> sloc_field(structuredcolumns_->make_leveled_view<int>(sloc));
            parallel::Field<int> glb_field(structuredcolumns_->make_leveled_view<int>(glb));
            structuredcolumns_->gather().gather(&sloc_field, &glb_field, nb_fields, root);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<long>()) {
            rev_block_copy<long>(loc, sloc, *this);
            parallel::Field<long const> sloc_field(structuredcolumns_->make_leveled_view<long>(sloc));
            parallel::Field<long> glb_field(structuredcolumns_->make_leveled_view<long>(glb));
            structuredcolumns_->gather().gather(&sloc_field, &glb_field, nb_fields, root);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<float>()) {
            rev_block_copy<float>(loc, sloc, *this);
            parallel::Field<float const> sloc_field(structuredcolumns_->make_leveled_view<float>(sloc));
            parallel::Field<float> glb_field(structuredcolumns_->make_leveled_view<float>(glb));
            structuredcolumns_->gather().gather(&sloc_field, &glb_field, nb_fields, root);
        }
        else if (sloc.datatype().kind() == array::DataType::kind<double>()) {
            rev_block_copy<double>(loc, sloc, *this);
            parallel::Field<double const> sloc_field(structuredcolumns_->make_leveled_view<double>(sloc));
            parallel::Field<double> glb_field(structuredcolumns_->make_leveled_view<double>(glb));
            structuredcolumns_->gather().gather(&sloc_field, &glb_field, nb_fields, root);
        }
        else {
            throw_Exception("datatype not supported", Here());
        }
    }
}
// ----------------------------------------------------------------------------

void BlockStructuredColumns::setup(const eckit::Configuration &config) {
    nproma_ = 1;
    idx_t tmp_nproma;
    if (config.get("nproma", tmp_nproma)) {
        ATLAS_ASSERT(tmp_nproma > 0);
        nproma_ = tmp_nproma;
    }

    nblks_ = std::floor( structuredcolumns_->size() / nproma_ );
    endblk_size_ = nproma_;
    if (structuredcolumns_->size() % nproma_ > 0) {
        endblk_size_ = structuredcolumns_->size() - nblks_ * nproma_;
        nblks_++;
    }
}

// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void BlockStructuredColumns::gather(const Field& local, Field& global) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields, global_fields);
}

std::string BlockStructuredColumns::checksum(const Field&) const {
    ATLAS_NOTIMPLEMENTED;
}

std::string BlockStructuredColumns::checksum(const FieldSet&) const {
    ATLAS_NOTIMPLEMENTED;
}


// ----------------------------------------------------------------------------


}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
