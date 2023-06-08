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

#include <functional>
#include <type_traits>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/detail/BlockStructuredColumns.h"

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class BlockStructuredColumns : public FunctionSpace {
public:
    class Block;
public:
    BlockStructuredColumns();
    BlockStructuredColumns(const FunctionSpace&);
    BlockStructuredColumns(const Grid&, const eckit::Configuration& = util::NoConfig());
    BlockStructuredColumns(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());
    BlockStructuredColumns(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());
    BlockStructuredColumns(const Grid&, const Vertical&, const eckit::Configuration& = util::NoConfig());
    BlockStructuredColumns(const Grid&, const Vertical&, const grid::Partitioner&,
                      const eckit::Configuration& = util::NoConfig());
    BlockStructuredColumns(const Grid&, const grid::Distribution&, const Vertical&,
                      const eckit::Configuration& = util::NoConfig());

    static std::string type() { return detail::BlockStructuredColumns::static_type(); }

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    idx_t size() const { return functionspace_->size(); }
    idx_t levels() const { return functionspace_->levels(); }

    const Vertical& vertical() const { return functionspace_->vertical(); }

    const StructuredGrid& grid() const { return functionspace_->grid(); }

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;

    idx_t index(idx_t blk, idx_t rof) const { return functionspace_->index(blk, rof); }
    idx_t k_begin() const { return functionspace_->k_begin(); }
    idx_t k_end() const { return functionspace_->k_end(); }
    idx_t nproma() const { return functionspace_->nproma(); }
    idx_t nblks() const { return functionspace_->nblks(); }

    Field xy() const { return functionspace_->xy(); }
    Field partition() const { return functionspace_->partition(); }
    Field global_index() const { return functionspace_->global_index(); }
    Field remote_index() const { return functionspace_->remote_index(); }
    Field index_i() const { return functionspace_->index_i(); }
    Field index_j() const { return functionspace_->index_j(); }
    Field ghost() const { return functionspace_->ghost(); }

    const Block block(idx_t jblk) const {
        return Block(functionspace_->block_begin(jblk), functionspace_->block_size(jblk));
    }

    size_t footprint() const { return functionspace_->footprint(); }

    class Block {
    public:
        Block(idx_t begin, idx_t size) : begin_(begin), size_(size) {}
        idx_t index(idx_t j) const { return begin_ + j; };
        idx_t size() const { return size_; }
    private:
        idx_t begin_;
        idx_t size_;
    };

private:
    const detail::BlockStructuredColumns* functionspace_;
    void setup(const Grid& grid, const Vertical& vertical, const grid::Distribution& distribution,
               const eckit::Configuration& config);
};

// -------------------------------------------------------------------


}  // namespace functionspace
}  // namespace atlas
