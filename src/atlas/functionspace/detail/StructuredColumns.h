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

#include "atlas/array/DataType.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
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
class HaloExchange;
class Checksum;
}  // namespace parallel
}  // namespace atlas

namespace atlas {
class Field;
class FieldSet;
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

class StructuredColumnsHaloExchangeCache;
class StructuredColumnsGatherScatterCache;
class StructuredColumnsChecksumCache;


// -------------------------------------------------------------------

class StructuredColumns : public FunctionSpaceImpl {
public:
    StructuredColumns(const Grid&, const eckit::Configuration& = util::NoConfig());

    StructuredColumns(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());

    StructuredColumns(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());

    StructuredColumns(const Grid&, const grid::Distribution&, const Vertical&,
                      const eckit::Configuration& = util::NoConfig());

    StructuredColumns(const Grid&, const Vertical&, const eckit::Configuration& = util::NoConfig());

    StructuredColumns(const Grid&, const Vertical&, const grid::Partitioner&,
                      const eckit::Configuration& = util::NoConfig());

    ~StructuredColumns() override;

    static std::string static_type() { return "StructuredColumns"; }
    std::string type() const override { return static_type(); }
    std::string distribution() const override;

    /// @brief Create a Structured field
    Field createField(const eckit::Configuration&) const override;

    Field createField(const Field&, const eckit::Configuration&) const override;

    void gather(const FieldSet&, FieldSet&) const override;
    void gather(const Field&, Field&) const override;

    void scatter(const FieldSet&, FieldSet&) const override;
    void scatter(const Field&, Field&) const override;

    void haloExchange(const FieldSet&, bool on_device = false) const override;
    void haloExchange(const Field&, bool on_device = false) const override;

    void adjointHaloExchange(const FieldSet&, bool on_device = false) const override;
    void adjointHaloExchange(const Field&, bool on_device = false) const override;

    idx_t sizeOwned() const { return size_owned_; }
    idx_t sizeHalo() const { return size_halo_; }
    idx_t size() const override { return size_halo_; }

    idx_t levels() const { return nb_levels_; }

    idx_t halo() const { return halo_; }

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;


    const Vertical& vertical() const { return vertical_; }

    const StructuredGrid& grid() const override;

    const Projection& projection() const override { return grid().projection(); }

    idx_t i_begin(idx_t j) const { return i_begin_[j]; }
    idx_t i_end(idx_t j) const { return i_end_[j]; }

    idx_t i_begin_halo(idx_t j) const { return i_begin_halo_[j]; }
    idx_t i_end_halo(idx_t j) const { return i_end_halo_[j]; }

    idx_t j_begin() const { return j_begin_; }
    idx_t j_end() const { return j_end_; }

    idx_t j_begin_halo() const { return j_begin_halo_; }
    idx_t j_end_halo() const { return j_end_halo_; }

    idx_t k_begin() const { return vertical_.k_begin(); }
    idx_t k_end() const { return vertical_.k_end(); }

    idx_t index(idx_t i, idx_t j) const {
        check_bounds(i, j);
        return ij2gp_(i, j);
    }

    Field lonlat() const override { return field_xy_; }
    Field xy() const { return field_xy_; }
    Field z() const { return field_z_; }
    Field partition() const override { return field_partition_; }
    Field global_index() const override { return field_global_index_; }
    Field remote_index() const override {
        if (not field_remote_index_) {
            create_remote_index();
        }
        return field_remote_index_;
    }
    Field index_i() const { return field_index_i_; }
    Field index_j() const { return field_index_j_; }
    Field ghost() const override { return field_ghost_; }

    void compute_xy(idx_t i, idx_t j, PointXY& xy) const;
    PointXY compute_xy(idx_t i, idx_t j) const {
        PointXY xy;
        compute_xy(i, j, xy);
        return xy;
    }

    virtual size_t footprint() const override;

    const util::PartitionPolygon& polygon(idx_t halo = 0) const override;

    const util::PartitionPolygons& polygons() const override;

    idx_t part() const override { return part_; }
    idx_t nb_parts() const override { return nb_partitions_; }
    std::string mpi_comm() const override { return mpi_comm_; }

private:  // methods
    idx_t config_size(const eckit::Configuration& config) const;
    array::DataType config_datatype(const eckit::Configuration&) const;
    std::string config_name(const eckit::Configuration&) const;
    idx_t config_levels(const eckit::Configuration&) const;
    array::ArrayShape config_shape(const eckit::Configuration&) const;
    array::ArrayAlignment config_alignment(const eckit::Configuration&) const;
    array::ArraySpec config_spec(const eckit::Configuration&) const;
    void set_field_metadata(const eckit::Configuration&, Field&) const;

    void check_bounds(idx_t i, idx_t j) const {
#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
        if (j < j_begin_halo() || j >= j_end_halo()) {
            throw_outofbounds(i, j);
        }
        if (i < i_begin_halo(j) || i >= i_end_halo(j)) {
            throw_outofbounds(i, j);
        }
#endif
    }
    [[noreturn]] void throw_outofbounds(idx_t i, idx_t j) const;

    const parallel::GatherScatter& gather() const override;
    const parallel::GatherScatter& scatter() const override;
    const parallel::Checksum& checksum() const;
    const parallel::HaloExchange& halo_exchange() const;

    void create_remote_index() const;

private:  // data
    std::string distribution_;

    const Vertical vertical_;
    idx_t nb_levels_;

    idx_t size_owned_;
    idx_t size_halo_;
    idx_t halo_;

    friend class StructuredColumnsHaloExchangeCache;
    friend class StructuredColumnsGatherScatterCache;
    friend class StructuredColumnsChecksumCache;
    bool periodic_points_{false};

    const StructuredGrid* grid_;
    mutable util::ObjectHandle<parallel::GatherScatter> gather_scatter_;
    mutable util::ObjectHandle<parallel::Checksum> checksum_;
    mutable util::ObjectHandle<parallel::HaloExchange> halo_exchange_;
    mutable std::vector<util::ObjectHandle<util::PartitionPolygon>> polygons_;
    mutable util::PartitionPolygons all_polygons_;

    Field field_xy_;
    Field field_z_;
    Field field_partition_;
    Field field_global_index_;
    mutable Field field_remote_index_;
    Field field_index_i_;
    Field field_index_j_;
    Field field_ghost_;

    class Map2to1 {
    public:
        Map2to1() { resize({1, 0}, {1, 0}); }

        Map2to1(std::array<idx_t, 2> i_range, std::array<idx_t, 2> j_range) { resize(i_range, j_range); }

        void resize(std::array<idx_t, 2> i_range, std::array<idx_t, 2> j_range);

        atlas::vector<idx_t> data_;
        idx_t i_min_;
        idx_t i_max_;
        idx_t j_min_;
        idx_t j_max_;
        idx_t j_stride_;

        idx_t operator()(idx_t i, idx_t j) const { return data_[(i - i_min_) + (j - j_min_) * j_stride_] - 1; }

        void set(idx_t i, idx_t j, idx_t n) { data_[(i - i_min_) + (j - j_min_) * j_stride_] = n + 1; }

        static idx_t missing() { return std::numeric_limits<idx_t>::max() - 1; }

        size_t footprint() const;

    private:
        void print(std::ostream&) const;

        friend std::ostream& operator<<(std::ostream& s, const Map2to1& p) {
            p.print(s);
            return s;
        }
    };

    class IndexRange {
    public:
        IndexRange() { resize(1, 0); }

        IndexRange(idx_t min, idx_t max) { resize(min, max); }

        std::vector<idx_t> data_;
        idx_t min_;
        idx_t max_;

        idx_t operator()(idx_t i) const { return data_[i - min_]; }

        idx_t& operator()(idx_t i) { return data_[i - min_]; }

        idx_t operator[](idx_t i) const { return data_[i - min_]; }

        idx_t& operator[](idx_t i) { return data_[i - min_]; }

        idx_t missing() const { return std::numeric_limits<idx_t>::max() - 1; }

        idx_t size() const { return idx_t(data_.size()); }

        void resize(idx_t min, idx_t max) {
            min_ = min;
            max_ = max;
            data_.resize(max_ - min_ + 1, missing() + 1);
        }

    private:
        void print(std::ostream&) const;

        friend std::ostream& operator<<(std::ostream& s, const IndexRange& p) {
            p.print(s);
            return s;
        }
    };

    idx_t j_begin_;
    idx_t j_end_;
    std::vector<idx_t> i_begin_;
    std::vector<idx_t> i_end_;
    idx_t j_begin_halo_;
    idx_t j_end_halo_;
    IndexRange i_begin_halo_;
    IndexRange i_end_halo_;

    idx_t north_pole_included_;
    idx_t south_pole_included_;
    idx_t ny_;

    idx_t part_;
    idx_t nb_partitions_;
    std::string mpi_comm_;

    friend struct StructuredColumnsFortranAccess;
    friend struct BlockStructuredColumnsFortranAccess;
    Map2to1 ij2gp_;

    friend class BlockStructuredColumns;
    void setup(const grid::Distribution& distribution, const eckit::Configuration& config);

    friend class BlockStructuredColumns;
};

// -------------------------------------------------------------------

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
