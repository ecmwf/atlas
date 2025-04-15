/*
 * (C) Copyright 2013 ECMWF
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include <memory>

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"
#include "atlas/grid/Partitioner.h"

namespace atlas {
namespace parallel {
class HaloExchange;
class GatherScatter;
}  // namespace parallel
}  // namespace atlas

namespace atlas {
class Grid;

namespace functionspace {

//------------------------------------------------------------------------------------------------------

namespace detail {

class PointCloud : public functionspace::FunctionSpaceImpl {
public:
    template <typename Point>
    PointCloud(const std::vector<Point>&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Field& lonlat, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Field& lonlat, const Field& ghost, const eckit::Configuration& = util::NoConfig());
    PointCloud(const FieldSet&, const eckit::Configuration& = util::NoConfig());  // assuming lonlat ghost ridx and partition present
    PointCloud(const Grid&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());
    ~PointCloud() override {}
    std::string type() const override { return "PointCloud"; }
    operator bool() const override { return true; }
    size_t footprint() const override { return sizeof(*this); }
    std::string distribution() const override;
    const Grid& grid() const override;
    Field lonlat() const override { return lonlat_; }
    const Field& vertical() const { return vertical_; }
    Field ghost() const override;
    Field remote_index() const override { return remote_index_; }
    Field global_index() const override { return global_index_; }
    Field partition() const override { return partition_; }
    idx_t size() const override { return lonlat_.shape(0); }
    idx_t part() const override { return part_; }
    idx_t nb_parts() const override { return nb_partitions_; }
    std::string mpi_comm() const override { return mpi_comm_; }

    using FunctionSpaceImpl::createField;
    Field createField(const eckit::Configuration&) const override;
    Field createField(const Field&, const eckit::Configuration&) const override;

    void haloExchange(const FieldSet&, bool on_device = false) const override;
    void haloExchange(const Field&, bool on_device = false) const override;

    void adjointHaloExchange(const FieldSet&, bool on_device = false) const override;
    void adjointHaloExchange(const Field&, bool on_device = false) const override;

    const parallel::HaloExchange& halo_exchange() const;

    void gather(const FieldSet&, FieldSet&) const override;
    void gather(const Field&, Field&) const override;
    const parallel::GatherScatter& gather() const override;

    void scatter(const FieldSet&, FieldSet&) const override;
    void scatter(const Field&, Field&) const override;
    const parallel::GatherScatter& scatter() const override;


    template <typename Point>
    class IteratorT {
    public:
        using difference_type   = std::ptrdiff_t;
        using value_type        = Point;
        using pointer           = Point*;
        using reference         = Point&;
        using iterator_category = std::output_iterator_tag;

        IteratorT(const PointCloud& fs, bool begin = true);

        bool next(Point&);

        const Point operator*() const;

        const IteratorT& operator++() {
            ++n_;
            return *this;
        }

        bool operator==(const IteratorT& other) const { return n_ == other.n_; }
        bool operator!=(const IteratorT& other) const { return n_ != other.n_; }

    private:
        const PointCloud& fs_;
        const array::ArrayView<const double, 2> xy_;
        const array::ArrayView<const double, 1> z_;
        idx_t n_;
        idx_t size_;
    };

    template <typename Point>
    class IterateT {
    public:
        using value_type     = Point;
        using iterator       = IteratorT<Point>;
        using const_iterator = iterator;

    public:
        IterateT(const PointCloud& fs): fs_(fs) {}
        iterator begin() const { return IteratorT<Point>(fs_); }
        iterator end() const { return IteratorT<Point>(fs_, false); }
        idx_t size() const { return fs_.size(); }

    private:
        const PointCloud& fs_;
    };


    class Iterate {
    public:
        Iterate(const PointCloud& fs): fs_(fs) {}
        IterateT<PointXYZ> xyz() const {
            ATLAS_ASSERT(fs_.vertical());
            return IterateT<PointXYZ>{fs_};
        }
        IterateT<PointXY> xy() const { return IterateT<PointXY>{fs_}; }
        IterateT<PointLonLat> lonlat() const { return IterateT<PointLonLat>{fs_}; }

    private:
        const PointCloud& fs_;
    };

    Iterate iterate() const { return Iterate(*this); }

private:
    array::ArrayShape config_shape(const eckit::Configuration& config) const;

    array::ArrayAlignment config_alignment(const eckit::Configuration& config) const;

    array::ArraySpec config_spec(const eckit::Configuration& config) const;

    array::DataType config_datatype(const eckit::Configuration& config) const;

    std::string config_name(const eckit::Configuration& config) const;

    void set_field_metadata(const eckit::Configuration& config, Field& field) const;

    void create_remote_index() const;

private:
    mutable Grid grid_;
    Field lonlat_;
    Field vertical_;
    mutable Field ghost_;
    mutable Field remote_index_;
    Field global_index_;
    Field partition_;
    idx_t size_owned_;
    idx_t size_global_{0};
    idx_t max_glb_idx_{0};

    idx_t levels_{0};
    idx_t part_{0};
    idx_t nb_partitions_{1};
    std::string mpi_comm_;

    mutable std::unique_ptr<parallel::HaloExchange> halo_exchange_;
    mutable std::unique_ptr<parallel::GatherScatter> gather_scatter_;

    void setupHaloExchange();
    void setupGatherScatter();

};

//------------------------------------------------------------------------------------------------------

}  // namespace detail

//------------------------------------------------------------------------------------------------------

class PointCloud : public FunctionSpace {
public:
    PointCloud(const FunctionSpace&);
    PointCloud(const Field& points, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Field&, const Field&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const FieldSet& flds, const eckit::Configuration& = util::NoConfig());
    PointCloud(const std::vector<PointXY>&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const std::vector<PointXYZ>&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const std::initializer_list<std::initializer_list<double>>&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Grid&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());
    PointCloud(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    const Field& vertical() const { return functionspace_->vertical(); }

    detail::PointCloud::Iterate iterate() const { return functionspace_->iterate(); }


private:
    const detail::PointCloud* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas
