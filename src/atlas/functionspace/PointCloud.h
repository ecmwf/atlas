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

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

namespace atlas {
class Grid;

namespace functionspace {

//------------------------------------------------------------------------------------------------------

namespace detail {

class PointCloud : public functionspace::FunctionSpaceImpl {
public:
    template <typename Point>
    PointCloud(const std::vector<Point>&);
    PointCloud(const Field& lonlat);
    PointCloud(const Field& lonlat, const Field& ghost);
    PointCloud(const Grid&);
    virtual ~PointCloud() override {}
    virtual std::string type() const override { return "PointCloud"; }
    virtual operator bool() const override { return true; }
    virtual size_t footprint() const override { return sizeof(*this); }
    virtual std::string distribution() const override;
    Field lonlat() const override { return lonlat_; }
    const Field& vertical() const { return vertical_; }
    Field ghost() const override;
    virtual idx_t size() const override { return lonlat_.shape(0); }

    using FunctionSpaceImpl::createField;
    virtual Field createField(const eckit::Configuration&) const override;
    virtual Field createField(const Field&, const eckit::Configuration&) const override;

    void haloExchange(const FieldSet&, bool /*on_device*/ = false) const override {}
    void haloExchange(const Field&, bool /*on_device*/ = false) const override {}

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

private:
    Field lonlat_;
    Field vertical_;
    mutable Field ghost_;
    idx_t levels_{0};
};

//------------------------------------------------------------------------------------------------------

}  // namespace detail

//------------------------------------------------------------------------------------------------------

class PointCloud : public FunctionSpace {
public:
    PointCloud(const FunctionSpace&);
    PointCloud(const Field& points);
    PointCloud(const std::vector<PointXY>&);
    PointCloud(const std::vector<PointXYZ>&);
    PointCloud(const std::initializer_list<std::initializer_list<double>>&);
    PointCloud(const Grid& grid);

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    const Field& vertical() const { return functionspace_->vertical(); }

    detail::PointCloud::Iterate iterate() const { return functionspace_->iterate(); }


private:
    const detail::PointCloud* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas
