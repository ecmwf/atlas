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

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "atlas/domain/Domain.h"
#include "atlas/library/config.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Object.h"

namespace eckit {
class Hash;
}
namespace atlas {
class PointXY;
class PointLonLat;
namespace util {
class Config;
};
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class GridObserver;

class Grid : public util::Object {
public:  // types
    using Projection = atlas::Projection;
    using Domain     = atlas::Domain;
    using Config     = atlas::util::Config;
    using Spec       = atlas::util::Config;
    using uid_t      = std::string;
    using hash_t     = std::string;

    template <typename Derived, typename Point>
    class IteratorT {
    public:
        using difference_type   = size_t;
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = Point;
        using pointer           = const Point*;
        using reference         = const Point&;

    public:
        virtual bool next(value_type&)                         = 0;
        virtual reference operator*() const                    = 0;
        virtual const Derived& operator++()                    = 0;
        virtual const Derived& operator+=(difference_type)     = 0;
        virtual bool operator==(const Derived&) const          = 0;
        virtual bool operator!=(const Derived&) const          = 0;
        virtual difference_type distance(const Derived&) const = 0;
        virtual std::unique_ptr<Derived> clone() const         = 0;
        virtual ~IteratorT() {}
    };
    class IteratorXY : public IteratorT<IteratorXY, PointXY> {};
    class IteratorLonLat : public IteratorT<IteratorLonLat, PointLonLat> {};

public:  // methods
    static const Grid* create(const Config&);

    static const Grid* create(const std::string& name);

    static const Grid* create(const std::string& name, const Config&);

    static const Grid* create(const Grid&, const Domain&);

    /// ctor (default)
    Grid();

    /// dtor
    virtual ~Grid();

    /// Human readable name (may not be unique)
    virtual std::string name() const = 0;
    virtual std::string type() const = 0;

    /// Unique grid id
    /// Computed from the hash. Can be used to compare 2 grids.
    virtual uid_t uid() const;

    /// Adds to the hash the information that makes this Grid unique
    virtual void hash(eckit::Hash&) const = 0;

    /// @returns the hash of the information that makes this Grid unique
    std::string hash() const;

    /// @return area represented by the grid
    const Domain& domain() const { return domain_; }

    /// @return parallel/meridian limits containing the grid
    virtual RectangularLonLatDomain lonlatBoundingBox() const = 0;

    virtual size_t footprint() const { return 0; }

    /// @return projection (mapping between geographic coordinates and grid
    /// coordinates)
    const Projection& projection() const { return projection_; }

    /// @return number of grid points
    /// @note This methods should have constant access time, if necessary derived
    //        classes should compute it at construction
    virtual idx_t size() const = 0;

    virtual Spec spec() const = 0;

    virtual std::unique_ptr<IteratorXY> xy_begin() const         = 0;
    virtual std::unique_ptr<IteratorXY> xy_end() const           = 0;
    virtual std::unique_ptr<IteratorLonLat> lonlat_begin() const = 0;
    virtual std::unique_ptr<IteratorLonLat> lonlat_end() const   = 0;

    void attachObserver(GridObserver&) const;
    void detachObserver(GridObserver&) const;

    virtual Config meshgenerator() const;
    virtual Config partitioner() const;

protected:  // methods
    /// Fill provided me
    virtual void print(std::ostream&) const = 0;

private:  // methods
    friend std::ostream& operator<<(std::ostream& s, const Grid& p) {
        p.print(s);
        return s;
    }

private:  // members
    /// Cache the unique ID
    mutable uid_t uid_;

    /// Cache the hash
    mutable hash_t hash_;

protected:  // members
    Projection projection_;
    Domain domain_;

    mutable std::vector<GridObserver*> grid_observers_;
};

class GridObserver {
private:
    std::vector<const Grid*> registered_grids_;

public:
    void registerGrid(const Grid& grid) {
        if (std::find(registered_grids_.begin(), registered_grids_.end(), &grid) == registered_grids_.end()) {
            registered_grids_.push_back(&grid);
            grid.attachObserver(*this);
        }
    }
    void unregisterGrid(const Grid& grid) {
        auto found = std::find(registered_grids_.begin(), registered_grids_.end(), &grid);
        if (found != registered_grids_.end()) {
            registered_grids_.erase(found);
            grid.detachObserver(*this);
        }
    }
    virtual ~GridObserver() {
        for (auto grid : registered_grids_) {
            grid->detachObserver(*this);
        }
    }
    virtual void onGridDestruction(Grid&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

extern "C" {
void atlas__grid__Grid__name(const Grid* This, char*& uid, int& size);
idx_t atlas__grid__Grid__size(Grid* This);
Grid::Spec* atlas__grid__Grid__spec(Grid* This);
void atlas__grid__Grid__uid(const Grid* This, char*& uid, int& size);
Grid::Domain::Implementation* atlas__grid__Grid__lonlat_bounding_box(const Grid* This);
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
