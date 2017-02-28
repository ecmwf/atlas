/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <vector>
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/utils/MD5.h"
#include "eckit/value/Properties.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Grid : public eckit::Owned {

public:  // types

    using Config    = atlas::util::Config;
    using builder_t = eckit::BuilderT1<Grid>;
    using ARG1      = const Config&;
    using Ptr       = eckit::SharedPtr<Grid>;
    using Point     = PointLonLat; // must be sizeof(double)*2
    using uid_t     = std::string;

public:  // methods

    static std::string className();

    static Grid* create( const Config& );

    static Grid* create( const Grid::uid_t& );

    /// ctor (default)
    Grid();

    /// dtor
    virtual ~Grid();

    /// Human readable name (may not be unique)
    virtual std::string shortName() const = 0;
    virtual std::string gridType() const=0;

    /// Unique grid id
    /// Computed from the shortName and the hash
    uid_t uniqueId() const;

    /// Adds to the MD5 the information that makes this Grid unique
    virtual void hash(eckit::MD5&) const = 0;

    /// @returns the hash of the information that makes this Grid unique
    eckit::MD5::digest_t hash() const;

    /// @return area represented by the grid
    const Domain& domain() const { return domain_; }

    /// @return projection (mapping between geographic coordinates and grid coordinates)
    const Projection& projection() const { return projection_; }

    /// @return number of grid points
    /// @note This methods should have constant access time, if necessary derived
    //        classes should compute it at construction
    virtual size_t npts() const = 0;

    /// Fill provided parameter with grid points, as (lon,lat) values
    /// @post resizes the vector
    virtual void lonlat(std::vector<Point>&) const = 0;
    //virtual void xy(std::vector<PointXY>&) const = 0;

    /// Fills the provided vector with the (lon,lat) values
    /// @post resizes the vector
    void fillLonLat(std::vector<double>&) const;

    /// Fills the provided array with the (lon,lat) values
    /// @note Assumes that the input array has been allocated with correct size
    /// @param array is an array already allocated with enough size to store all the latlon values
    /// @param arraySize is the size of the array
    void fillLonLat(double array[], size_t arraySize) const;

    //virtual std::string gridType() const = 0;

    virtual std::string getOptimalMeshGenerator() const;

    virtual eckit::Properties spec() const;

    virtual bool same(const grid::Grid&) const;

protected:  // methods

    /// Fill provided memory buffer with the grid points, as (lon,lat) values
    /// This implementation in the base Grid class is not optimal as it incurs in double copy
    /// Derived classes should reimplement more optimised versions.
    ///
    /// @note Assumes that the input buffer has been allocated with correct size,
    ///       possibly from calling method npts()
    ///
    /// @param array to be filled in with the (lon,lat) values
    /// @param size number of doubles in array
    ///
    /// @return the size of bytes copyied in
    virtual size_t copyLonLatMemory(double* pts, size_t size) const;

    virtual void print(std::ostream&) const = 0;

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const grid::Grid& p) {
        p.print(s);
        return s;
    }

private:  // members

    /// Cache the unique ID
    mutable uid_t uid_;

    /// Cache the hash
    mutable eckit::MD5::digest_t hash_;

protected: // members

    Projection projection_;
    Domain     domain_;
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
