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

    using Projection = atlas::grid::Projection;
    using Domain     = atlas::grid::Domain;
    using Config     = atlas::util::Config;
    using Spec       = eckit::Properties;
    using builder_t  = eckit::BuilderT1<Grid>;
    using ARG1       = const Config&;
    using Point      = PointXY; // must be sizeof(double)*2
    using uid_t      = std::string;


    class Iterator {
    public:
      virtual bool next(PointXY&) =0;
      virtual const PointXY operator *() const =0;
      virtual const Iterator& operator ++() =0;
      virtual bool operator ==(const Iterator &other) const =0;
      virtual bool operator !=(const Iterator &other) const =0;
    };


public:  // methods

    static std::string className();

    static const Grid* create( const Config& );

    static const Grid* create( const std::string& name, const Config& = Config() );

    /// ctor (default)
    Grid();

    /// dtor
    virtual ~Grid();

    /// Human readable name (may not be unique)
    virtual std::string name() const = 0;
    virtual std::string type() const=0;

    /// Unique grid id
    /// Computed from the hash. Can be used to compare 2 grids.
    uid_t uid() const;

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

    virtual std::string getOptimalMeshGenerator() const;

    virtual Spec spec() const;

    virtual bool same(const grid::Grid&) const;

    virtual Iterator* begin() const=0;
    virtual Iterator* end()   const=0;

protected:  // methods

    /// Fill provided me
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
    atlas::grid::Domain     domain_;
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
