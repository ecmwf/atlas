/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Peter Bispham
/// @author Tiago Quintino
/// @date Oct 2013

#ifndef atlas_grid_Grid_H
#define atlas_grid_Grid_H

#include <cstddef>
#include <vector>
#include <cmath>

#include "eckit/memory/NonCopyable.h"
#include "eckit/exception/Exceptions.h"

#include "eckit/geometry/Point2.h"

#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Field.hpp"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/// Interface to a grid of points in a 2d cartesian space
/// For example a LatLon grid or a Reduced Graussian grid

class Grid : private eckit::NonCopyable {
public: // types

    typedef atlas::Mesh                         Mesh;      ///< mesh type
    typedef eckit::geometry::LLPoint2           Point;     ///< point type
    typedef eckit::geometry::BoundBox2<Point>   BoundBox;  ///< boundbox type

    typedef std::shared_ptr<Grid> Ptr;

//    class Iterator {
//    public:
//        virtual ~GridIterator() {}
//        virtual bool next( double& lat, double& lon ) = 0;
//    };
//    virtual Iterator* makeIterator() const = 0;

    class Coords {
    public:

        Coords( std::vector<double>& v ) : coords_(v) { ASSERT( v.size() && v.size()%2 == 0 ); }
        size_t size() const { return coords_.size() / 2; }
        double& lat( size_t i ) { return coords_[i];   }
        double& lon( size_t i ) { return coords_[i+1]; }

    private:

        std::vector<double>& coords_;
    };


public: // methods

    Grid();

    virtual ~Grid();

    virtual std::string hash() const = 0;

    virtual const char* gridType() const = 0;

    virtual BoundBox boundingBox() const = 0;

    virtual size_t nPoints() const = 0;

    virtual void coordinates( Grid::Coords& ) const = 0;

    /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
    virtual const std::vector<Point>& coordinates() const = 0;

//    const Mesh& mesh() const;
    Mesh& mesh();

protected:

    virtual void make_mesh();

    static int scanningMode(long iScansNegatively, long jScansPositively);


    std::unique_ptr< Mesh > mesh_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
