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

#ifndef atlas_Grid_H
#define atlas_Grid_H

#include <cstddef>
#include <vector>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/value/Params.h"
#include "eckit/memory/Builder.h"

#include "eckit/geometry/Point2.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class Mesh;



class GridSpec;

//------------------------------------------------------------------------------------------------------

/// Interface to a grid of points in a 2d cartesian space
/// For example a LatLon grid or a Reduced Graussian grid
///
///      DECODE                       ATLAS                      ENCODE
///      NetCDFBuilder ---->|-------|         |----------|------>NetCDFWrite
///                         | Grid  |<------> | GridSpec |
///      GribBuilder ------>|-------|         |----------|------>GribWrite

class Grid : public eckit::Owned {

public: // types

	typedef eckit::BuilderT1<Grid> builder_t;
	typedef const eckit::Params&   ARG1;

    typedef eckit::geometry::LLPoint2           Point;     ///< point type
	typedef eckit::geometry::LLBoundBox2        BoundBox;  ///< bounding box type

    typedef eckit::SharedPtr<Grid> Ptr;

public: // methods

    static std::string className() { return "atlas.grid.Grid"; }

    static double degrees_eps();

    static Grid::Ptr create( const eckit::Params& );
    static Grid::Ptr create( const GridSpec& );

    Grid();

    virtual ~Grid();

	Ptr self() { return Ptr(this); }

    virtual std::string uid() const = 0;
    virtual std::string hash() const = 0;

    /// Assumes north > south, and east > west.
    /// and hence independent of scanning mode(since that is GRIB specific)
    /// Assumes lat -90 --> +90
    /// Assumes lon   0 --> +360
    /// When the bounding box is not on the grid, the co-ordinate values, will be
    /// on the enclosing grid points
    virtual BoundBox boundingBox() const = 0;

	/// Returns the number of points
	/// This methods should have constant access time
	/// If necessary derived classes should compute it at cosntruction
	virtual size_t nPoints() const = 0;

    /// Assumes we start at NORTH,WEST --> SOUTH,EAST
    /// Assumes that the input vectors have the correct size.
    /// Points represent latitude and longitude values
    virtual void coordinates( std::vector<double>& ) const = 0;
    virtual void coordinates( std::vector<Point>& ) const = 0;

    virtual std::string gridType() const = 0;

    virtual GridSpec spec() const = 0;

    virtual bool same(const Grid&) const = 0;

    Mesh& mesh();
    const Mesh& mesh() const;

protected: // methods

	/// helper function to initialize global grids, with a global area (BoundBox)
	static BoundBox makeGlobalBBox();

	/// helper function to create bounding boxes (for non-global grids)
	static BoundBox makeBBox( const eckit::Params& );

private: // methods

	void buildMesh() const;

private: // members

    mutable eckit::SharedPtr< Mesh > mesh_;

};

//------------------------------------------------------------------------------------------------------


} // namespace atlas

#endif
