/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Unstructured_H
#define atlas_Unstructured_H

/// @author Tiago Quintino
/// @date April 2013

#include <cstddef>
#include <vector>

#include "eckit/memory/ScopedPtr.h"

#include "atlas/Grid.h"


//-----------------------------------------------------------------------------

namespace atlas {



//-----------------------------------------------------------------------------

class Unstructured : public Grid {

public: // methods

	static std::string className() { return "atlas.grid.Unstructured"; }

	Unstructured( const eckit::Params& p );

    /// @warning temporary constructor taking a list of points
    Unstructured( std::vector< Point >* pts, const std::string& hash );

    virtual ~Unstructured();

	virtual std::string uid() const;
	virtual std::string hash() const;

    virtual BoundBox bounding_box() const;

    virtual size_t npts() const;

	virtual void coordinates( std::vector<double>& ) const;
	virtual void coordinates( std::vector<Point>& ) const;

    virtual std::string grid_type() const { return std::string("unstructured"); }

	virtual GridSpec spec() const;

	virtual bool same(const Grid&) const;

protected:

    eckit::ScopedPtr< std::vector< Point > > points_; ///< storage of coordinate points

    BoundBox bound_box_;              ///< bounding box for the domain

    std::string hash_;

};

//-----------------------------------------------------------------------------


} // namespace eckit

#endif
