/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Mesh_h
#define atlas_Mesh_h

#include <map>
#include <vector>
#include <string>

#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/Owned.h"

#include "atlas/mesh/Metadata.h"
#include "atlas/grid/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class FunctionSpace;

//------------------------------------------------------------------------------------------------------

class Mesh : public eckit::Owned {

public: // types

	typedef grid::Grid Grid;
	typedef eckit::SharedPtr<Mesh> Ptr;

public: // methods

//	static Mesh& create( const eckit::Params& );

	Mesh();

	virtual ~Mesh();

	Metadata& metadata() { return metadata_; }

	/// checks if function space exists
	bool has_function_space(const std::string& name) const;

	/// Takes ownership, and will be deleted automatically
	FunctionSpace& add_function_space( FunctionSpace* function_space );

	/// accessor by name
	FunctionSpace& function_space(const std::string& name) const;

	/// accessor by index
	FunctionSpace& function_space(int idx) const;

	/// number of functional spaces
	int nb_function_spaces() const { return function_spaces_.size(); }

	/// checks if has a Grid
	bool has_grid() const { return grid_; }

	/// assign a Grid to this Mesh
	void grid( grid::Grid& p )
	{
		DEBUG_VAR(grid_);
		grid_ = &p;
		DEBUG_VAR(grid_);
	}

	/// accessor of the Grid
	const grid::Grid& grid() const {  ASSERT( grid_ ); return *grid_; }

	/// accessor of the Grid
	grid::Grid& grid() { ASSERT( grid_ ); return *grid_; }

private: // members

	Metadata      metadata_;

	grid::Grid* grid_;

	std::map< std::string, size_t > index_; ///< index of function spaces

	std::vector< FunctionSpace* > function_spaces_; ///< function spaces

};

//------------------------------------------------------------------------------------------------------

typedef grid::Grid Grid;

// C wrapper interfaces to C++ routines
extern "C"
{
	Mesh* atlas__Mesh__new ();
	void atlas__Mesh__delete (Mesh* This);
	void atlas__Mesh__add_function_space (Mesh* This, FunctionSpace* function_space);
	FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name);
	Grid* atlas__Mesh__grid (Mesh* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_Mesh_h
