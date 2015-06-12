/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_h
#define atlas_h

#include <iosfwd>
#include <string>
#include <vector>

#include "eckit/container/DenseMap.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/Metadata.h"
#include "atlas/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class FunctionSpace;
class GridDistribution;

//------------------------------------------------------------------------------------------------------

class Mesh : public eckit::Owned {

public: // types

	typedef eckit::SharedPtr<Mesh> Ptr;

public: // methods

	static Mesh::Ptr create();

  /*!
   * @brief Construct a empty Mesh
   */
  Mesh();

  /*!
   * @brief Construct mesh from grid.
   * The mesh is global and only has a "nodes" FunctionSpace
   */
  Mesh(const Grid&);

  virtual ~Mesh();

        Metadata& metadata() { return metadata_; }
        const Metadata& metadata() const { return metadata_; }

  /// checks if function space exists
	bool has_function_space(const std::string& name) const;

	/// Takes ownership, and will be deleted automatically
	FunctionSpace& create_function_space(const std::string& name, const std::string& shape_func, const std::vector<int>& shape);

	/// accessor by name
	FunctionSpace& function_space(const std::string& name) const;

	/// accessor by index
	FunctionSpace& function_space( size_t ) const;

	/// number of functional spaces
	int nb_function_spaces() const { return function_spaces_.size(); }

	/// checks if has a Grid
	bool has_grid() const { return grid_; }

	/// assign a Grid to this Mesh
	void set_grid( const Grid& p ) { grid_ = &p; }

	/// accessor of the Grid
	const Grid& grid() const {  ASSERT( grid_ ); return *grid_; }

	friend std::ostream& operator<<(std::ostream&, const Mesh&);

  FunctionSpace& add_nodes(const Grid& g);

  FunctionSpace& add_nodes(size_t nb_nodes);

private: // members

	Metadata      metadata_;

	const Grid* grid_;

	typedef eckit::DenseMap< std::string, eckit::SharedPtr<FunctionSpace> > StoreFS_t;

	StoreFS_t function_spaces_;

};

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
	Mesh* atlas__Mesh__new ();
	void atlas__Mesh__delete (Mesh* This);
	void atlas__Mesh__create_function_space (Mesh* This, char* name,char* shape_func,int shape[], int shape_size);
	FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name);
	Grid* atlas__Mesh__grid (Mesh* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif // atlas_h
