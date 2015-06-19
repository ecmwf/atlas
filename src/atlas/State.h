/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date June 2015

#ifndef atlas_State_H
#define atlas_State_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/container/DenseMap.h"

namespace eckit { class Parametrisation; }
namespace atlas { class Field; }
namespace atlas { class Mesh; }
namespace atlas { class Grid; }
namespace atlas { class Metadata; }


namespace atlas {

/**
 * \brief State class as ultimate data holder class
 * Fields are owned by a State,
 * whereas Mesh or Grid can be shared between different States
 * Multiple meshes or grids can be added for e.g coupling models
 * or multigrid methods.
 */
class State : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr< State > Ptr;
  typedef Metadata Parameters;

private:

    /// @TODO Please consider **not** using DenseMap anymore
    ///       This class needs serious (re)considering.
    ///       ATM cannot safely remove entries from the collecition
    ///       I have a solution but don't know if it is good one.
    ///       DenseMap should probably not exist in eckit but in atlas

  typedef eckit::DenseMap< std::string, eckit::SharedPtr<Field> >  FieldMap;
  typedef eckit::DenseMap< std::string, eckit::SharedPtr<Grid>  >  GridMap;
  typedef eckit::DenseMap< std::string, eckit::SharedPtr<Mesh>  >  MeshMap;

public: // methods

//-- Constructors

  State();

  State( const eckit::Parametrisation& );

//-- Accessors

  const Field& field(const std::string& name) const;
        Field& field(const std::string& name);
  bool has_field(const std::string& name) const { return fields_.has(name); }
  std::vector< std::string > field_names() const;

  const Field& field(const size_t idx) const;
        Field& field(const size_t idx);
  size_t nb_fields() const { return fields_.size(); }

  const Mesh& mesh(const std::string& name = "") const;
        Mesh& mesh(const std::string& name = "");
  bool has_mesh(const std::string& name) const { return meshes_.has(name); }

  const Mesh& mesh(const size_t idx = 0) const;
        Mesh& mesh(const size_t idx = 0);
  size_t nb_meshes() const { return meshes_.size(); }

  const Grid& grid(const std::string& name = "") const;
        Grid& grid(const std::string& name = "");
  bool has_grid(const std::string& name) const { return grids_.has(name); }

  const Grid& grid(const size_t idx = 0) const;
        Grid& grid(const size_t idx = 0);
  size_t nb_grids() const { return grids_.size(); }

// -- Modifiers

  Field& add( Field* ); // Take shared ownership!
  Mesh&  add( Mesh*  ); // Take shared ownership!
  Grid&  add( Grid*  ); // Take shared ownership!

//  TODO: needs DenseMap to implement erase()
//  void remove_field(const std::string& name);
//  void remove_mesh(const std::string& name = "");
//  void remove_grid(const std::string& name = "");

private:

  void set_name( Field& );

  FieldMap fields_;
  MeshMap meshes_;
  GridMap grids_;

};

// ------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
}


} // namespace atlas


#endif
