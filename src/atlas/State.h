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

namespace eckit { class Parametrisation; }
namespace atlas { class Parametrisation; }
namespace atlas { class Field; }
namespace atlas { class Mesh; }
namespace atlas { class Grid; }


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
  typedef atlas::Parametrisation Parameters;

public: // methods

//-- Static methods

  static State* create(const std::string& factory, const eckit::Parametrisation& = Parameters() );

//-- Constructors

  State();

//-- Accessors

  const Field& field(const std::string& name) const;
        Field& field(const std::string& name);
  bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }
  std::vector< std::string > field_names() const;

  const Field& field(const size_t idx) const;
        Field& field(const size_t idx);
  size_t nb_fields() const { return fields_.size(); }

  const Mesh& mesh(const std::string& name = "") const;
        Mesh& mesh(const std::string& name = "");
  bool has_mesh(const std::string& name) const { return (meshes_.find(name) != meshes_.end()); }

  const Mesh& mesh(const size_t idx = 0) const;
        Mesh& mesh(const size_t idx = 0);
  size_t nb_meshes() const { return meshes_.size(); }

  const Grid& grid(const std::string& name = "") const;
        Grid& grid(const std::string& name = "");
  bool has_grid(const std::string& name) const { return (grids_.find(name) != grids_.end()); }

  const Grid& grid(const size_t idx = 0) const;
        Grid& grid(const size_t idx = 0);
  size_t nb_grids() const { return grids_.size(); }

// -- Modifiers

  Field& add( Field* ); // Take shared ownership!
  Mesh&  add( Mesh*  ); // Take shared ownership!
  Grid&  add( Grid*  ); // Take shared ownership!

  void remove_field(const std::string& name);
  void remove_mesh(const std::string& name = "");
  void remove_grid(const std::string& name = "");

private:

  typedef std::map< std::string, eckit::SharedPtr<Field> >  FieldMap;
  typedef std::map< std::string, eckit::SharedPtr<Grid>  >  GridMap;
  typedef std::map< std::string, eckit::SharedPtr<Mesh>  >  MeshMap;

private:

  FieldMap fields_;
  MeshMap meshes_;
  GridMap grids_;

};

//------------------------------------------------------------------------------------------------------

class StateFactory {
  public:

    /*!
     * \brief build StateCreator with options specified in parametrisation
     * \return mesh generator
     */
    static State* build(const std::string& state_type);
    static State* build(const std::string& state_type, const eckit::Parametrisation&);

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);

  private:
    std::string name_;
    virtual State* make() = 0 ;
    virtual State* make(const eckit::Parametrisation&) = 0 ;

  protected:

    StateFactory(const std::string&);
    virtual ~StateFactory();
};


template<class T>
class StateBuilder : public StateFactory {
  virtual State* make() {
      return new T();
  }
  virtual State* make(const eckit::Parametrisation& param) {
        return new T(param);
  }
  public:
    StateBuilder(const std::string& name) : StateFactory(name) {}
};

// ------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define Parametrisation eckit::Parametrisation
extern "C"
{
  State* atlas__State__new ();
  State* atlas__State__create (const char* factory, const Parametrisation* params);
  void atlas__State__delete (State* This);
  void atlas__State__add_field (State* This, Field* field);
  void atlas__State__remove_field (State* This, const char* name);
  Field* atlas__State__field_by_name (State* This, const char* name);
  Field* atlas__State__field_by_index (State* This, int index);
  int atlas__State__nb_fields(const State* This);
  void atlas__State__add_grid (State* This, Grid* grid);
  void atlas__State__remove_grid (State* This, const char* name);
  Grid* atlas__State__grid_by_name (State* This, const char* name);
  Grid* atlas__State__grid_by_index (State* This, int index);
  int atlas__State__nb_grids(const State* This);
  void atlas__State__add_mesh (State* This, Mesh* grid);
  void atlas__State__remove_mesh (State* This, const char* name);
  Mesh* atlas__State__mesh_by_name (State* This, const char* name);
  Mesh* atlas__State__mesh_by_index (State* This, int index);
  int atlas__State__nb_meshes(const State* This);
}
#undef Parametrisation

} // namespace atlas


#endif
