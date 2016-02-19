/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/util/Config.h"
#include "atlas/util/Metadata.h"

namespace eckit { class Parametrisation; }
namespace atlas { class Field; }

namespace atlas {

/**
 * \brief State class that owns a collection of fields
 */
class State : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr< State > Ptr;

public: // methods

//-- Constructors

  State();

  State(const std::string& generator, const eckit::Parametrisation& = Config());

//-- Accessors

  const Field& field(const std::string& name) const;
        Field& field(const std::string& name);
  bool has(const std::string& name) const { return (fields_.find(name) != fields_.end()); }
  std::vector< std::string > field_names() const;

  const Field& field(const size_t idx) const;
        Field& field(const size_t idx);
  size_t size() const { return fields_.size(); }

  const Field& operator[](const size_t idx) const { return field(idx); }
        Field& operator[](const size_t idx)       { return field(idx); }

  const Field& operator[](const std::string& name) const { return field(name); }
        Field& operator[](const std::string& name)       { return field(name); }

  const Metadata& metadata() const;
        Metadata& metadata();

// -- Modifiers

  void initialize(const std::string& generator, const eckit::Parametrisation& = Config() );

  Field& add( const eckit::SharedPtr<Field>& );
  Field& add( Field* ); // Take shared ownership!

  void remove(const std::string& name);

private:

  typedef std::map< std::string, eckit::SharedPtr<Field> >  FieldMap;

private:

  FieldMap fields_;
  Metadata metadata_;

};

//------------------------------------------------------------------------------------------------------

class StateGenerator : public eckit::Owned {

public:

    StateGenerator( const eckit::Parametrisation& = Config() );

    virtual ~StateGenerator();

    virtual void generate( State&, const eckit::Parametrisation& = Config() ) const = 0;

};

//------------------------------------------------------------------------------------------------------

class StateGeneratorFactory {
  public:

    /*!
     * \brief build StateCreator with options specified in parametrisation
     * \return mesh generator
     */
    static StateGenerator* build(const std::string& state_generator,
                                 const eckit::Parametrisation& = Config() );

    /*!
     * \brief list all registered field creators
     */
    static void list(std::ostream &);
    static bool has(const std::string& name);

  private:

    virtual StateGenerator* make(const eckit::Parametrisation& = Config() ) = 0 ;

    std::string name_;

  protected:

    StateGeneratorFactory(const std::string&);
    virtual ~StateGeneratorFactory();
};


template<class T>
class StateGeneratorBuilder : public StateGeneratorFactory {

  virtual StateGenerator* make(const eckit::Parametrisation& param = Config() ) {
        return new T(param);
  }
  public:
    StateGeneratorBuilder(const std::string& name) : StateGeneratorFactory(name) {}
};

// ------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define eckit_Parametrisation eckit::Parametrisation
extern "C"
{
  State* atlas__State__new ();
  void atlas__State__initialize (State* This, const char* generator, const eckit_Parametrisation* params);
  void atlas__State__delete (State* This);
  void atlas__State__add (State* This, Field* field);
  void atlas__State__remove (State* This, const char* name);
  int atlas__State__has (State* This, const char* name);
  Field* atlas__State__field_by_name (State* This, const char* name);
  Field* atlas__State__field_by_index (State* This, int index);
  int atlas__State__size(const State* This);
  Metadata* atlas__State__metadata (State* This);
}
#undef eckit_Parametrisation

} // namespace atlas


#endif
