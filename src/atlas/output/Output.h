/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_output_Output_h
#define atlas_output_Output_h

#include <iosfwd>
#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/serialisation/FileStream.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/util/Config.h"

namespace eckit {
  class Parametrisation;
  class PathName;
}

namespace atlas {
namespace mesh {
    class Mesh;
} }

namespace atlas {
namespace field {
  class Field;
  class FieldSet;
} }

namespace atlas {
namespace functionspace {
    class FunctionSpace;
} }

namespace atlas {
namespace output {
  
typedef std::ostream Stream;
typedef eckit::PathName PathName;

// -----------------------------------------------------------------------------

class Output : public eckit::Owned {

public:

  typedef eckit::SharedPtr<Output> Ptr;
  typedef atlas::util::Config Parameters;
  static Output* create(const std::string&, Stream&, const eckit::Parametrisation & = util::NoConfig() );

public:

    Output();

    virtual ~Output();

    /// Write mesh file
    virtual void write(
        const mesh::Mesh&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write field to file
    virtual void write(
        const field::Field&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(
        const field::FieldSet&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write field to file using Functionspace
    virtual void write(
        const field::Field&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(
        const field::FieldSet&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

};



class OutputFactory {
  public:
    /*!
     * \brief build Output with factory key, and default options
     * \return mesh generator
     */
    static Output* build(const std::string&, Stream&);

    /*!
     * \brief build Output with factory key inside parametrisation,
     * and options specified in parametrisation as well
     * \return mesh generator
     */
    static Output* build(const std::string&, Stream&, const eckit::Parametrisation&);

    /*!
     * \brief list all registered mesh generators
     */
    static void list(std::ostream &);

  private:
    std::string name_;
    virtual Output* make(Stream&) = 0 ;
    virtual Output* make(Stream&, const eckit::Parametrisation&) = 0 ;

  protected:

    OutputFactory(const std::string&);
    virtual ~OutputFactory();

};


template<class T>
class OutputBuilder : public OutputFactory {
  virtual Output* make(Stream& stream) {
      return new T(stream);
  }
  virtual Output* make(Stream& stream, const eckit::Parametrisation& param) {
        return new T(stream,param);
  }
  public:
    OutputBuilder(const std::string& name) : OutputFactory(name) {}
};

// -----------------------------------------------------------------------------

extern "C" {
void atlas__Output__delete(Output* This);
Output* atlas__Output__create(const char* factory_key, Stream* stream, const eckit::Parametrisation* params);
void atlas__Output__write_mesh(const Output* This, mesh::Mesh::mesh_t* mesh, const eckit::Parametrisation* params);
void atlas__Output__write_fieldset(const Output* This, const field::FieldSet* fieldset, const eckit::Parametrisation* params);
void atlas__Output__write_field(const Output* This, const field::Field* field, const eckit::Parametrisation* params);
void atlas__Output__write_fieldset_fs(const Output* This, const field::FieldSet* fieldset, const functionspace::FunctionSpace* functionspace, const eckit::Parametrisation* params);
void atlas__Output__write_field_fs(const Output* This, const field::Field* field, const functionspace::FunctionSpace* functionspace, const eckit::Parametrisation* params);
}

} // namespace output
} // namespace atlas

#endif
