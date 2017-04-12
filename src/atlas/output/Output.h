/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>
#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
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
  class Field;
  class FieldSet;
namespace field {
  class FieldImpl;
  class FieldSetImpl;
} }

namespace atlas {
namespace functionspace {
    class FunctionSpace;
    class FunctionSpaceImpl;
} }

namespace atlas {
namespace output {

typedef std::ostream Stream;
typedef eckit::PathName PathName;

// -----------------------------------------------------------------------------

class OutputImpl : public eckit::Owned {

public:

  typedef atlas::util::Config Parameters;

public:

    OutputImpl();

    virtual ~OutputImpl();

    /// Write mesh file
    virtual void write(
        const mesh::Mesh&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write field to file
    virtual void write(
        const Field&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(
        const FieldSet&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write field to file using Functionspace
    virtual void write(
        const Field&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(
        const FieldSet&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const = 0;

};

class Output {

public:

     using output_t = OutputImpl;

private:

     eckit::SharedPtr<const output_t> output_;

public:

    Output();
    Output(const output_t*);
    Output(const Output&);
    Output(const std::string&, Stream&, const eckit::Parametrisation & = util::NoConfig() );

    /// Write mesh file
    void write(
        const mesh::Mesh&,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// Write field to file
    void write(
        const Field&,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// Write fieldset to file using FunctionSpace
    void write(
        const FieldSet&,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// Write field to file using Functionspace
    void write(
        const Field&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// Write fieldset to file using FunctionSpace
    void write(
        const FieldSet&,
        const functionspace::FunctionSpace&,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    const output_t* get() const { return output_.get(); }
};





class OutputFactory {
  public:
    /*!
     * \brief build Output with factory key, and default options
     * \return mesh generator
     */
    static const OutputImpl* build(const std::string&, Stream&);

    /*!
     * \brief build Output with factory key inside parametrisation,
     * and options specified in parametrisation as well
     * \return mesh generator
     */
    static const OutputImpl* build(const std::string&, Stream&, const eckit::Parametrisation&);

    /*!
     * \brief list all registered mesh generators
     */
    static void list(std::ostream &);

  private:
    std::string name_;
    virtual const OutputImpl* make(Stream&) = 0 ;
    virtual const OutputImpl* make(Stream&, const eckit::Parametrisation&) = 0 ;

  protected:

    OutputFactory(const std::string&);
    virtual ~OutputFactory();

};


template<class T>
class OutputBuilder : public OutputFactory {
  virtual const OutputImpl* make(Stream& stream) {
      return new T(stream);
  }
  virtual const OutputImpl* make(Stream& stream, const eckit::Parametrisation& param) {
        return new T(stream,param);
  }
  public:
    OutputBuilder(const std::string& name) : OutputFactory(name) {}
};

// -----------------------------------------------------------------------------

extern "C" {
void atlas__Output__delete(OutputImpl* This);
const OutputImpl* atlas__Output__create(const char* factory_key, Stream* stream, const eckit::Parametrisation* params);
void atlas__Output__write_mesh(const OutputImpl* This, mesh::Mesh::Implementation* mesh, const eckit::Parametrisation* params);
void atlas__Output__write_fieldset(const OutputImpl* This, const field::FieldSetImpl* fieldset, const eckit::Parametrisation* params);
void atlas__Output__write_field(const OutputImpl* This, const field::FieldImpl* field, const eckit::Parametrisation* params);
void atlas__Output__write_fieldset_fs(const OutputImpl* This, const field::FieldSetImpl* fieldset, const functionspace::FunctionSpaceImpl* functionspace, const eckit::Parametrisation* params);
void atlas__Output__write_field_fs(const OutputImpl* This, const field::FieldImpl* field, const functionspace::FunctionSpaceImpl* functionspace, const eckit::Parametrisation* params);
}

} // namespace output
} // namespace atlas
