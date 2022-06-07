/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>
#include <string>

#include "eckit/config/Parametrisation.h"
#include "eckit/serialisation/FileStream.h"

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Parametrisation;
class PathName;
}  // namespace eckit

namespace atlas {
class Mesh;
namespace mesh {
namespace detail {
class MeshImpl;
}
}  // namespace mesh
}  // namespace atlas

namespace atlas {
class Field;
class FieldSet;
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field
}  // namespace atlas

namespace atlas {
class FunctionSpace;
namespace functionspace {
class FunctionSpaceImpl;
}
}  // namespace atlas

namespace atlas {
namespace output {

using PathName = eckit::PathName;

// -----------------------------------------------------------------------------

namespace detail {
class OutputImpl : public util::Object {
public:
    typedef atlas::util::Config Parameters;

public:
    OutputImpl();

    virtual ~OutputImpl();

    /// Write mesh file
    virtual void write(const Mesh&, const eckit::Parametrisation& = util::NoConfig()) const = 0;

    /// Write field to file
    virtual void write(const Field&, const eckit::Parametrisation& = util::NoConfig()) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(const FieldSet&, const eckit::Parametrisation& = util::NoConfig()) const = 0;

    /// Write field to file using Functionspace
    virtual void write(const Field&, const FunctionSpace&, const eckit::Parametrisation& = util::NoConfig()) const = 0;

    /// Write fieldset to file using FunctionSpace
    virtual void write(const FieldSet&, const FunctionSpace&,
                       const eckit::Parametrisation& = util::NoConfig()) const = 0;
};

}  // namespace detail

class Output : DOXYGEN_HIDE(public util::ObjectHandle<detail::OutputImpl>) {
public:
    using Handle::Handle;
    Output() = default;
    Output(const std::string&, std::ostream&, const eckit::Parametrisation& = util::NoConfig());

    /// Write mesh file
    const Output& write(const Mesh&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write field to file
    const Output& write(const Field&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write fieldset to file using FunctionSpace
    const Output& write(const FieldSet&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write field to file using Functionspace
    const Output& write(const Field&, const FunctionSpace&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write fieldset to file using FunctionSpace
    const Output& write(const FieldSet&, const FunctionSpace&, const eckit::Parametrisation& = util::NoConfig()) const;
};

namespace detail {

class OutputFactory {
public:
    /*!
   * \brief build Output with factory key, and default options
   * \return mesh generator
   */
    static const OutputImpl* build(const std::string&, std::ostream&);

    /*!
   * \brief build Output with factory key inside parametrisation,
   * and options specified in parametrisation as well
   * \return mesh generator
   */
    static const OutputImpl* build(const std::string&, std::ostream&, const eckit::Parametrisation&);

    /*!
   * \brief list all registered mesh generators
   */
    static void list(std::ostream&);

private:
    std::string name_;
    virtual const OutputImpl* make(std::ostream&)                                = 0;
    virtual const OutputImpl* make(std::ostream&, const eckit::Parametrisation&) = 0;

protected:
    OutputFactory(const std::string&);
    virtual ~OutputFactory();
};

template <class T>
class OutputBuilder : public OutputFactory {
    virtual const OutputImpl* make(std::ostream& stream) { return new T(stream); }
    virtual const OutputImpl* make(std::ostream& stream, const eckit::Parametrisation& param) {
        return new T(stream, param);
    }

public:
    OutputBuilder(const std::string& name): OutputFactory(name) {}
};

// -----------------------------------------------------------------------------

extern "C" {
void atlas__Output__delete(OutputImpl* This);
const OutputImpl* atlas__Output__create(const char* factory_key, std::ostream* stream,
                                        const eckit::Parametrisation* params);
void atlas__Output__write_mesh(const OutputImpl* This, mesh::detail::MeshImpl* mesh,
                               const eckit::Parametrisation* params);
void atlas__Output__write_fieldset(const OutputImpl* This, const field::FieldSetImpl* fieldset,
                                   const eckit::Parametrisation* params);
void atlas__Output__write_field(const OutputImpl* This, const field::FieldImpl* field,
                                const eckit::Parametrisation* params);
void atlas__Output__write_fieldset_fs(const OutputImpl* This, const field::FieldSetImpl* fieldset,
                                      const functionspace::FunctionSpaceImpl* functionspace,
                                      const eckit::Parametrisation* params);
void atlas__Output__write_field_fs(const OutputImpl* This, const field::FieldImpl* field,
                                   const functionspace::FunctionSpaceImpl* functionspace,
                                   const eckit::Parametrisation* params);
}

}  // namespace detail
}  // namespace output
}  // namespace atlas
