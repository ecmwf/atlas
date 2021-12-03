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

#include <fstream>

#include "atlas/output/Output.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace output {
namespace detail {
class GmshIO;
}
}  // namespace output
}  // namespace atlas

namespace atlas {
namespace output {
namespace detail {

// -----------------------------------------------------------------------------

class GmshImpl : public OutputImpl {
public:
    GmshImpl(std::ostream&);
    GmshImpl(std::ostream&, const eckit::Parametrisation&);

    GmshImpl(const eckit::PathName&, const std::string& mode);
    GmshImpl(const eckit::PathName&, const std::string& mode, const eckit::Parametrisation&);

    GmshImpl(const eckit::PathName&);
    GmshImpl(const eckit::PathName&, const eckit::Parametrisation&);

    virtual ~GmshImpl();

    /// Write mesh file
    virtual void write(const Mesh&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write field to file
    virtual void write(const Field&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write fieldset to file using FunctionSpace
    virtual void write(const FieldSet&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write field to file using Functionspace
    virtual void write(const Field&, const FunctionSpace&, const eckit::Parametrisation& = util::NoConfig()) const;

    /// Write fieldset to file using FunctionSpace
    virtual void write(const FieldSet&, const FunctionSpace&, const eckit::Parametrisation& = util::NoConfig()) const;

public:
    struct Configuration {
        bool binary;
        bool edges;
        bool elements;
        bool gather;
        bool ghost;
        bool info;
        std::vector<long> levels;
        std::string nodes;
        std::string file;
        std::string openmode;
        std::string coordinates;
    };

    static void setGmshConfiguration(detail::GmshIO&, const Configuration&);

private:
    mutable Configuration config_;

    void defaults();
};

// -----------------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
