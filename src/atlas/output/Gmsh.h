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

namespace eckit {
class Parametrisation;
class PathName;
}  // namespace eckit

namespace atlas {
namespace output {
namespace detail {
class GmshIO;
}
}  // namespace output
}  // namespace atlas

namespace atlas {
namespace output {

// -----------------------------------------------------------------------------

class GmshFileStream : public std::ofstream {
public:
    static std::string parallelPathName(const eckit::PathName& path, int part = mpi::rank());
    GmshFileStream(const eckit::PathName& file_path, const char* mode, int part = mpi::rank());
};

// -----------------------------------------------------------------------------

class Gmsh : public Output {
public:
    Gmsh(const Output& output);
    Gmsh(std::ostream&);
    Gmsh(std::ostream&, const eckit::Parametrisation&);

    Gmsh(const eckit::PathName&, const std::string& mode);
    Gmsh(const eckit::PathName&, const std::string& mode, const eckit::Parametrisation&);

    Gmsh(const eckit::PathName&);
    Gmsh(const eckit::PathName&, const eckit::Parametrisation&);
};

// -----------------------------------------------------------------------------

}  // namespace output
}  // namespace atlas
