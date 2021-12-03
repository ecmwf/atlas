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

namespace atlas {
class Mesh;
namespace mesh {
namespace actions {

void write_load_balance_report(const Mesh& mesh, std::ostream& ofs);
void write_load_balance_report(const Mesh& mesh, const std::string& filename);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
void atlas__write_load_balance_report(Mesh::Implementation* mesh, char* filename);
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
