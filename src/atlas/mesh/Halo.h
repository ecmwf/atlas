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

#include <cstddef>

namespace atlas {
class Mesh;
namespace mesh {
namespace detail {
class MeshImpl;
}
}  // namespace mesh
}  // namespace atlas

namespace atlas {
namespace mesh {

// -------------------------------------------------------------------

class Halo {
public:
    Halo() {}
    Halo(const Mesh& mesh);
    Halo(const detail::MeshImpl& mesh);
    Halo(const int size): size_(size) {}
    int size() const;

private:
    int size_{-1};
};

}  // namespace mesh
}  // namespace atlas
