/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_actions_AddVirtualNodes_h
#define atlas_actions_AddVirtualNodes_h


namespace atlas {
namespace mesh { class Mesh; }
namespace grid { class Grid; }
}


namespace atlas {
namespace mesh {
namespace actions {


/// Adds virtual nodes to the mesh that aren't contained in the Grid Domain
class AddVirtualNodes {
public:

    void operator()(const atlas::grid::Grid&, atlas::mesh::Mesh&) const;

};


} // namespace actions
} // namespace mesh
} // namespace atlas


#endif
