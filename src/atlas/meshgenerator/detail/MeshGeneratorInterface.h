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

namespace eckit {
class Parametrisation;
}  // namespace eckit

namespace atlas {
namespace grid {
class DistributionImpl;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
class Grid;
}  // namespace grid
}  // namespace detail
}  // namespace grid
using GridImpl = grid::detail::grid::Grid;
}  // namespace atlas
namespace atlas {
namespace mesh {
namespace detail {
class MeshImpl;
}  // namespace detail
}  // namespace mesh
}  // namespace atlas
namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {
class Partitioner;
}  // namespace partitioner
}  // namespace detail
}  // namespace grid
using PartitionerImpl = grid::detail::partitioner::Partitioner;
}  // namespace atlas


namespace atlas {
namespace meshgenerator {

class MeshGeneratorImpl;

//----------------------------------------------------------------------------------------------------------------------

extern "C" {
void atlas__MeshGenerator__delete(MeshGeneratorImpl* This);
const MeshGeneratorImpl* atlas__MeshGenerator__create_noconfig(const char* name);
const MeshGeneratorImpl* atlas__MeshGenerator__create(const char* name, const eckit::Parametrisation* params);
mesh::detail::MeshImpl* atlas__MeshGenerator__generate__grid_griddist(const MeshGeneratorImpl* This,
                                                                      const GridImpl* grid,
                                                                      const grid::DistributionImpl* distribution);
mesh::detail::MeshImpl* atlas__MeshGenerator__generate__grid(const MeshGeneratorImpl* This, const GridImpl* grid);
mesh::detail::MeshImpl* atlas__MeshGenerator__generate__grid_partitioner(const MeshGeneratorImpl* This,
                                                                         const GridImpl* grid,
                                                                         const PartitionerImpl* partitioner);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
