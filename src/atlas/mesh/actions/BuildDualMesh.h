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

#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {

class Grid;

namespace mesh {
namespace actions {

/*
 * Build median-dual mesh by connecting cell centres with edge centres
 * - dual_volumes field
 * - dual_normals field
 */
void build_median_dual_mesh(Mesh& mesh);

/*
 * Build centroid-dual mesh by connecting cell centres with cell centres
 * - dual_volumes field
 * - dual_normals field
 * - skewness field value 1 at ip1, value -1 at ip2
 * - alpha field  value 1 at ip1, value 0 at ip2
 */
void build_centroid_dual_mesh(Mesh& mesh);

/*
 * Build brick-dual mesh
 * - dual_volumes field
 * @note Only works for reduced grids with 1 MPI task
 *       This was urgently coded to show the resolution as IFS sees it.
 *
 *  #   |   #   |   #   |   #   |
 *      |       |       |       |
 *  ----+-+-----+-+-----+-+-----+-+
 *        |       |       |       |
 *    #   |   #   |   #   |   #   |
 *        |       |       |       |
 *  -+----+--+----+--+----+--+----+
 *   |       |       |       |
 *   |   #   |   #   |   #   |   #
 */
void build_brick_dual_mesh(const Grid& grid, Mesh& mesh);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {
void atlas__build_median_dual_mesh(Mesh::Implementation* mesh);
void atlas__build_centroid_dual_mesh(Mesh::Implementation* mesh);
}
// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
