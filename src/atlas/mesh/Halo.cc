/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/Halo.h"
#include "atlas/Mesh.h"
#include "atlas/Metadata.h"

namespace atlas {
namespace mesh {

Halo::Halo(const Mesh& mesh)
{
  size_=0;
  mesh.metadata().get("halo",size_);
}

} // namespace mesh
} // namespace atlas

