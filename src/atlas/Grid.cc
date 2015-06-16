/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "Grid.h"

#include <vector>

#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/config/Resource.h"

#include "atlas/Mesh.h"
#include "atlas/GridSpec.h"
//#include "atlas/Tesselation.h"
#include "atlas/grids/grids.h"

using eckit::Params;
using eckit::Factory;

namespace atlas {

Grid& Grid::from_id(Id id)
{
  return *registry().get(id);
}

Grid::Registry& Grid::registry()
{
  static Registry r;
  return r;
}

//------------------------------------------------------------------------------------------------------

Grid* Grid::create(const Params& p) {
  if (p.has("shortName")) {
    if (Factory<Grid>::instance().exists(p["shortName"])) {
      return Factory<Grid>::instance().get(p["shortName"]).create(p);
    }
  }
  return Factory<Grid>::instance().get(p["grid_type"]).create(p);
}

Grid* Grid::create(const Grid::uid_t& uid) { return grids::grid_from_uid(uid); }

Grid* Grid::create(const GridSpec& g) { return Grid::create(Params(g)); }

Grid::Grid() { registry_id_ = registry().add(*this); }

Grid::~Grid() { registry().remove(registry_id_); }

Grid::Domain Grid::domain() const { return bounding_box(); }

void Grid::mask(const Domain&) { NOTIMP; }

void Grid::mask(const Params&) { NOTIMP; }

Grid::uid_t Grid::unique_id() const {
  if (uid_.empty()) {
    std::ostringstream s;
    s << shortName() << "." << hash();
    uid_ = s.str();
  }
  return uid_;
}

eckit::MD5::digest_t Grid::hash() const {
  if (hash_.empty()) {
    eckit::MD5 md5;
    hash(md5);
    hash_ = md5.digest();
  }
  return hash_;
}

Grid* Grid::masked(const Domain&) const {
  NOTIMP;
  return NULL;
}

Grid* Grid::masked(const eckit::Params&) const {
  NOTIMP;
  return NULL;
}

void Grid::lonlat(std::vector<double>& crd) const {
  crd.resize(npts() * 2);
  lonlat(crd.data());
}

void Grid::set_mesh(const Mesh& mesh)
{
  mesh_ = eckit::SharedPtr<Mesh>( const_cast<Mesh*>(&mesh) );
}

Mesh& Grid::mesh() const {
  if( !mesh_ )
  {
    mesh_.reset( new Mesh() );
    mesh_->add_nodes(*this);
  }
  return *mesh_;
}

bool Grid::same(const Grid& g) const { return unique_id() == g.unique_id(); }

//------------------------------------------------------------------------------------------------------

BoundBox Grid::make_global_bounding_box() { return BoundBox(90., -90., 360. - degrees_eps(), 0.); }

BoundBox Grid::make_bounding_box(const Params& p) {
  if (!p.has("bbox_s")) return Grid::make_global_bounding_box();

  return BoundBox(p["bbox_n"], p["bbox_s"], p["bbox_e"], p["bbox_w"]);
}

void Grid::print(std::ostream &) const
{
  NOTIMP;
}

double Grid::degrees_eps() {
  /// default is 1E-3 because
  /// some bugs in IFS means we need a lower resolution epsilon when decoding from grib2

  static double eps = eckit::Resource<double>("$ATLAS_DEGREES_EPSILON;atlas.degrees_epsilon", 1E-3);
  return eps;
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
