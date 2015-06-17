/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Grid.h"

#include <vector>

#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/config/Resource.h"

#include "atlas/Mesh.h"
#include "atlas/GridSpec.h"
#include "atlas/grids/grids.h"

using eckit::Params;
using eckit::Factory;

namespace atlas {

//------------------------------------------------------------------------------------------------------

static void checkSizeOfPoint()
{
    // compilen time check support C++11
    #if __cplusplus >= 201103L
        static_assert( sizeof(Grid::Point)==2*double, "Grid requires size of Point to be 2*double" );
    #endif

    // runtime check
    ASSERT( sizeof(Grid::Point) == 2*sizeof(double) );
}


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

Grid::Grid() : domain_( Domain::makeGlobal() )
{
    checkSizeOfPoint();
}

Grid::Grid(const Domain& domain) : domain_(domain)
{
    checkSizeOfPoint();
}

Grid::~Grid() {}

const atlas::Domain& Grid::domain() const
{
  return domain_;
}

Grid::uid_t Grid::uniqueId() const {
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

BoundBox Grid::boundingBox() const
{
    NOTIMP;
    return BoundBox();
}

void Grid::fillLonLat(double array[], size_t arraySize) const
{
    const size_t size = npts()*2;
    ASSERT(arraySize >= size);
    copyLonLatMemory(array, size_t(sizeof(double)*size));
}

void Grid::fillLonLat(std::vector<double>& v) const {
    v.resize(npts()*2);
    copyLonLatMemory(&v[0], size_t(sizeof(double)*v.size()));
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

size_t Grid::copyLonLatMemory(double* pts, size_t size) const
{
    std::vector<Grid::Point> gpts;
    lonlat(gpts);

    size_t sizePts = 2*npts();

    ASSERT(size >= sizePts);

    for(size_t c = 0, i = 0; i < gpts.size(); ++i) {
        pts[c++] = gpts[i].lon();
        pts[c++] = gpts[i].lat();
    }

    return sizePts;
}

bool Grid::same(const Grid& g) const { return uniqueId() == g.uniqueId(); }

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
