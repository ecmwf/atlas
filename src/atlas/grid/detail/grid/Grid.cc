/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/detail/grid/Grid.h"

#include <vector>
#include "eckit/memory/Factory.h"
#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {


static void checkSizeOfPoint() {
    // compile time check support C++11
#if __cplusplus >= 201103L
    static_assert( sizeof(Grid::Point)==2*sizeof(double), "Grid requires size of Point to be 2*double" );
#endif

    // runtime check
    ASSERT( sizeof(Grid::Point) == 2*sizeof(double) );
}


std::string Grid::className() {
    return "atlas.Grid";
}


const Grid* Grid::create(const Config& p) {

    eckit::Factory<Grid>& fact = eckit::Factory<Grid>::instance();

    std::string name;
    if (p.get("name",name)) {
        return create(name);
    }
    std::string type;
    if (p.get("type",type) && fact.exists(type)) {
        return fact.get(type).create(p);
    }

    if( name.size() ) {
      Log::info() << "name provided: " << name << std::endl;
    }
    if( type.size() ) {
      Log::info() << "type provided: " << type << std::endl;
    }
    if( name.empty() && type.empty() ) {
      throw eckit::BadParameter("no name or type in configuration",Here());
    } else {
      throw eckit::BadParameter("name or type in configuration don't exist",Here());
    }

    return NULL;
}


const Grid* Grid::create( const std::string& name ) {

      const GridCreator::Registry& registry = GridCreator::registry();
      for( GridCreator::Registry::const_iterator it = registry.begin(); it!=registry.end(); ++it ) {
        const Grid* grid = it->second->create(name);
        if( grid ) {
          return grid;
        }
      }

      // Throw exception
      std::ostringstream log;
      log << "Could not construct Grid from the name \""<< name<< "\"\n";
      log << "Accepted names are: \n";
      for( GridCreator::Registry::const_iterator it = registry.begin(); it!=registry.end(); ++it ) {
        log << "  -  " << *it->second << "\n";
      }
      throw eckit::BadParameter(log.str());
//    return GridCreator::createNamed(name);
}


Grid::Grid()
{
    checkSizeOfPoint();
}

Grid::~Grid() {
}


Grid::uid_t Grid::uid() const {
    if (uid_.empty()) {
        std::ostringstream s;
        s << hash();
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


std::string Grid::getOptimalMeshGenerator() const {
    return "Delaunay";
}


bool Grid::same(const grid::Grid& g) const {
    return uid() == g.uid();
}

Grid::Spec Grid::spec() const {
    eckit::Properties grid_spec, dom_spec, proj_spec;

    // grid type
    grid_spec.set("name", name());

    // add domain specs
    dom_spec=domain().spec();
    grid_spec.set("domain",dom_spec);

    // add projection specs
    proj_spec=projection().spec();
    grid_spec.set("projection",proj_spec);

    return grid_spec;
}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
