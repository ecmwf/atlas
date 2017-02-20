/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/CustomStructured.h"
#include "atlas/internals/Debug.h"


namespace atlas {
namespace grid {


register_BuilderT1(Grid, CustomStructured, CustomStructured::grid_type_str());


std::string CustomStructured::className() {
    return "atlas.grid.CustomStructured";
}


std::string CustomStructured::grid_type_str() {
    return "structured";
}


CustomStructured::CustomStructured(const util::Config& config) :
    Structured()
{

  util::Config config_proj;
  if( dynamic_cast<const util::Config&>(config).get("projection",config_proj) )
    projection_.reset( projection::Projection::create(config_proj) );
  else
    projection_.reset( projection::Projection::create() );


  util::Config config_yspace;
  if( not config.get("yspace",config_yspace) )
    throw eckit::BadParameter("yspace missing in configuration");

  eckit::SharedPtr<spacing::Spacing> yspace( spacing::Spacing::create(config_yspace) );

  const size_t ny = yspace->size();
  std::vector<long>   nx;       nx.reserve(ny);
  std::vector<double> xmin;    xmin.reserve(ny);
  std::vector<double> xmax;    xmax.reserve(ny);


  double dom_xmin =  std::numeric_limits<double>::max();
  double dom_xmax = -std::numeric_limits<double>::max();

  std::vector<util::Config> config_xspace_list;
  if( config.get("xspace[]",config_xspace_list) ) {
    
    ASSERT( config_xspace_list.size() == ny );
    std::string xspace_type;

    for( size_t j=0; j<ny; ++j ) {
      config_xspace_list[j].get("type",xspace_type);
      ASSERT( xspace_type == "linear" );
      eckit::SharedPtr<spacing::Spacing> xspace( spacing::Spacing::create(config_xspace_list[j]) );
      nx.push_back(xspace->size());
      xmin.push_back(xspace->front());
      xmax.push_back(xspace->back());
      dom_xmin = std::min(dom_xmin,xspace->min());
      dom_xmax = std::max(dom_xmax,xspace->max());
    }

  } else {

    util::Config config_xspace;
    if( not config.get("xspace",config_xspace) )
      throw eckit::BadParameter("xspace missing in configuration");

    std::string xspace_type;
    config_xspace.get("type",xspace_type);
    ASSERT( xspace_type == "linear" );

    std::vector<long>   v_N;
    std::vector<double> v_start;
    std::vector<double> v_end;
    std::vector<double> v_length;
    config_xspace.get("N[]",      v_N     );
    config_xspace.get("start[]",  v_start );
    config_xspace.get("end[]",    v_end   );
    config_xspace.get("length[]", v_length);

    if( not v_N.     empty() ) ASSERT(v_N.     size() == ny);
    if( not v_start. empty() ) ASSERT(v_start. size() == ny);
    if( not v_end.   empty() ) ASSERT(v_end.   size() == ny);
    if( not v_length.empty() ) ASSERT(v_length.size() == ny);
    
    for( size_t j=0; j<ny; ++j ) {
      if( not v_N.     empty() ) config_xspace.set("N",     v_N[j]);
      if( not v_start. empty() ) config_xspace.set("start", v_start[j]);
      if( not v_end.   empty() ) config_xspace.set("end",   v_end[j]);
      if( not v_length.empty() ) config_xspace.set("length",v_length[j]);
      eckit::SharedPtr<spacing::Spacing> xspace( spacing::Spacing::create(config_xspace) );
      nx.push_back(xspace->size());
      xmin.push_back(xspace->front());
      xmax.push_back(xspace->back());
      dom_xmin = std::min(dom_xmin,xspace->min());
      dom_xmax = std::max(dom_xmax,xspace->max());
    }
  }
  
  std::vector<double> y; y.assign(yspace->begin(),yspace->end());
  
  
  util::Config config_domain;
  if( dynamic_cast<const util::Config&>(config).get("domain",config_domain) )
    domain_.reset( domain::Domain::create(config_domain) );
  else {
    config_domain.set("type","rectangular");
    config_domain.set("ymin",yspace->min());
    config_domain.set("ymax",yspace->max());
    config_domain.set("xmin",dom_xmin);
    config_domain.set("xmax",dom_xmax);
    domain_.reset( domain::Domain::create(config_domain) );
  }
  
  Structured::setup(ny, y.data(), nx.data(), xmin.data(), xmax.data());
}



CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long pl[]) :
    Structured()
{
  NOTIMP;
    ASSERT(nlat);

    util::Config config_domain;
    config_domain.set("type","global");
    domain_.reset( domain::Domain::create(config_domain) );

    util::Config config_proj;
    config_proj.set("type","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // assign longitude limits
    std::vector<double> lonmin(nlat,std::numeric_limits<double>::max());
    std::vector<double> lonmax(nlat,-std::numeric_limits<double>::max());
    ASSERT(lonmin.size() == nlat);
    ASSERT(lonmax.size() == nlat);

    setup_cropped(nlat, lats, pl, lonmin.data(), lonmax.data(), *domain_ );

}


CustomStructured::CustomStructured(
    size_t nlat,
    const double latitudes[],
    const long pl[],
    const double lonmin[],
    const double lonmax[] ) :
    Structured()
{
  NOTIMP;
  
    Structured::setup(nlat, latitudes, pl, lonmin, lonmax);
}

extern "C" {


    Structured* atlas__grid__CustomStructured_int(size_t nlat, double lats[], int pl[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl,pl+nlat);
        return new CustomStructured(nlat, lats, pl_vector.data());
    }


    Structured* atlas__grid__CustomStructured_long(size_t nlat, double lats[], long pl[]) {
        return new CustomStructured(nlat, lats, pl);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(size_t nlat, double lats[], int pl[], double lonmin[], double lonmax[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl, pl+nlat);
        return new CustomStructured(nlat, lats, pl_vector.data(), lonmin, lonmax);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(size_t nlat, double lats[], long pl[], double lonmin[], double lonmax[]) {
        return new CustomStructured(nlat, lats, pl, lonmin, lonmax);
    }


}


}  // namespace grid
}  // namespace atlas

