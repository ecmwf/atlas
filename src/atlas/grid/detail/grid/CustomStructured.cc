/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/detail/grid/CustomStructured.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/internals/Debug.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {


register_BuilderT1(Grid, CustomStructured, CustomStructured::static_type());

CustomStructured::CustomStructured(const util::Config& config) :
    Structured()
{

  util::Config config_proj;
  if( dynamic_cast<const util::Config&>(config).get("projection",config_proj) )
    projection_ = Projection(config_proj);
  else
    projection_ = Projection();


  util::Config config_yspace;
  if( not config.get("yspace",config_yspace) )
    throw eckit::BadParameter("yspace missing in configuration");

  YSpace yspace(config_yspace);

  const size_t ny = yspace.size();
  std::vector<long>   nx;      nx.  reserve(ny);
  std::vector<double> xmin;    xmin.reserve(ny);
  std::vector<double> xmax;    xmax.reserve(ny);
  std::vector<double> dx;      dx  .reserve(ny);


  double dom_xmin =  std::numeric_limits<double>::max();
  double dom_xmax = -std::numeric_limits<double>::max();

  std::vector<util::Config> config_xspace_list;
  if( config.get("xspace[]",config_xspace_list) ) {

    ASSERT( config_xspace_list.size() == ny );
    std::string xspace_type;

    for( size_t j=0; j<ny; ++j ) {
      config_xspace_list[j].get("type",xspace_type);
      ASSERT( xspace_type == "linear" );
      spacing::LinearSpacing::Params xspace( config_xspace_list[j] );
      xmin.push_back(xspace.start);
      xmax.push_back(xspace.end);
      nx.push_back(xspace.N);
      dx.push_back(xspace.step);
      dom_xmin = std::min(dom_xmin,xspace.start);
      dom_xmax = std::max(dom_xmax,xspace.end);
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
      spacing::LinearSpacing::Params xspace( config_xspace );
      xmin.push_back(xspace.start);
      xmax.push_back(xspace.end);
      nx.push_back(xspace.N);
      dx.push_back(xspace.step);
      dom_xmin = std::min(dom_xmin,xspace.start);
      dom_xmax = std::max(dom_xmax,xspace.end);
    }
  }

  util::Config config_domain;
  if( dynamic_cast<const util::Config&>(config).get("domain",config_domain) )
    domain_ = Domain(config_domain);
  else {
    config_domain.set("type","rectangular");
    config_domain.set("ymin",yspace.min());
    config_domain.set("ymax",yspace.max());
    config_domain.set("xmin",dom_xmin);
    config_domain.set("xmax",dom_xmax);
    config_domain.set("units",projection_.units());
    domain_ = Domain(config_domain);
  }

  Structured::setup(yspace, nx, xmin, xmax, dx);
}



CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long pl[]) :
    Structured()
{
    ASSERT(nlat);
    YSpace yspace( new spacing::CustomSpacing(nlat, lats) );

    util::Config config_domain;
    config_domain.set("type","global");
    domain_ = Domain(config_domain);

    util::Config config_proj;
    config_proj.set("type","lonlat");
    projection_ = Projection(config_proj);

    std::vector<long> nx;
    nx.assign(pl,pl+nlat);

    std::vector<double> xmin(nlat, 0.);
    std::vector<double> xmax(nlat, 360.);
    std::vector<double> dx(nlat);
    for( size_t j=0; j<nlat; ++j) {
      if( not nx[j] ) dx[j] = 0;
      else dx[j] = 360./double(nx[j]);
    }
    setup(yspace,nx,xmin,xmax,dx);
}

extern "C" {


    Structured* atlas__grid__CustomStructured_int(long nlat, double lats[], int pl[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl,pl+nlat);
        return new CustomStructured(nlat, lats, pl_vector.data());
    }


    Structured* atlas__grid__CustomStructured_long(long nlat, double lats[], long pl[]) {
        return new CustomStructured(nlat, lats, pl);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(long nlat, double lats[], int pl[], double lonmin[], double lonmax[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl, pl+nlat);
        NOTIMP;
        return 0;
        // return new CustomStructured(nlat, lats, pl_vector.data(), lonmin, lonmax);
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(long nlat, double lats[], long pl[], double lonmin[], double lonmax[]) {
        NOTIMP;
        return 0;
        //return new CustomStructured(nlat, lats, pl, lonmin, lonmax);
    }


}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

