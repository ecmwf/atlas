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


namespace atlas {
namespace grid {


register_BuilderT1(Grid, CustomStructured, CustomStructured::grid_type_str());


std::string CustomStructured::className() {
    return "atlas.grid.CustomStructured";
}


std::string CustomStructured::grid_type_str() {
    return "custom_structured";
}


CustomStructured::CustomStructured(const eckit::Parametrisation& params) :
    Structured()
{
  //domain_ = domain::Domain::makeGlobal();

  if( not params.get("periodic_x",periodic_x_) ) periodic_x_ = true;
  if( not params.get("periodic_y",periodic_y_) ) periodic_y_ = false;

  util::Config config_domain;
  if( dynamic_cast<const util::Config&>(params).get("domain",config_domain) )
    domain_.reset( domain::Domain::create(config_domain) );
  else
    domain_.reset( domain::Domain::create() );

  util::Config config_proj;
  if( dynamic_cast<const util::Config&>(params).get("projection",config_proj) )
    projection_.reset( projection::Projection::create(config_proj) );
  else
    projection_.reset( projection::Projection::create() );

  // mandatory parameters: pl, latitudes
  std::vector<long> nx;
  std::vector<double> y;

  if (!params.get("nx", nx))
      throw eckit::BadParameter("nx missing in Params",Here());
  if (!params.get("y", y))
      throw eckit::BadParameter("y missing in Params",Here());

  ASSERT(y.size());
  ASSERT(y.size() == y.size());
  const size_t ny = nx.size();

  // optional parameters: N identifier, longitude limits (lon_min, lon_max)
  std::vector<double> xmin(ny,std::numeric_limits<double>::max());
  std::vector<double> xmax(ny,-std::numeric_limits<double>::max());

  params.get("N", Structured::N_);
  if (params.has("xmin"))
      params.get("xmin", xmin);
  if (params.has("xmax"))
      params.get("xmax", xmax);

  ASSERT(xmin.size() == ny);
  ASSERT(xmax.size() == ny);

  setup_cropped(ny, y.data(), nx.data(), xmin.data(), xmax.data(), *domain_ );

  // common (base) class setup
  //Structured::setup(latitudes.size(), latitudes.data(), pl.data(), lonmin.data(), lonmax.data());
}



CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long pl[]) :
    Structured()
{
    ASSERT(nlat);

    periodic_x_ = true;
    periodic_y_ = false;

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

