/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/global/lonlat/CustomLonLat.h"

using eckit::BadParameter;

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,CustomLonLat,CustomLonLat::grid_type_str());

std::string CustomLonLat::grid_type_str()
{
  return "custom_lonlat";
}

std::string CustomLonLat::className()
{
  return "atlas.grid.global.lonlat.CustomLonLat";
}

void CustomLonLat::set_typeinfo()
{
  std::stringstream s;
  s << "custom_lonlat.N" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

CustomLonLat::CustomLonLat()
  : ReducedGrid()
{
}

CustomLonLat::CustomLonLat( const size_t nlat, const long nlons[], bool poles, const Domain& domain)
 : ReducedGrid(domain)
{
  ReducedGrid::N_ = nlat;
  poles_ = poles;
  setup(nlat,nlons,poles_);
  set_typeinfo();
}

CustomLonLat::CustomLonLat( const eckit::Parametrisation& params )
{
  setup(params);
  set_typeinfo();
}

void CustomLonLat::setup( const eckit::Parametrisation& params )
{
  if( ! params.has("nlat") )         throw BadParameter("N missing in Params",Here());
  if( ! params.has("npts_per_lat") ) throw BadParameter("npts_per_lat missing in Params",Here());

  size_t nlat;
  params.get("nlat",nlat);

  N_ = nlat;
  params.get("N",N_);

  poles_ = defaults::poles();
  params.get("poles",poles_);


  if( params.has("latitudes") )
  {
    ReducedGrid::setup(params);
  }
  else
  {
    std::vector<long> nlons;
    params.get("npts_per_lat",nlons);
    setup(nlat,nlons.data(),poles_);
  }
}

void CustomLonLat::setup( const size_t nlat, const long nlons[], bool poles )
{
  std::vector<double> lats (nlat);

  double delta = domain_.north() - domain_.south();
  double latmax;

  if( poles )
  {
    delta = delta / static_cast<double>(nlat-1);
    latmax = domain_.north();
  }
  else
  {
    delta = delta / static_cast<double>(nlat);
    latmax = domain_.north() - 0.5*delta;
  }

  for(size_t jlat = 0; jlat < nlat; ++jlat)
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }
  ReducedGrid::setup(lats.size(),lats.data(),nlons);
}


eckit::Properties CustomLonLat::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type", grid_type_str() );

  grid_spec.set("N", N() );
  grid_spec.set("nlat", nlat() );

  grid_spec.set("npts_per_lat",eckit::makeVectorValue(npts_per_lat()));

  if( nlat() != N() )
    grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());
  grid_spec.set("poles",poles_);

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas
