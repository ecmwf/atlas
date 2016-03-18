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
#include <string>
#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"
#include "atlas/grid/global/CustomStructured.h"

using eckit::Factory;
using eckit::MD5;
using eckit::BadParameter;

namespace atlas {
namespace grid {
namespace global {

//------------------------------------------------------------------------------

register_BuilderT1(Grid, CustomStructured, CustomStructured::grid_type_str());

std::string CustomStructured::className()
{ 
  return "atlas.global.CustomStructured";
}

CustomStructured::CustomStructured(const eckit::Parametrisation& params)
  : Structured()
{
  setup(params);

  if( ! params.get("grid_type",grid_type_) ) throw BadParameter("grid_type missing in Params",Here());
  if( ! params.get("shortName",shortName_) ) throw BadParameter("shortName missing in Params",Here());
  //if( ! params.has("hash") ) throw BadParameter("hash missing in Params",Here());
}

void CustomStructured::setup(const eckit::Parametrisation& params)
{
  eckit::ValueList list;

  std::vector<long> npts_per_lat;
  std::vector<double> latitudes;

  if( ! params.get("npts_per_lat",npts_per_lat) ) throw BadParameter("npts_per_lat missing in Params",Here());
  if( ! params.get("latitudes",latitudes) ) throw BadParameter("latitudes missing in Params",Here());

  params.get("N",N_);

  Structured::setup(latitudes.size(),latitudes.data(),npts_per_lat.data());
}

CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long nlons[])
  : Structured()
{
  Structured::setup(nlat,lats,nlons);
}

CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long nlons[],
    const double lonmin[] )
  : Structured()
{
  setup(nlat,lats,nlons,lonmin);
}

void CustomStructured::setup(
    const size_t nlat,
    const double lats[],
    const long nlon[],
    const double lonmin[] )
{
  std::vector<double> lonmax(nlat);
  for( size_t jlat=0; jlat<nlat; ++jlat )
  {
    lonmax[jlat] = lonmin[jlat] + 360. - 360./static_cast<double>(nlon[jlat]);
  }
  Structured::setup(nlat,lats,nlon,lonmin,lonmax.data());
}

eckit::Properties CustomStructured::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type",gridType());

  grid_spec.set("nlat",nlat());

  grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));
  grid_spec.set("npts_per_lat",eckit::makeVectorValue(pl()));
  grid_spec.set("first_longitude_per_latitude",eckit::makeVectorValue(lonmin_));

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  if( N_ != 0 )
    grid_spec.set("N", N_ );

  return grid_spec;
}

extern "C" 
{

Structured* atlas__grid__global__CustomStructured_int(size_t nlat, double lats[], int nlon[])
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new CustomStructured(nlat,lats,nlon_vector.data());
}
Structured* atlas__grid__global__CustomStructured_long(size_t nlat, double lats[], long nlon[])
{
  return new CustomStructured(nlat,lats,nlon);
}
}


} // namespace global
} // namespace grid
} // namespace atlas
