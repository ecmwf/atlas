/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo>  // std::bad_cast
#include <string>

#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"

#include "atlas/grids/ReducedGrid.h"
#include "atlas/util/Debug.h"

using eckit::Factory;
using eckit::MD5;
using eckit::BadParameter;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid, ReducedGrid, ReducedGrid::grid_type_str());

ReducedGrid* ReducedGrid::create(const eckit::Parametrisation& p) {

  ReducedGrid* grid = dynamic_cast<ReducedGrid*>(Grid::create(p));
  if (!grid) throw BadParameter("Grid is not a reduced grid", Here());
  return grid;

}

ReducedGrid* ReducedGrid::create(const std::string& uid)
{
  ReducedGrid* grid = dynamic_cast<ReducedGrid*>( Grid::create(uid) );
  if( !grid )
    throw BadParameter("Grid "+uid+" is not a reduced grid",Here());
  return grid;
}

//ReducedGrid* ReducedGrid::create(const eckit::Properties& g)
//{
//  ReducedGrid* grid = dynamic_cast<ReducedGrid*>( Grid::create(g) );
//  if( !grid )
//    throw BadParameter("Grid is not a reduced grid",Here());
//  return grid;
//}

std::string ReducedGrid::className() { return "atlas.ReducedGrid"; }

ReducedGrid::ReducedGrid(const Domain& d) : Grid(d), N_(0)
{
}

ReducedGrid::ReducedGrid(const eckit::Parametrisation& params) : N_(0)
{
  setup(params);

  if( ! params.get("grid_type",grid_type_) ) throw BadParameter("grid_type missing in Params",Here());
  if( ! params.get("shortName",shortName_) ) throw BadParameter("shortName missing in Params",Here());
  //if( ! params.has("hash") ) throw BadParameter("hash missing in Params",Here());
}

void ReducedGrid::setup(const eckit::Parametrisation& params)
{
  eckit::ValueList list;

  std::vector<long> npts_per_lat;
  std::vector<double> latitudes;

  if( ! params.get("npts_per_lat",npts_per_lat) ) throw BadParameter("npts_per_lat missing in Params",Here());
  if( ! params.get("latitudes",latitudes) ) throw BadParameter("latitudes missing in Params",Here());

  params.get("N",N_);

  setup(latitudes.size(),latitudes.data(),npts_per_lat.data());
}

ReducedGrid::ReducedGrid(const std::vector<double>& _lats, const std::vector<size_t>& _nlons, const Domain& d)
  : Grid(d)
{
  size_t nlat = _nlons.size();
  std::vector<long> nlons(nlat);
  for( int j=0; j<nlat; ++j )
  {
    nlons[j] = static_cast<long>(nlons[j]);
  }
  setup(nlat,_lats.data(),nlons.data());
}

ReducedGrid::ReducedGrid(size_t nlat, const double lats[], const long nlons[], const Domain& d)
  : Grid(d)
{
  setup(nlat,lats,nlons);
}

void ReducedGrid::setup( const size_t nlat, const double lats[], const long nlons[], const double lonmin[], const double lonmax[] )
{
  ASSERT(nlat > 1);  // can't have a grid with just one latitude

  nlons_.assign(nlons,nlons+nlat);

  lat_.assign(lats,lats+nlat);

  lonmin_.assign(lonmin,lonmin+nlat);
  lonmax_.assign(lonmax,lonmax+nlat);

  npts_ = 0;
  nlonmax_ = 0;
  double lon_min(1000), lon_max(-1000);

  for( int jlat=0; jlat<nlat; ++jlat )
  {
    //ASSERT( nlon(jlat) > 1 ); // can't have grid with just one longitude
    nlonmax_ = std::max(nlon(jlat),nlonmax_);

    lon_min = std::min(lon_min,lonmin_[jlat]);
    lon_max = std::max(lon_max,lonmax_[jlat]);

    npts_ += nlons_[jlat];
  }

  bounding_box_ = BoundBox(lat_[0]/*north*/, lat_[nlat-1]/*south*/, lon_max/*east*/, lon_min/*west*/ );
}


void ReducedGrid::setup( const size_t nlat, const double lats[], const long nlons[] )
{
  std::vector<double> lonmin(nlat,0.);
  std::vector<double> lonmax(nlat);
  for( int jlat=0; jlat<nlat; ++jlat )
  {
    if( nlons[jlat] )
      lonmax[jlat] = 360.-360./static_cast<double>(nlons[jlat]);
    else
      lonmax[jlat] = 0.;
  }
  setup(nlat,lats,nlons,lonmin.data(),lonmax.data());
}

void ReducedGrid::setup_lat_hemisphere(const size_t N, const double lat[], const long lon[], const AngleUnit unit)
{
  std::vector<long> nlons(2*N);
  std::copy( lon, lon+N, nlons.begin() );
  std::reverse_copy( lon, lon+N, nlons.begin()+N );
  std::vector<double> lats(2*N);
  std::copy( lat, lat+N, lats.begin() );
  std::reverse_copy( lat, lat+N, lats.begin()+N );
  double convert = (unit == RAD ? 180.*M_1_PI : 1.);
  for( int j=0; j<N; ++j )
    lats[j] *= convert;
  for( int j=N; j<2*N; ++j )
    lats[j] *= -convert;
  setup(2*N,lats.data(),nlons.data());
}

size_t ReducedGrid::N() const
{
  if( N_==0 )
  {
    throw eckit::Exception("N cannot be returned because grid of type "+gridType()+
                           " is not based on a global grid.", Here() );
  }
  return N_;
}

BoundBox ReducedGrid::boundingBox() const
{
  return bounding_box_;
}

size_t ReducedGrid::npts() const { return npts_; }

void ReducedGrid::lonlat( std::vector<Point>& pts ) const
{
  pts.resize(npts());
  int c(0);
  for( int jlat=0; jlat<nlat(); ++jlat )
  {
    double y = lat(jlat);
    for( int jlon=0; jlon<nlon(jlat); ++jlon )
    {
      pts[c++].assign(lon(jlat,jlon),y);
    }
  }
}

std::string ReducedGrid::gridType() const
{
  return grid_type_;
}

eckit::Properties ReducedGrid::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type",gridType());

  grid_spec.set("nlat",nlat());

  grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));
  grid_spec.set("npts_per_lat",eckit::makeVectorValue(npts_per_lat()));

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  if( N_ != 0 )
    grid_spec.set("N", N_ );

  return grid_spec;
}

size_t ReducedGrid::nlat() const
{
  return lat_.size();
}

size_t ReducedGrid::nlon(size_t jlat) const
{
  return nlons_[jlat];
}

size_t ReducedGrid::nlonmax() const
{
  return nlonmax_;
}

const std::vector<int>&  ReducedGrid::npts_per_lat() const
{
  return nlons_;
}

double ReducedGrid::lon(const size_t jlat, const size_t jlon) const
{
  return lonmin_[jlat] + (double)jlon * (lonmax_[jlat]-lonmin_[jlat]) / ( (double)nlon(jlat) - 1. );
}

double ReducedGrid::lat(const size_t jlat) const
{
  return lat_[jlat];
}

void ReducedGrid::lonlat( const size_t jlon, const size_t jlat, double crd[] ) const
{
  crd[0] = lon(jlat,jlon);
  crd[1] = lat(jlat);
}

size_t ReducedGrid::copyLonLatMemory(double* pts, size_t size) const
{
    size_t sizePts = 2*npts();

    ASSERT(size >= sizePts);

    for(size_t c = 0, jlat=0; jlat<nlat(); ++jlat )
    {
      double y = lat(jlat);
      for( size_t jlon=0; jlon<nlon(jlat); ++jlon )
      {
        pts[c++] = lon(jlat,jlon);
        pts[c++] = y;
      }
    }
    return sizePts;
}

void ReducedGrid::print(std::ostream& os) const
{
    os << "ReducedGrid(Name:" << shortName() << ")";
}

const std::vector<double>& ReducedGrid::latitudes() const
{
  return lat_;
}

std::string ReducedGrid::shortName() const {
  ASSERT(!shortName_.empty());
  return shortName_;
}

void ReducedGrid::hash(eckit::MD5& md5) const {
  // Through inheritance the grid_type_str() might differ while still being same grid
      //md5.add(grid_type_str());

  md5.add(latitudes().data(),    sizeof(double)*latitudes().size());
  md5.add(npts_per_lat().data(), sizeof(int)*npts_per_lat().size());
  bounding_box_.hash(md5);
}

//----------------------------------------------------------------------------------------------------------------------

int  atlas__ReducedGrid__nlat(ReducedGrid* This)
{
  return This->nlat();
}

void atlas__ReducedGrid__nlon(ReducedGrid* This, const int* &nlons, int &size)
{
  nlons = This->npts_per_lat().data();
  size  = This->npts_per_lat().size();
}

int atlas__ReducedGrid__npts(ReducedGrid* This)
{
  return This->npts();
}

double atlas__ReducedGrid__lon(ReducedGrid* This,int jlat,int jlon)
{
  return This->lon(jlat,jlon);
}

double atlas__ReducedGrid__lat(ReducedGrid* This,int jlat)
{
  return This->lat(jlat);
}

void atlas__ReducedGrid__latitudes(ReducedGrid* This, const double* &lat, int &size)
{
  lat  = This->latitudes().data();
  size = This->latitudes().size();
}

} // namespace grids
} // namespace atlas
