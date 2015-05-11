/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/parser/JSON.h"
#include "atlas/grids/ReducedGrid.h"
#include "atlas/grids/LonLatGrid.h"
#include "atlas/ErrorHandling.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace trans {

Trans::Trans(const grids::ReducedGrid& g, const Trans::Options& p)
{
  int nsmax = 0;
  ctor_rgg(g.nlat(),g.npts_per_lat().data(), nsmax, p);
}

Trans::Trans(const int N, const Trans::Options& p)
{
  int nsmax = 0;
  std::vector<int> npts_per_lat(2*N,4*N);
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::Trans(const grids::ReducedGrid& g, const int nsmax, const Trans::Options& p )
{
  const grids::LonLatGrid* lonlat = dynamic_cast<const grids::LonLatGrid*>(&g);
  if( lonlat )
    ctor_lonlat( lonlat->nlon(), lonlat->nlat(), nsmax, p );
  else
    ctor_rgg(g.nlat(),g.npts_per_lat().data(), nsmax, p);
}


Trans::Trans(const int N, const int nsmax, const Trans::Options& p)
{
  std::vector<int> npts_per_lat(2*N,4*N);
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::Trans( const std::vector<int>& npts_per_lat, const int nsmax, const Trans::Options& p )
{
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::~Trans()
{
  ::trans_delete(&trans_);
}

void Trans::ctor_rgg(const int ndgl, const int nloen[], int nsmax, const Trans::Options& p )
{
  ::trans_new(&trans_);
  ::trans_set_resol(&trans_,ndgl,nloen);
  ::trans_set_trunc(&trans_,nsmax);

  trans_.fft = p.fft();
  trans_.lsplit = p.split_latitudes();
  trans_.luseflt = p.flt();

  ::trans_setup(&trans_);
}

void Trans::ctor_lonlat(const int nlon, const int nlat, int nsmax, const Trans::Options& p )
{
  ::trans_new(&trans_);
  ::trans_set_resol_lonlat(&trans_,nlon,nlat);
  ::trans_set_trunc(&trans_,nsmax);

  trans_.fft = p.fft();
  trans_.lsplit = p.split_latitudes();
  trans_.luseflt = p.flt();

  ::trans_setup(&trans_);
}

Trans::Options::Options() : eckit::Properties()
{
  set_split_latitudes(true);
  set_fft(FFTW);
  set_flt(false);
}

void Trans::Options::print( std::ostream& s) const
{
  eckit::JSON js(s);
  js.precision(16);
  js << *this;
}

void Trans::Options::set_fft( FFT fft )
{
  if( fft == FFTW )
  {
    set( "fft", "FFTW" );
  }
  else if( fft == FFT992 )
  {
    set( "fft", "FFT992" );
  }
  else
  {
    NOTIMP;
  }
}

void Trans::Options::set_split_latitudes( bool split )
{
  set("split_latitudes",split);
}

void Trans::Options::set_flt( bool flt )
{
  set("flt",flt);
}

bool Trans::Options::split_latitudes() const
{
  return get("split_latitudes");
}

FFT Trans::Options::fft() const
{
  std::string fftstr = get( "fft" );
  if( fftstr == "FFTW" )
    return FFTW;
  else if( fftstr == "FFT992" )
    return FFT992;
  else
    NOTIMP;
  return FFTW;
}

bool Trans::Options::flt() const
{
  return get("flt");
}

eckit::Params::value_t get( const Trans::Options& p, const eckit::Params::key_t& key )
{
  return p.get(key);
}

void print( const Trans::Options& p, std::ostream& s )
{
  p.print(s);
}

void encode( const Trans::Options& p, eckit::Stream& s )
{
  s << p;
}


Trans* atlas__Trans__new (grids::ReducedGrid* grid)
{
  Trans* trans;
  ATLAS_ERROR_HANDLING(
    ASSERT( grid != NULL );
    trans = new Trans(*grid);
  );
  return trans;
}

void atlas__Trans__delete (Trans* trans)
{
  ATLAS_ERROR_HANDLING( delete trans );
}

int atlas__Trans__handle (Trans* trans)
{
  ATLAS_ERROR_HANDLING( return trans->handle() );
}

}
}
