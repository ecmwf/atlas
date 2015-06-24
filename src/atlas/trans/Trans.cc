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
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/trans/Trans.h"
#include "atlas/FieldSet.h"
#include "atlas/FunctionSpace.h"
#include "atlas/util/Array.h"
#include "atlas/util/Bitflags.h"
#include "eckit/exception/Exceptions.h"

#define TRANS_CHECK( CALL ) do {\
  int errcode = CALL;\
  if( errcode != TRANS_SUCCESS) {\
    std::stringstream msg; msg << "ERROR: " << #CALL << " failed:\n"\
                               <<  ::trans_error_msg(errcode);\
    throw eckit::Exception(msg.str(),Here());\
  }\
} while(0)

using atlas::util::Topology;

namespace atlas {
namespace trans {

Trans::Trans(const Grid& grid, const Trans::Options& p)
{
  const grids::ReducedGrid* reduced = dynamic_cast<const grids::ReducedGrid*>(&grid);
  if( !reduced )
    throw eckit::BadCast("Grid is not a grids::ReducedGrid type. Cannot partition using IFS trans",Here());
  size_t nsmax = 0;
  ctor_rgg(reduced->nlat(),reduced->npts_per_lat().data(), nsmax, p);
}

Trans::Trans(const size_t N, const Trans::Options& p)
{
  size_t nsmax = 0;
  std::vector<int> npts_per_lat(2*N,4*N);
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::Trans(const Grid& grid, const size_t nsmax, const Trans::Options& p )
{
  const grids::ReducedGrid* reduced = dynamic_cast<const grids::ReducedGrid*>(&grid);
  if( !reduced )
    throw eckit::BadCast("Grid is not a grids::ReducedGrid type. Cannot partition using IFS trans",Here());

  const grids::LonLatGrid* lonlat = dynamic_cast<const grids::LonLatGrid*>(reduced);
  if( lonlat )
    ctor_lonlat( lonlat->nlon(), lonlat->nlat(), nsmax, p );
  else
    ctor_rgg(reduced->nlat(),reduced->npts_per_lat().data(), nsmax, p);
}


Trans::Trans(const size_t N, const size_t nsmax, const Trans::Options& p)
{
  std::vector<int> npts_per_lat(2*N,4*N);
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::Trans( const std::vector<int>& npts_per_lat, const size_t nsmax, const Trans::Options& p )
{
  ctor_rgg(npts_per_lat.size(),npts_per_lat.data(), nsmax, p);
}

Trans::~Trans()
{
  ::trans_delete(&trans_);
}

void Trans::ctor_rgg(const size_t ndgl, const int nloen[], size_t nsmax, const Trans::Options& p )
{
  TRANS_CHECK(::trans_new(&trans_));
  TRANS_CHECK(::trans_set_resol(&trans_,ndgl,nloen));
  TRANS_CHECK(::trans_set_trunc(&trans_,nsmax));
  TRANS_CHECK(::trans_set_cache(&trans_,p.cache(),p.cachesize()));

  if( !p.read().empty() )
  {
    if( eckit::PathName(p.read()).exists() )
    {
      std::stringstream msg; msg << "File " << p.read() << "doesn't exist";
      throw eckit::CantOpenFile(msg.str(),Here());
    }
    TRANS_CHECK(::trans_set_read(&trans_,p.read().c_str()));
  }
  if( !p.write().empty() )
    TRANS_CHECK(::trans_set_write(&trans_,p.write().c_str()));

  trans_.fft = p.fft();
  trans_.lsplit = p.split_latitudes();
  trans_.flt = p.flt();

  TRANS_CHECK(::trans_setup(&trans_));
}

void Trans::ctor_lonlat(const size_t nlon, const size_t nlat, size_t nsmax, const Trans::Options& p )
{
  TRANS_CHECK(::trans_new(&trans_));
  TRANS_CHECK(::trans_set_resol_lonlat(&trans_,nlon,nlat));
  TRANS_CHECK(::trans_set_trunc(&trans_,nsmax));
  TRANS_CHECK(::trans_set_cache(&trans_,p.cache(),p.cachesize()));

  if( ! p.read().empty() )
  {
    if( eckit::PathName(p.read()).exists() )
    {
      std::stringstream msg; msg << "File " << p.read() << "doesn't exist";
      throw eckit::CantOpenFile(msg.str(),Here());
    }
    TRANS_CHECK(::trans_set_read(&trans_,p.read().c_str()));
  }
  if( !p.write().empty() )
    TRANS_CHECK(::trans_set_write(&trans_,p.write().c_str()));

  trans_.fft = p.fft();
  trans_.lsplit = p.split_latitudes();
  trans_.flt = p.flt();

  TRANS_CHECK(::trans_setup(&trans_));
}

Trans::Options::Options() : eckit::Properties()
{
  set_cache(NULL,0);
  set_split_latitudes(true);
  set_fft(FFTW);
  set_flt(false);
}

void Trans::Options::set_cache(const void* buffer, size_t size)
{
  cacheptr_=buffer;
  cachesize_=size;
}

const void* Trans::Options::cache() const
{
  return cacheptr_;
}

size_t Trans::Options::cachesize() const
{
  return cachesize_;
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


void Trans::Options::set_read(const std::string& file)
{
  set("read",file);
}

std::string Trans::Options::read() const
{
  if( has("read") )
    return get("read");
  else
    return std::string();
}

void Trans::Options::set_write(const std::string& file)
{
  set("write",file);
}

std::string Trans::Options::write() const
{
  if( has("write") )
    return get("write");
  else
    return std::string();
}

eckit::Params::value_t getValue( const Trans::Options& p, const eckit::Params::key_t& key )
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


// --------------------------------------------------------------------------------------------



void Trans::dirtrans(const Field& gpfield, Field& spfield, const TransParameters& context) const
{
  FieldSet gpfields; gpfields.add_field(Field::Ptr( const_cast<Field*>( &gpfield )) );
  FieldSet spfields; spfields.add_field(spfield.self());
  dirtrans(gpfields,spfields,context);
}


// --------------------------------------------------------------------------------------------


void Trans::dirtrans(const FieldSet& gpfields, FieldSet& spfields, const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  FunctionSpace::Ptr gp;
  for( int jfld=0; jfld<gpfields.size(); ++jfld )
  {
    const Field& f = gpfields[jfld];
    nfld += f.nb_vars();

    if( jfld == 0 ) { gp = FunctionSpace::Ptr( &f.function_space() ); }
    else {
      if( &f.function_space() != gp.get() )
        throw eckit::SeriousBug("dirtrans: fields within gridpoint fieldset don't match");
    }
  }

  int trans_spnfld(0);
  FunctionSpace::Ptr sp;
  for( int jfld=0; jfld<spfields.size(); ++jfld )
  {
    const Field& f = spfields[jfld];
    trans_spnfld += f.nb_vars();
    if( jfld == 0 ) { sp = FunctionSpace::Ptr( &f.function_space() ); }
    else {
      if( &f.function_space() != sp.get() )
        throw eckit::SeriousBug("dirtrans: fields within spectral fieldset don't match");
    }
  }

  if( nfld != trans_spnfld )
    throw eckit::SeriousBug("dirtrans: different number of gridpoint fields than spectral fields");

  if( sp->shape(0) != nspec2() )
    throw eckit::SeriousBug("dirtrans: spectral fields have wrong dimension");

  // Arrays Trans expects
  Array<double> rgp(nfld,ngptot());
  Array<double> rspec(nspec2(),nfld);

  ArrayView<double,2> rgpview (rgp);
  ArrayView<double,2> rspecview (rspec);


  // Pack gridpoints
  {
    ArrayView<int,1> flags  ( gp->field( "flags" ) );

    int f=0;
    for( int jfld=0; jfld<gpfields.size(); ++jfld )
    {
      const ArrayView<double,2> gpfield ( gpfields[jfld] );
      const int nvars = gpfield.shape(1);

      for( int jvar=0; jvar<nvars; ++jvar )
      {
        int n=0;
        for( int jnode=0; jnode<gp->shape(0); ++jnode )
        {
          bool ghost = Topology::check(flags(jnode),Topology::GHOST);
          if( !ghost )
          {
            rgpview(f,n) = gpfield(jnode,jvar);
            ++n;
          }
        }
        ASSERT( n == ngptot() );
        ++f;
      }
    }
  }

  // Do transform
  {
    struct ::DirTrans_t transform = ::new_dirtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data();
    transform.rspscalar  = rspec.data();

    TRANS_CHECK(::trans_dirtrans(&transform));
  }

  // Unpack the spectral fields
  {
    int f=0;
    for( int jfld=0; jfld<spfields.size(); ++jfld )
    {
      ArrayView<double,2> spfield ( spfields[jfld] );
      const int nvars = spfield.shape(1);

      for( int jvar=0; jvar<nvars; ++jvar )
      {
        for( int jwave=0; jwave<nspec2(); ++jwave )
        {
          spfield(jwave,jvar) = rspecview(jwave,f);
        }
        ++f;
      }
    }
  }

}


// --------------------------------------------------------------------------------------------


void Trans::invtrans(const Field& spfield, Field& gpfield, const TransParameters& context) const
{
  FieldSet spfields; spfields.add_field(Field::Ptr( const_cast<Field*>( &spfield )) );
  FieldSet gpfields; gpfields.add_field(gpfield.self());
  invtrans(spfields,gpfields,context);
}


// --------------------------------------------------------------------------------------------


void Trans::invtrans(const FieldSet& spfields, FieldSet& gpfields, const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  FunctionSpace::Ptr gp;
  for( int jfld=0; jfld<gpfields.size(); ++jfld )
  {
    const Field& f = gpfields[jfld];
    nfld += f.nb_vars();

    if( jfld == 0 ) { gp = FunctionSpace::Ptr( &f.function_space() ); }
    else {
      if( &f.function_space() != gp.get() )
        throw eckit::SeriousBug("invtrans: fields within gridpoint fieldset don't match",Here());
    }
  }

  int nb_spectral_fields(0);
  FunctionSpace::Ptr sp;
  for( int jfld=0; jfld<spfields.size(); ++jfld )
  {
    const Field& f = spfields[jfld];
    nb_spectral_fields += f.nb_vars();
    if( jfld == 0 ) { sp = FunctionSpace::Ptr( &f.function_space() ); }
    else {
      if( &f.function_space() != sp.get() )
        throw eckit::SeriousBug("invtrans: fields within spectral fieldset don't match",Here());
    }
  }

  if( nfld != nb_spectral_fields )
    throw eckit::SeriousBug("invtrans: different number of gridpoint fields than spectral fields",Here());

  if( sp->shape(0) != nspec2() ) {
    std::stringstream msg;
    msg << "invtrans: spectral fields have wrong dimension: nspec2 "<<sp->shape(0)<<" should be "<<nspec2();
    throw eckit::SeriousBug(msg.str(),Here());
  }
  // Arrays Trans expects
  Array<double> rgp(nfld,ngptot());
  Array<double> rspec(nspec2(),nfld);

  ArrayView<double,2> rgpview (rgp);
  ArrayView<double,2> rspecview (rspec);

  // Pack spectral fields
  {
    int f=0;
    for( int jfld=0; jfld<spfields.size(); ++jfld )
    {
      const ArrayView<double,2> field ( spfields[jfld] );
      const int nvars = field.shape(1);

      for( int jvar=0; jvar<nvars; ++jvar )
      {
        for( int jwave=0; jwave<nspec2(); ++jwave )
        {
          rspecview(jwave,f) = field(jwave,jvar);
        }
        ++f;
      }
    }
  }

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data();
    transform.rspscalar  = rspec.data();

    TRANS_CHECK(::trans_invtrans(&transform));
  }

  // Unpack the gridpoint fields
  {
    ArrayView<int,1> flags  ( gp->field( "flags" ) );

    int f=0;
    for( int jfld=0; jfld<gpfields.size(); ++jfld )
    {
      ArrayView<double,2> field ( gpfields[jfld] );
      const int nvars = field.shape(1);

      for( int jvar=0; jvar<nvars; ++jvar )
      {
        int n=0;
        for( int jnode=0; jnode<gp->shape(0); ++jnode )
        {
          bool ghost = Topology::check(flags(jnode),Topology::GHOST);
          if( !ghost )
          {
            field(jnode,jvar) = rgpview(f,n);
            ++n;
          }
        }
        ASSERT( n == ngptot() );
        ++f;
      }
    }
  }

}


Trans* atlas__Trans__new (const Grid* grid, int nsmax)
{
  Trans* trans;
  ATLAS_ERROR_HANDLING(
    ASSERT( grid != NULL );
    trans = new Trans(*grid,nsmax);
  );
  return trans;
}

void atlas__Trans__delete (Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( delete This );
}

int atlas__Trans__handle (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->handle() );
  return 0;
}

///////////////////////////////////////////////////////////////////////////////


void Trans::distspec( const int nb_fields, const int origin[], const double global_spectra[], double spectra[] ) const
{
  struct ::DistSpec_t args = new_distspec(&trans_);
    args.nfld = nb_fields;
    args.rspecg = global_spectra;
    args.nfrom = origin;
    args.rspec = spectra;
  if( int errcode = ::trans_distspec(&args) != TRANS_SUCCESS )
    throw eckit::Exception("distspec failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::gathspec( const int nb_fields, const int destination[], const double spectra[], double global_spectra[] ) const
{
  struct ::GathSpec_t args = new_gathspec(&trans_);
    args.nfld = nb_fields;
    args.rspecg = global_spectra;
    args.nto = destination;
    args.rspec = spectra;
  if( int errcode = ::trans_gathspec(&args) != TRANS_SUCCESS )
    throw eckit::Exception("gathspec failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::distgrid( const int nb_fields, const int origin[], const double global_fields[], double fields[] ) const
{
  struct ::DistGrid_t args = new_distgrid(&trans_);
    args.nfld  = nb_fields;
    args.nfrom = origin;
    args.rgpg  = global_fields;
    args.rgp   = fields;
  if( int errcode = ::trans_distgrid(&args) != TRANS_SUCCESS )
    throw eckit::Exception("distgrid failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::gathgrid( const int nb_fields, const int destination[], const double fields[], double global_fields[] ) const
{
  struct ::GathGrid_t args = new_gathgrid(&trans_);
    args.nfld = nb_fields;
    args.nto  = destination;
    args.rgp  = fields;
    args.rgpg = global_fields;
  if( int errcode = ::trans_gathgrid(&args) != TRANS_SUCCESS )
    throw eckit::Exception("gathgrid failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::invtrans( const int nb_fields, const double scalar_spectra[], double scalar_fields[] ) const
{
  struct ::InvTrans_t args = new_invtrans(&trans_);
    args.nscalar = nb_fields;
    args.rspscalar = scalar_spectra;
    args.rgp = scalar_fields;
  if( int errcode = ::trans_invtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("invtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::invtrans( const int nb_fields, const double vorticity_spectra[], const double divergence_spectra[], double wind_fields[] ) const
{
  struct ::InvTrans_t args = new_invtrans(&trans_);
    args.nvordiv = nb_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp = wind_fields;
  if( int errcode = ::trans_invtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("invtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[] ) const
{
  struct ::DirTrans_t args = new_dirtrans(&trans_);
    args.nscalar = nb_fields;
    args.rgp = scalar_fields;
    args.rspscalar = scalar_spectra;
  if( int errcode = ::trans_dirtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("dirtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void Trans::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[], double divergence_spectra[] ) const
{
  struct ::DirTrans_t args = new_dirtrans(&trans_);
    args.nvordiv = nb_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp    = wind_fields;
  if( int errcode = ::trans_dirtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("dirtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

///////////////////////////////////////////////////////////////////////////////

void atlas__Trans__distspec( const Trans* t, int nb_fields, int origin[], double global_spectra[], double spectra[] )
{
  return t->distspec(nb_fields,origin,global_spectra,spectra);
}

void atlas__Trans__gathspec( const Trans* t, int nb_fields, int destination[], double spectra[], double global_spectra[] )
{
  return t->gathspec(nb_fields,destination,spectra,global_spectra);
}

void atlas__Trans__distgrid( const Trans* t, int nb_fields, int origin[], double global_fields[], double fields[] )
{
  return t->distgrid(nb_fields,origin,global_fields,fields);
}

void atlas__Trans__gathgrid( const Trans* t, int nb_fields, int destination[], double fields[], double global_fields[] )
{
  return t->gathgrid(nb_fields,destination,fields,global_fields);
}

void atlas__Trans__invtrans_scalar( const Trans* t, int nb_fields, double scalar_spectra[], double scalar_fields[] )
{
  return t->invtrans(nb_fields,scalar_spectra,scalar_fields);
}

void atlas__Trans__invtrans_vordiv2wind( const Trans* t, int nb_fields, double vorticity_spectra[], double divergence_spectra[], double wind_fields[] )
{
  return t->invtrans(nb_fields,vorticity_spectra,divergence_spectra,wind_fields);
}

void atlas__Trans__dirtrans_scalar( const Trans* t, int nb_fields, double scalar_fields[], double scalar_spectra[] )
{
  return t->dirtrans(nb_fields,scalar_fields,scalar_spectra);
}

void atlas__Trans__dirtrans_wind2vordiv( const Trans* t, int nb_fields, double wind_fields[], double vorticity_spectra[], double divergence_spectra[] )
{
  return t->dirtrans(nb_fields,wind_fields,vorticity_spectra,divergence_spectra);
}

int atlas__Trans__nproc (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nproc() );
  return 0;
}

int atlas__Trans__myproc (const Trans* This, int proc0)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->myproc(proc0) );
  return 0;
}

int atlas__Trans__ndgl (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->ndgl() );
  return 0;
}

int atlas__Trans__nsmax (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nsmax() );
  return 0;
}

int atlas__Trans__ngptot (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->ngptot() );
  return 0;
}

int atlas__Trans__ngptotg (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->ngptotg() );
  return 0;
}

int atlas__Trans__ngptotmx (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->ngptotmx() );
  return 0;
}

int atlas__Trans__nspec (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nspec() );
  return 0;
}

int atlas__Trans__nspec2 (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nspec2() );
  return 0;
}

int atlas__Trans__nspec2g (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nspec2g() );
  return 0;
}

int atlas__Trans__nspec2mx (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nspec2mx() );
  return 0;
}

int atlas__Trans__n_regions_NS (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->n_regions_NS() );
  return 0;
}

int atlas__Trans__n_regions_EW (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->n_regions_EW() );
  return 0;
}

int atlas__Trans__nump (const Trans* This)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nump() );
  return 0;
}

const int* atlas__Trans__nloen(const Trans* This, int& size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nloen(size) );
  return NULL;
}

const int* atlas__Trans__n_regions (const Trans* This, int& size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->n_regions(size) );
  return NULL;
}

const int* atlas__Trans__nfrstlat(const Trans* This, int& size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nfrstlat(size) );
  return NULL;
}

const int* atlas__Trans__nlstlat (const Trans* This, int& size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nlstlat(size) );
  return NULL;
}

const int* atlas__Trans__nptrfrstlat (const Trans* This, int& size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nptrfrstlat(size) );
  return NULL;
}

const int* atlas__Trans__nsta (const Trans* This, int& sizef2, int& sizef1)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nsta(sizef2,sizef1) );
  return NULL;
}

const int* atlas__Trans__nonl (const Trans* This, int& sizef2, int& sizef1)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nonl(sizef2,sizef1) );
  return NULL;
}

const int* atlas__Trans__nmyms (const Trans* This, int &size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nmyms(size) );
  return NULL;
}

const int* atlas__Trans__nasm0 (const Trans* This, int &size)
{
  ASSERT( This != NULL );
  ATLAS_ERROR_HANDLING( return This->nasm0(size) );
  return NULL;
}

void atlas__Trans__dirtrans_fieldset (const Trans* This, const FieldSet* gpfields, FieldSet* spfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING( This->dirtrans(*gpfields,*spfields,*parameters) );
}

void atlas__Trans__invtrans_fieldset (const Trans* This, const FieldSet* spfields, FieldSet* gpfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING( This->invtrans(*spfields,*gpfields,*parameters) );
}

void atlas__Trans__dirtrans_field (const Trans* This, const Field* gpfield, Field* spfield, const TransParameters* parameters)
{
  ASSERT( This != NULL );
  ASSERT( spfield != NULL );
  ASSERT( gpfield != NULL );
  ASSERT( parameters != NULL );
  ATLAS_ERROR_HANDLING( This->dirtrans(*gpfield,*spfield,*parameters) );
}

void atlas__Trans__invtrans_field (const Trans* This, const Field* spfield, Field* gpfield, const TransParameters* parameters)
{
  ASSERT( This != NULL );
  ASSERT( spfield != NULL );
  ASSERT( gpfield != NULL );
  ASSERT( parameters != NULL );
  ATLAS_ERROR_HANDLING( This->invtrans(*spfield,*gpfield,*parameters) );
}

TransParameters* atlas__TransParameters__new ()
{
  return new TransParameters();
}

void atlas__TransParameters__delete (TransParameters* This)
{
  ASSERT( This != NULL );
  delete This;
}

}
}
