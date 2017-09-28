/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/trans/Trans.h"

#include "eckit/parser/JSON.h"
#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"

using Topology = atlas::mesh::Nodes::Topology;
using atlas::mesh::IsGhostNode;
using atlas::functionspace::StructuredColumns;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::Spectral;
using atlas::Field;
using atlas::array::ArrayView;
using atlas::array::make_view;

// anonymous namespace
namespace {

void trans_check(const int code, const char* msg, const eckit::CodeLocation& location) {
  if(code != TRANS_SUCCESS) {
    std::stringstream errmsg;
    errmsg << "atlas::trans ERROR: " << msg << " failed: \n";
    errmsg << ::trans_error_msg(code);
    throw eckit::Exception(errmsg.str(),location);
  }
}

struct PackNodeColumns
{
  ArrayView<double,2>& rgpview_;
  IsGhostNode is_ghost;
  size_t f;

  PackNodeColumns( ArrayView<double,2>& rgpview, const NodeColumns& fs ) :
    rgpview_(rgpview), is_ghost( fs.nodes() ), f(0) {}

  void operator()(const Field& field, int components = 0) {
    switch (field.rank()) {
      case 1:
        pack_1(field,components);
        break;
      case 2:
        pack_2(field,components);
        break;
      case 3:
        pack_3(field,components);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
      break;
    }
  }

  void pack_1(const Field& field, int)
  {
    const ArrayView<double,1> gpfield = make_view<double,1>( field );
    size_t n=0;
    for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
    {
      if( !is_ghost(jnode) )
      {
        rgpview_(f,n) = gpfield(jnode);
        ++n;
      }
    }
    ++f;
  }
  void pack_2(const Field& field, int)
  {
    const ArrayView<double,2> gpfield = make_view<double,2>( field );
    const size_t nvars = gpfield.shape(1);
    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      size_t n=0;
      for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
      {
        if( !is_ghost(jnode) )
        {
          rgpview_(f,n) = gpfield(jnode,jvar);
          ++n;
        }
      }
      ++f;
    }
  }
  void pack_3(const Field& field, int components)
  {
    const ArrayView<double,3> gpfield = make_view<double,3>( field );
    if( not components ) components = gpfield.shape(2);
    ATLAS_DEBUG_VAR( components );
    for( size_t jcomp=0; jcomp<size_t(components); ++jcomp )
    {
      for( size_t jlev=0; jlev<gpfield.shape(1); ++jlev )
      {
        size_t n = 0;
        for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
        {
          if( !is_ghost(jnode) )
          {
            rgpview_(f,n) = gpfield(jnode,jlev,jcomp);
            ++n;
          }
        }
        ++f;
      }
    }
  }
};


struct PackStructuredColumns
{
  ArrayView<double,2>& rgpview_;
  size_t f;

  PackStructuredColumns( ArrayView<double,2>& rgpview ) :
    rgpview_(rgpview), f(0) {}

  void operator()(const Field& field) {
    switch (field.rank()) {
      case 1:
        pack_1(field);
        break;
      case 2:
        pack_2(field);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
        break;
    }
  }

  void pack_1(const Field& field)
  {
    const ArrayView<double,1> gpfield = make_view<double,1>( field );
    size_t n=0;
    for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
    {
      rgpview_(f,n) = gpfield(jnode);
      ++n;
    }
    ++f;
  }
  void pack_2(const Field& field)
  {
    const ArrayView<double,2> gpfield = make_view<double,2>( field );
    const size_t nvars = gpfield.shape(1);
    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      size_t n=0;
      for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
      {
        rgpview_(f,n) = gpfield(jnode,jvar);
        ++n;
      }
      ++f;
    }
  }
};

struct PackSpectral
{
  ArrayView<double,2>& rspecview_;
  size_t f;
  PackSpectral( ArrayView<double,2>& rspecview ) :
    rspecview_(rspecview), f(0) {}

  void operator()(const Field& field) {
    switch (field.rank()) {
      case 1:
        pack_1(field);
        break;
      case 2:
        pack_2(field);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
        break;
    }
  }

  void pack_1(const Field& field)
  {
    const ArrayView<double,1> spfield = make_view<double,1>( field );

    for( size_t jwave=0; jwave<spfield.shape(0); ++jwave )
    {
      rspecview_(jwave,f) = spfield(jwave);
    }
    ++f;
  }
  void pack_2(const Field& field)
  {
    const ArrayView<double,2> spfield = make_view<double,2>( field );

    const size_t nvars = spfield.shape(1);

    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      for( size_t jwave=0; jwave<spfield.shape(0); ++jwave )
      {
        rspecview_(jwave,f) = spfield(jwave,jvar);
      }
      ++f;
    }
  }

};

struct UnpackNodeColumns
{
  const ArrayView<double,2>& rgpview_;
  IsGhostNode is_ghost;
  size_t f;

  UnpackNodeColumns( const ArrayView<double,2>& rgpview, const NodeColumns& fs ) :
    rgpview_(rgpview), is_ghost( fs.nodes() ), f(0) {}

  void operator()(Field& field, int components = 0) {
    switch (field.rank()) {
      case 1:
        unpack_1(field,components);
        break;
      case 2:
        unpack_2(field,components);
        break;
      case 3:
        unpack_3(field,components);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
        break;
    }
  }

  void unpack_1(Field& field, int)
  {
    ArrayView<double,1> gpfield = make_view<double,1>( field );
    size_t n(0);
    for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
    {
      if( !is_ghost(jnode) )
      {
        gpfield(jnode) = rgpview_(f,n);
        ++n;
      }
    }
    ++f;
  }
  void unpack_2(Field& field, int)
  {
    ArrayView<double,2> gpfield = make_view<double,2>( field );
    const size_t nvars = gpfield.shape(1);
    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      int n=0;
      for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
      {
        if( !is_ghost(jnode) )
        {
          gpfield(jnode,jvar) = rgpview_(f,n);
          ++n;
        }
      }
      ++f;
    }
  }
  void unpack_3(Field& field, int components)
  {
    ArrayView<double,3> gpfield = make_view<double,3>( field );
    if( not components ) components = gpfield.shape(2);
    for( size_t jcomp=0; jcomp<size_t(components); ++jcomp )
    {
      for( size_t jlev=0; jlev<gpfield.shape(1); ++jlev )
      {
        size_t n = 0;
        for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
        {
          if( !is_ghost(jnode) )
          {
            gpfield(jnode,jlev,jcomp) = rgpview_(f,n);
            ++n;
          }
        }
        ++f;
      }
    }
  }
};

struct UnpackStructuredColumns
{
  const ArrayView<double,2>& rgpview_;
  size_t f;

  UnpackStructuredColumns( const ArrayView<double,2>& rgpview ) :
    rgpview_(rgpview), f(0) {}

  void operator()(Field& field) {
    switch (field.rank()) {
      case 1:
        unpack_1(field);
        break;
      case 2:
        unpack_2(field);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
        break;
    }
  }

  void unpack_1(Field& field)
  {
    ArrayView<double,1> gpfield = make_view<double,1>( field );
    size_t n=0;
    for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
    {
      gpfield(jnode) = rgpview_(f,n);
      ++n;
    }
    ++f;
  }
  void unpack_2(Field& field)
  {
    ArrayView<double,2> gpfield = make_view<double,2>( field );
    const size_t nvars = gpfield.shape(1);
    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      size_t n=0;
      for( size_t jnode=0; jnode<gpfield.shape(0); ++jnode )
      {
        gpfield(jnode,jvar) = rgpview_(f,n);
        ++n;
      }
      ++f;
    }
  }
};

struct UnpackSpectral
{
  const ArrayView<double,2>& rspecview_;
  size_t f;
  UnpackSpectral( const ArrayView<double,2>& rspecview ) :
    rspecview_(rspecview), f(0) {}

  void operator()(Field& field) {
    switch (field.rank()) {
      case 1:
        unpack_1(field);
        break;
      case 2:
        unpack_2(field);
        break;
      default:
        ATLAS_DEBUG_VAR(field.rank());
        NOTIMP;
        break;
    }
  }

  void unpack_1(Field& field)
  {
    ArrayView<double,1> spfield = make_view<double,1>( field );

    for( size_t jwave=0; jwave<spfield.shape(0); ++jwave )
    {
      spfield(jwave) = rspecview_(jwave,f);
    }
    ++f;
  }
  void unpack_2(Field& field)
  {
    ArrayView<double,2> spfield = make_view<double,2>( field );

    const size_t nvars = spfield.shape(1);

    for( size_t jvar=0; jvar<nvars; ++jvar )
    {
      for( size_t jwave=0; jwave<spfield.shape(0); ++jwave )
      {
        spfield(jwave,jvar) = rspecview_(jwave,f);
      }
      ++f;
    }
  }

};


} // end anonymous namespace

#define TRANS_CHECK( CALL ) trans_check(CALL, #CALL, Here() )


namespace atlas {
namespace trans {

Trans::Trans(const Grid& grid, const Trans::Options& p) :
  grid_(grid) {
  ASSERT( grid.domain().global() );
  ASSERT( not grid.projection() );
  size_t grid_only = -1; // triggers grid-only setup
  ctor(grid,grid_only,p);
}

Trans::Trans(const Grid& grid, const long truncation, const Trans::Options& p ) :
  grid_(grid) {
  ASSERT( grid.domain().global() );
  ASSERT( not grid.projection() );
  ctor( grid,truncation,p);
}

Trans::Trans(const long truncation, const Trans::Options& p )
{
  ctor_spectral_only(truncation,p);
}


Trans::~Trans()
{
  ::trans_delete(&trans_);
}

void Trans::ctor( const Grid& grid, long truncation, const Trans::Options& p ) {

  if( auto gg = grid::GaussianGrid(grid) ) {
    ctor_rgg(gg.ny(), gg.nx().data(), truncation, p);
    return;
  }
  if( auto ll = grid::RegularLonLatGrid(grid) ) {
    if( ll.standard() || ll.shifted() ) {
      ctor_lonlat( ll.nx(), ll.ny(), truncation, p );
      return;
    }
  }
  throw eckit::NotImplemented("Grid type not supported for Spectral Transforms",Here());
}

void Trans::ctor_spectral_only(long truncation, const Trans::Options& p )
{
  TRANS_CHECK(::trans_new(&trans_));
  TRANS_CHECK(::trans_set_trunc(&trans_,truncation));
  TRANS_CHECK(::trans_use_mpi(parallel::mpi::comm().size()>1));
  TRANS_CHECK(::trans_setup(&trans_));
}

void Trans::ctor_rgg(const long nlat, const long pl[], long truncation, const Trans::Options& p )
{
  std::vector<int> nloen(nlat);
  for( long jlat=0; jlat<nlat; ++jlat )
    nloen[jlat] = pl[jlat];
  TRANS_CHECK(::trans_new(&trans_));
  TRANS_CHECK(::trans_set_resol(&trans_,nlat,nloen.data()));
  if( truncation >= 0 )
    TRANS_CHECK(::trans_set_trunc(&trans_,truncation));
  TRANS_CHECK(::trans_set_cache(&trans_,p.cache(),p.cachesize()));

  if( not p.read().empty() )
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

  TRANS_CHECK(::trans_use_mpi(parallel::mpi::comm().size()>1));
  TRANS_CHECK(::trans_setup(&trans_));
}

void Trans::ctor_lonlat(const long nlon, const long nlat, long truncation, const Trans::Options& p )
{
  TRANS_CHECK(::trans_new(&trans_));
  TRANS_CHECK(::trans_set_resol_lonlat(&trans_,nlon,nlat));
  if( truncation >= 0 )
    TRANS_CHECK(::trans_set_trunc(&trans_,truncation));
  TRANS_CHECK(::trans_set_cache(&trans_,p.cache(),p.cachesize()));

  if( not p.read().empty() )
  {
    if( eckit::PathName(p.read()).exists() )
    {
      std::stringstream msg; msg << "File " << p.read() << "doesn't exist";
      throw eckit::CantOpenFile(msg.str(),Here());
    }
    TRANS_CHECK(::trans_set_read(&trans_,p.read().c_str()));
  }
  if( not p.write().empty() )
    TRANS_CHECK(::trans_set_write(&trans_,p.write().c_str()));

  trans_.fft = p.fft();
  trans_.lsplit = p.split_latitudes();
  trans_.flt = p.flt();

  TRANS_CHECK(::trans_use_mpi(parallel::mpi::comm().size()>1));
  TRANS_CHECK(::trans_setup(&trans_));
}

Trans::Options::Options() : util::Config()
{
  set_cache(nullptr,0);
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
  return getBool("split_latitudes");
}

FFT Trans::Options::fft() const
{
  std::string fftstr = getString( "fft" );
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
  return getBool("flt");
}


void Trans::Options::set_read(const std::string& file)
{
  set("read",file);
}

std::string Trans::Options::read() const
{
  return getString("read","");
}

void Trans::Options::set_write(const std::string& file)
{
  set("write",file);
}

std::string Trans::Options::write() const
{
  return getString("write","");
}


// --------------------------------------------------------------------------------------------



void Trans::dirtrans(const functionspace::NodeColumns& gp, const Field& gpfield,
                     const Spectral& sp, Field& spfield, const TransParameters& context) const
{
  FieldSet gpfields; gpfields.add(gpfield);
  FieldSet spfields; spfields.add(spfield);
  dirtrans(gp,gpfields,sp,spfields,context);
}


// --------------------------------------------------------------------------------------------

void Trans::dirtrans(const functionspace::NodeColumns& gp,const FieldSet& gpfields,
                     const Spectral& sp, FieldSet& spfields, const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
  {
    const Field& f = gpfields[jfld];
    nfld += f.stride(0);
  }

  int trans_spnfld(0);
  for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
  {
    const Field& f = spfields[jfld];
    trans_spnfld += f.stride(0);
  }

  if( nfld != trans_spnfld )
  {
    throw eckit::SeriousBug("dirtrans: different number of gridpoint fields than spectral fields",Here());
  }
  // Arrays Trans expects
  array::ArrayT<double> rgp(nfld,ngptot());
  array::ArrayT<double> rspec(nspec2(),nfld);

  array::ArrayView<double,2> rgpview   = array::make_view<double,2>(rgp);
  array::ArrayView<double,2> rspecview = array::make_view<double,2>(rspec);

  // Pack gridpoints
  {
    PackNodeColumns pack(rgpview,gp);
    for( size_t jfld=0; jfld<gpfields.size(); ++jfld )
      pack( gpfields[jfld]);
  }

  // Do transform
  {
    struct ::DirTrans_t transform = ::new_dirtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data<double>();
    transform.rspscalar  = rspec.data<double>();

    TRANS_CHECK( ::trans_dirtrans(&transform) );
  }

  // Unpack the spectral fields
  {
    UnpackSpectral unpack(rspecview);
    for( size_t jfld=0; jfld<spfields.size(); ++jfld )
      unpack(spfields[jfld]);
  }

}

// --------------------------------------------------------------------------------------------

void Trans::dirtrans(
    const Field& gpfield,
          Field& spfield,
    const TransParameters& context) const
{
  ASSERT( gpfield.functionspace() == 0 ||
          functionspace::StructuredColumns(gpfield.functionspace()) );
  ASSERT( spfield.functionspace() == 0 ||
          functionspace::Spectral(spfield.functionspace()) );
  if ( gpfield.stride(0) != spfield.stride(0) )
  {
    throw eckit::SeriousBug("dirtrans: different number of gridpoint fields than spectral fields",Here());
  }
  if ( (int)gpfield.shape(0) != ngptot() )
  {
    throw eckit::SeriousBug("dirtrans: slowest moving index must be ngptot",Here());
  }
  const int nfld = gpfield.stride(0);

  // Do transform
  {
    struct ::DirTrans_t transform = ::new_dirtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = gpfield.data<double>();
    transform.rspscalar  = spfield.data<double>();
    transform.ngpblks    = gpfield.shape(0);
    transform.nproma     = 1;
    TRANS_CHECK( ::trans_dirtrans(&transform) );
  }
}

void Trans::dirtrans(
    const FieldSet& gpfields,
          FieldSet& spfields,
    const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
  {
    const Field& f = gpfields[jfld];
    nfld += f.stride(0);
    ASSERT( f.functionspace() == 0 ||
            functionspace::StructuredColumns(f.functionspace()) );
  }

  int trans_spnfld(0);
  for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
  {
    const Field& f = spfields[jfld];
    trans_spnfld += f.stride(0);
  }

  if( nfld != trans_spnfld )
  {
    throw eckit::SeriousBug("dirtrans: different number of gridpoint fields than spectral fields",Here());
  }
  // Arrays Trans expects
  array::ArrayT<double> rgp(nfld,ngptot());
  array::ArrayT<double> rspec(nspec2(),nfld);

  array::ArrayView<double,2> rgpview   = array::make_view<double,2>(rgp);
  array::ArrayView<double,2> rspecview = array::make_view<double,2>(rspec);

  // Pack gridpoints
  {
    PackStructuredColumns pack(rgpview);
    for( size_t jfld=0; jfld<gpfields.size(); ++jfld )
      pack(gpfields[jfld]);
  }

  // Do transform
  {
    struct ::DirTrans_t transform = ::new_dirtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data<double>();
    transform.rspscalar  = rspec.data<double>();

    TRANS_CHECK( ::trans_dirtrans(&transform) );
  }

  // Unpack the spectral fields
  {
    UnpackSpectral unpack(rspecview);
    for( size_t jfld=0; jfld<spfields.size(); ++jfld )
      unpack(spfields[jfld]);
  }
}

// --------------------------------------------------------------------------------------------

void Trans::invtrans_grad(const Spectral& sp, const Field& spfield,
                          const functionspace::NodeColumns& gp, Field& gradfield) const
{
  FieldSet spfields;   spfields.  add(spfield);
  FieldSet gradfields; gradfields.add(gradfield);
  invtrans_grad(sp,spfields,gp,gradfields);
}

void Trans::invtrans_grad(const Spectral& sp, const FieldSet& spfields,
                          const functionspace::NodeColumns& gp, FieldSet& gradfields) const
{
  // Count total number of fields and do sanity checks
  int nb_gridpoint_field(0);
  for(size_t jfld = 0; jfld < gradfields.size(); ++jfld)
  {
    const Field& f = gradfields[jfld];
    nb_gridpoint_field += f.stride(0);
  }

  int nfld(0);
  for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
  {
    const Field& f = spfields[jfld];
    nfld += f.stride(0);
    ASSERT( std::max<size_t>(1,f.levels()) == f.stride(0) );
  }

  if( nb_gridpoint_field != 2*nfld ) // factor 2 because N-S and E-W derivatives
    throw eckit::SeriousBug("invtrans_grad: different number of gridpoint fields than spectral fields",Here());

  // Arrays Trans expects
  // Allocate space for
  array::ArrayT<double> rgp(3*nfld,ngptot()); // (scalars) + (NS ders) + (EW ders)
  array::ArrayT<double> rspec(nspec2(),nfld);

  array::ArrayView<double,2> rgpview   = array::make_view<double,2>(rgp);
  array::ArrayView<double,2> rspecview = array::make_view<double,2>(rspec);

  // Pack spectral fields
  {
    PackSpectral pack(rspecview);
    for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
      pack(spfields[jfld]);
  }

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nscalar     = nfld;
    transform.rgp         = rgp.data<double>();
    transform.rspscalar   = rspec.data<double>();
    transform.lscalarders = true;

    TRANS_CHECK(::trans_invtrans(&transform));
  }

  // Unpack the gridpoint fields
  {
    mesh::IsGhostNode is_ghost( gp.nodes());
    int f=nfld; // skip to where derivatives start
    for(size_t dim=0; dim<2; ++dim) {
      for(size_t jfld = 0; jfld < gradfields.size(); ++jfld)
      {
        const size_t nlev = std::max<size_t>(1,gradfields[jfld].levels());
        const size_t nb_nodes  = gradfields[jfld].shape(0);

        array::LocalView<double,3> field ( gradfields[jfld].data<double>(),
           array::make_shape(nb_nodes, nlev, 2 ) );

        for( size_t jlev=0; jlev<nlev; ++jlev )
        {
          int n=0;
          for( size_t jnode=0; jnode<nb_nodes; ++jnode )
          {
            if( !is_ghost(jnode) )
            {
              field(jnode,jlev,1-dim) = rgpview(f,n);
              ++n;
            }
          }
          ASSERT( n == ngptot() );
          ++f;
        }
      }
    }
  }
}

// --------------------------------------------------------------------------------------------

void Trans::invtrans(const Spectral& sp, const Field& spfield,
                     const functionspace::NodeColumns& gp, Field& gpfield, const TransParameters& context) const
{
  FieldSet spfields; spfields.add(spfield);
  FieldSet gpfields; gpfields.add(gpfield);
  invtrans(sp,spfields,gp,gpfields,context);
}


// --------------------------------------------------------------------------------------------


void Trans::invtrans(const Spectral& sp, const FieldSet& spfields,
                     const functionspace::NodeColumns& gp, FieldSet& gpfields, const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
  {
    const Field& f = gpfields[jfld];
    nfld += f.stride(0);
  }

  int nb_spectral_fields(0);
  for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
  {
    const Field& f = spfields[jfld];
    nb_spectral_fields += f.stride(0);
  }

  if( nfld != nb_spectral_fields )
    throw eckit::SeriousBug("invtrans: different number of gridpoint fields than spectral fields",Here());

  // Arrays Trans expects
  array::ArrayT<double> rgp(nfld,ngptot());
  array::ArrayT<double> rspec(nspec2(),nfld);

  array::ArrayView<double,2> rgpview   = array::make_view<double,2>(rgp);
  array::ArrayView<double,2> rspecview = array::make_view<double,2>(rspec);

  // Pack spectral fields
  {
    PackSpectral pack(rspecview);
    for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
      pack(spfields[jfld]);
  }

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data<double>();
    transform.rspscalar  = rspec.data<double>();

    TRANS_CHECK(::trans_invtrans(&transform));
  }

  // Unpack the gridpoint fields
  {
    UnpackNodeColumns unpack(rgpview,gp);
    for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
      unpack(gpfields[jfld]);
  }

}

// --------------------------------------------------------------------------------------------


void Trans::invtrans(const  Field& spfield,
                            Field& gpfield,
                     const TransParameters& context) const
{
  ASSERT( gpfield.functionspace() == 0 ||
          functionspace::StructuredColumns( gpfield.functionspace() ) );
  ASSERT( spfield.functionspace() == 0 ||
          functionspace::Spectral( spfield.functionspace() ) );
  if ( gpfield.stride(0) != spfield.stride(0) )
  {
    throw eckit::SeriousBug("dirtrans: different number of gridpoint fields than spectral fields",Here());
  }
  if ( (int)gpfield.shape(0) != ngptot() )
  {
    throw eckit::SeriousBug("dirtrans: slowest moving index must be ngptot",Here());
  }
  const int nfld = gpfield.stride(0);

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = gpfield.data<double>();
    transform.rspscalar  = spfield.data<double>();
    transform.ngpblks    = gpfield.shape(0);
    transform.nproma     = 1;
    TRANS_CHECK( ::trans_invtrans(&transform) );
  }
}


// --------------------------------------------------------------------------------------------


void Trans::invtrans(const  FieldSet& spfields,
                            FieldSet& gpfields,
                     const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  int nfld(0);
  for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
  {
    const Field& f = gpfields[jfld];
    nfld += f.stride(0);
    ASSERT( f.functionspace() == 0 ||
            functionspace::StructuredColumns( f.functionspace() ) );
  }

  int nb_spectral_fields(0);
  for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
  {
    const Field& f = spfields[jfld];
    nb_spectral_fields += f.stride(0);
  }

  if( nfld != nb_spectral_fields ) {
    std::stringstream msg;
    msg << "invtrans: different number of gridpoint fields than spectral fields"
        << "[ " << nfld << " != " << nb_spectral_fields << " ]";
    throw eckit::SeriousBug(msg.str(),Here());
  }

  // Arrays Trans expects
  array::ArrayT<double> rgp(nfld,ngptot());
  array::ArrayT<double> rspec(nspec2(),nfld);

  array::ArrayView<double,2> rgpview   = array::make_view<double,2>(rgp);
  array::ArrayView<double,2> rspecview = array::make_view<double,2>(rspec);

  // Pack spectral fields
  {
    PackSpectral pack(rspecview);
    for(size_t jfld = 0; jfld < spfields.size(); ++jfld)
      pack(spfields[jfld]);
  }

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nscalar    = nfld;
    transform.rgp        = rgp.data<double>();
    transform.rspscalar  = rspec.data<double>();

    TRANS_CHECK(::trans_invtrans(&transform));
  }

  // Unpack the gridpoint fields
  {
    UnpackStructuredColumns unpack(rgpview);
    for(size_t jfld = 0; jfld < gpfields.size(); ++jfld)
      unpack(gpfields[jfld]);
  }
}

// -----------------------------------------------------------------------------------------------

void Trans::dirtrans_wind2vordiv(const functionspace::NodeColumns& gp, const Field& gpwind,
                                 const Spectral& sp, Field& spvor, Field&spdiv,
                                 const TransParameters& context) const
{
  // Count total number of fields and do sanity checks
  size_t nfld = spvor.stride(0);
  if( spdiv.shape(0) != spvor.shape(0) ) throw eckit::SeriousBug("invtrans: vorticity not compatible with divergence.",Here());
  if( spdiv.shape(1) != spvor.shape(1) ) throw eckit::SeriousBug("invtrans: vorticity not compatible with divergence.",Here());
  size_t nwindfld = gpwind.stride(0);
  if (nwindfld != 2*nfld && nwindfld != 3*nfld) throw eckit::SeriousBug("dirtrans: wind field is not compatible with vorticity, divergence.",Here());

  if( spdiv.shape(0) != size_t(nspec2()) ) {
    std::stringstream msg;
    msg << "dirtrans: Spectral vorticity and divergence have wrong dimension: nspec2 "<<spdiv.shape(0)<<" should be "<<nspec2();
    throw eckit::SeriousBug(msg.str(),Here());
  }

  if( spvor.size() == 0 ) throw eckit::SeriousBug("dirtrans: spectral vorticity field is empty.");
  if( spdiv.size() == 0 ) throw eckit::SeriousBug("dirtrans: spectral divergence field is empty.");

  // Arrays Trans expects
  array::ArrayT<double> rgp(2*nfld,size_t(ngptot()));
  array::ArrayView<double,2> rgpview = array::make_view<double,2>(rgp);

  ATLAS_DEBUG_VAR(gpwind.size());
  ATLAS_DEBUG_VAR(rgp.size());
  ATLAS_DEBUG_VAR(gpwind.stride(0));
  ATLAS_DEBUG_VAR(2*nfld);
  ATLAS_DEBUG_VAR(rgp.stride(0));

  ATLAS_DEBUG_HERE();
  // Pack gridpoints
  {
    PackNodeColumns pack( rgpview, gp );
    int wind_components = 2;
    pack(gpwind, wind_components);
  }
  ATLAS_DEBUG_HERE();

  // Do transform
  {
    struct ::DirTrans_t transform = ::new_dirtrans(&trans_);
    transform.nvordiv = nfld;
    transform.rgp     = rgp.data<double>();
    transform.rspvor  = spvor.data<double>();
    transform.rspdiv  = spdiv.data<double>();

    ASSERT( transform.rspvor );
    ASSERT( transform.rspdiv );
    TRANS_CHECK( ::trans_dirtrans(&transform) );
  }
  ATLAS_DEBUG_HERE();

}


void Trans::invtrans_vordiv2wind(const Spectral& sp, const Field& spvor, const Field& spdiv,
                                 const functionspace::NodeColumns& gp, Field& gpwind, const TransParameters&) const
{
  // Count total number of fields and do sanity checks
  size_t nfld = spvor.stride(0);
  if( spdiv.shape(0) != spvor.shape(0) ) throw eckit::SeriousBug("invtrans: vorticity not compatible with divergence.",Here());
  if( spdiv.shape(1) != spvor.shape(1) ) throw eckit::SeriousBug("invtrans: vorticity not compatible with divergence.",Here());
  size_t nwindfld = gpwind.stride(0);
  if (nwindfld != 2*nfld && nwindfld != 3*nfld) throw eckit::SeriousBug("invtrans: wind field is not compatible with vorticity, divergence.",Here());

  if( spdiv.shape(0) != size_t(nspec2()) ) {
    std::stringstream msg;
    msg << "invtrans: Spectral vorticity and divergence have wrong dimension: nspec2 "<<spdiv.shape(0)<<" should be "<<nspec2();
    throw eckit::SeriousBug(msg.str(),Here());
  }

  ASSERT( spvor.rank() == 2 );
  ASSERT( spdiv.rank() == 2 );
  if( spvor.size() == 0 ) throw eckit::SeriousBug("invtrans: spectral vorticity field is empty.");
  if( spdiv.size() == 0 ) throw eckit::SeriousBug("invtrans: spectral divergence field is empty.");

  // Arrays Trans expects
  array::ArrayT<double> rgp(2*nfld,size_t(ngptot()));
  array::ArrayView<double,2> rgpview = array::make_view<double,2>(rgp);

  // Do transform
  {
    struct ::InvTrans_t transform = ::new_invtrans(&trans_);
    transform.nvordiv = nfld;
    transform.rgp     = rgp.data<double>();
    transform.rspvor  = spvor.data<double>();
    transform.rspdiv  = spdiv.data<double>();

    ASSERT( transform.rspvor );
    ASSERT( transform.rspdiv );
    TRANS_CHECK(::trans_invtrans(&transform));
  }

  // Unpack the gridpoint fields
  {
    UnpackNodeColumns unpack( rgpview, gp );
    int wind_components = 2;
    unpack(gpwind,wind_components);
  }

}





Trans* atlas__Trans__new (const Grid::Implementation* grid, int nsmax)
{
  Trans* trans(0);
  ATLAS_ERROR_HANDLING(
    ASSERT( grid );
    trans = new Trans( Grid(grid) ,nsmax);
  );
  return trans;
}

void atlas__Trans__delete (Trans* This)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( delete This );
}

int atlas__Trans__handle (const Trans* This)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING(
    ::Trans_t* t = *This;
    return t->handle;
  );
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
  TRANS_CHECK( ::trans_distspec(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::gathspec( const int nb_fields, const int destination[], const double spectra[], double global_spectra[] ) const
{
  struct ::GathSpec_t args = new_gathspec(&trans_);
    args.nfld = nb_fields;
    args.rspecg = global_spectra;
    args.nto = destination;
    args.rspec = spectra;
  TRANS_CHECK( ::trans_gathspec(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::distgrid( const int nb_fields, const int origin[], const double global_fields[], double fields[] ) const
{
  struct ::DistGrid_t args = new_distgrid(&trans_);
    args.nfld  = nb_fields;
    args.nfrom = origin;
    args.rgpg  = global_fields;
    args.rgp   = fields;
  TRANS_CHECK( ::trans_distgrid(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::gathgrid( const int nb_fields, const int destination[], const double fields[], double global_fields[] ) const
{
  struct ::GathGrid_t args = new_gathgrid(&trans_);
    args.nfld = nb_fields;
    args.nto  = destination;
    args.rgp  = fields;
    args.rgpg = global_fields;
  TRANS_CHECK( ::trans_gathgrid(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::invtrans( const int nb_scalar_fields, const double scalar_spectra[],
                      const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                      double gp_fields[],
                      const TransParameters& context ) const
{
  struct ::InvTrans_t args = new_invtrans(&trans_);
    args.nscalar = nb_scalar_fields;
    args.rspscalar = scalar_spectra;
    args.nvordiv = nb_vordiv_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp = gp_fields;
  if( context.scalar_derivatives() ) {
    args.lscalarders = true;
  }
  if( context.wind_EW_derivatives() ) {
    args.luvder_EW = true;
  }
  if( context.vorticity_divergence_fields() ) {
    args.lvordivgp = true;
  }
  TRANS_CHECK( ::trans_invtrans(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::invtrans(const int nb_scalar_fields, const double scalar_spectra[],
                     double gp_fields[],
                     const TransParameters& context) const
{
  struct ::InvTrans_t args = new_invtrans(&trans_);
    args.nscalar     = nb_scalar_fields;
    args.rspscalar   = scalar_spectra;
    args.rgp         = gp_fields;
  if( context.scalar_derivatives() ) {
    args.lscalarders = true;
  }
  TRANS_CHECK( ::trans_invtrans(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                      double gp_fields[],
                      const TransParameters& context ) const
{
  struct ::InvTrans_t args = new_invtrans(&trans_);
    args.nvordiv = nb_vordiv_fields;
    args.rspvor  = vorticity_spectra;
    args.rspdiv  = divergence_spectra;
    args.rgp     = gp_fields;
  if( context.wind_EW_derivatives() ) {
      args.luvder_EW = true;
    }
  if( context.vorticity_divergence_fields() ) {
    args.lvordivgp = true;
  }
  TRANS_CHECK( ::trans_invtrans(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[] ) const
{
  struct ::DirTrans_t args = new_dirtrans(&trans_);
    args.nscalar = nb_fields;
    args.rgp = scalar_fields;
    args.rspscalar = scalar_spectra;
  TRANS_CHECK( ::trans_dirtrans(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void Trans::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[], double divergence_spectra[] ) const
{
  struct ::DirTrans_t args = new_dirtrans(&trans_);
    args.nvordiv = nb_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp    = wind_fields;
  TRANS_CHECK( ::trans_dirtrans(&args) );
}

///////////////////////////////////////////////////////////////////////////////


void Trans::specnorm( const int nb_fields, const double spectra[], double norms[], int rank ) const
{
  struct ::SpecNorm_t args = new_specnorm(&trans_);
    args.nfld = nb_fields;
    args.rspec = spectra;
    args.rnorm = norms;
    args.nmaster = rank+1;
  TRANS_CHECK( ::trans_specnorm(&args) );
}

///////////////////////////////////////////////////////////////////////////////

void atlas__Trans__distspec( const Trans* t, int nb_fields, int origin[], double global_spectra[], double spectra[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->distspec(nb_fields,origin,global_spectra,spectra);
  );
}

void atlas__Trans__gathspec( const Trans* t, int nb_fields, int destination[], double spectra[], double global_spectra[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->gathspec(nb_fields,destination,spectra,global_spectra);
  );
}

void atlas__Trans__distgrid( const Trans* t, int nb_fields, int origin[], double global_fields[], double fields[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->distgrid(nb_fields,origin,global_fields,fields);
  );
}

void atlas__Trans__gathgrid( const Trans* t, int nb_fields, int destination[], double fields[], double global_fields[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->gathgrid(nb_fields,destination,fields,global_fields);
  );
}

void atlas__Trans__invtrans_scalar( const Trans* t, int nb_fields, double scalar_spectra[], double scalar_fields[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->invtrans(nb_fields,scalar_spectra,scalar_fields);
  );
}

void atlas__Trans__invtrans_vordiv2wind( const Trans* t, int nb_fields, double vorticity_spectra[], double divergence_spectra[], double wind_fields[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->invtrans(nb_fields,vorticity_spectra,divergence_spectra,wind_fields);
  );
}

void atlas__Trans__dirtrans_scalar( const Trans* t, int nb_fields, double scalar_fields[], double scalar_spectra[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->dirtrans(nb_fields,scalar_fields,scalar_spectra);
  );
}

void atlas__Trans__dirtrans_wind2vordiv( const Trans* t, int nb_fields, double wind_fields[], double vorticity_spectra[], double divergence_spectra[] )
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->dirtrans(nb_fields,wind_fields,vorticity_spectra,divergence_spectra);
  );
}

void atlas__Trans__specnorm (const Trans* t, int nb_fields, double spectra[], double norms[], int rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( t );
    return t->specnorm(nb_fields, spectra, norms, rank);
  );
}

int atlas__Trans__nspec2 (const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nb_spectral_coefficients();
  );
  return 0;
}

int atlas__Trans__nspec2g (const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nb_spectral_coefficients_global();
  );
  return 0;
}

int atlas__Trans__ngptot (const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nb_gridpoints();
  );
  return 0;
}

int atlas__Trans__ngptotg (const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nb_gridpoints_global();
  );
  return 0;
}

int atlas__Trans__truncation (const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->truncation();
  );
  return 0;
}

const Grid::Implementation* atlas__Trans__grid(const Trans* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( This->grid() );
    return This->grid().get();
  );
  return nullptr;
}


void atlas__Trans__dirtrans_fieldset_nodes (const Trans* This, const functionspace::detail::NodeColumns* gp, const field::FieldSetImpl* gpfields, const functionspace::detail::Spectral* sp, field::FieldSetImpl* spfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( gp );
    ASSERT( gpfields );
    ASSERT( sp );
    ASSERT( spfields );
    ASSERT( parameters );
    FieldSet fspfields(spfields);
    This->dirtrans(FunctionSpace(gp),gpfields,FunctionSpace(sp),fspfields,*parameters);
  );
}

void atlas__Trans__dirtrans_fieldset (const Trans* This, const field::FieldSetImpl* gpfields, field::FieldSetImpl* spfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( gpfields );
    ASSERT( spfields );
    ASSERT( parameters );
    FieldSet fspfields(spfields);
    This->dirtrans(gpfields,fspfields,*parameters);
  );
}

void atlas__Trans__invtrans_fieldset_nodes (const Trans* This, const functionspace::detail::Spectral* sp, const field::FieldSetImpl* spfields, const functionspace::detail::NodeColumns* gp, field::FieldSetImpl* gpfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( sp );
    ASSERT( spfields );
    ASSERT( gp );
    ASSERT( gpfields );
    ASSERT( parameters );
    FieldSet fgpfields(gpfields);
    This->invtrans(FunctionSpace(sp),spfields,FunctionSpace(gp),fgpfields,*parameters);
  );
}


void atlas__Trans__dirtrans_field (const Trans* This, const field::FieldImpl* gpfield, field::FieldImpl* spfield, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( spfield );
    ASSERT( gpfield );
    ASSERT( parameters );
    Field fspfield(spfield);
    This->dirtrans(gpfield,fspfield,*parameters);
  );
}

void atlas__Trans__dirtrans_field_nodes (const Trans* This, const functionspace::detail::NodeColumns* gp, const field::FieldImpl* gpfield, const functionspace::detail::Spectral* sp, field::FieldImpl* spfield, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( spfield );
    ASSERT( gpfield );
    ASSERT( parameters );
    Field fspfield(spfield);
    This->dirtrans(FunctionSpace(gp),gpfield,FunctionSpace(sp),fspfield,*parameters);
  );
}

void atlas__Trans__invtrans_fieldset (const Trans* This, const field::FieldSetImpl* spfields, field::FieldSetImpl* gpfields, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( spfields );
    ASSERT( gpfields );
    ASSERT( parameters );
    FieldSet fgpfields(gpfields);
    This->invtrans(spfields,fgpfields,*parameters);
  );
}

void atlas__Trans__invtrans_field (const Trans* This, const field::FieldImpl* spfield, field::FieldImpl* gpfield, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( spfield );
    ASSERT( gpfield );
    ASSERT( parameters );
    Field fgpfield(gpfield);
    This->invtrans(spfield,fgpfield,*parameters);
  );
}

void atlas__Trans__invtrans_field_nodes (const Trans* This, const functionspace::detail::Spectral* sp, const field::FieldImpl* spfield, const functionspace::detail::NodeColumns* gp, field::FieldImpl* gpfield, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( sp );
    ASSERT( spfield );
    ASSERT( gp );
    ASSERT( gpfield );
    ASSERT( parameters );
    Field fgpfield(gpfield);
    This->invtrans(FunctionSpace(sp),spfield,FunctionSpace(gp),fgpfield,*parameters);
  );
}

void atlas__Trans__dirtrans_wind2vordiv_field_nodes (const Trans* This, const functionspace::detail::NodeColumns* gp, const field::FieldImpl* gpwind, const functionspace::detail::Spectral* sp, field::FieldImpl* spvor, field::FieldImpl* spdiv, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( gp );
    ASSERT( gpwind );
    ASSERT( sp );
    ASSERT( spvor );
    ASSERT( spdiv );
    ASSERT( parameters );
    Field fspvor(spvor);
    Field fspdiv(spdiv);
    This->dirtrans_wind2vordiv(FunctionSpace(gp),gpwind,FunctionSpace(sp),fspvor,fspdiv,*parameters);
  );
}

void atlas__Trans__invtrans_vordiv2wind_field_nodes (const Trans* This, const functionspace::detail::Spectral* sp, const field::FieldImpl* spvor, const field::FieldImpl* spdiv, const functionspace::detail::NodeColumns* gp, field::FieldImpl* gpwind, const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( sp );
    ASSERT( spvor );
    ASSERT( spdiv );
    ASSERT( gp );
    ASSERT( gpwind );
    ASSERT( parameters );
    Field fgpwind(gpwind);
    This->invtrans_vordiv2wind(FunctionSpace(sp),spvor,spdiv,FunctionSpace(gp),fgpwind,*parameters);
  );
}

void atlas__Trans__invtrans (const Trans* This, int nb_scalar_fields, double scalar_spectra[], int nb_vordiv_fields, double vorticity_spectra[], double divergence_spectra[], double gp_fields[], const TransParameters* parameters)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
      This->invtrans( nb_scalar_fields, scalar_spectra,
                      nb_vordiv_fields, vorticity_spectra, divergence_spectra,
                      gp_fields,
                      *parameters );
  );
}

void atlas__Trans__invtrans_grad_field_nodes (const Trans* This, const functionspace::detail::Spectral* sp, const field::FieldImpl* spfield, const functionspace::detail::NodeColumns* gp, field::FieldImpl* gpfield)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    ASSERT( sp );
    ASSERT( spfield );
    ASSERT( gp );
    ASSERT( gpfield );
    Field fgpfield(gpfield);
    This->invtrans_grad(FunctionSpace(sp),spfield,FunctionSpace(gp),fgpfield);
  );
}



TransParameters* atlas__TransParameters__new ()
{
  return new TransParameters();
}

void atlas__TransParameters__delete (TransParameters* This)
{
  ATLAS_ERROR_HANDLING( ASSERT( This ); );
  delete This;
}

} // namespace trans
} // namespace atlas
