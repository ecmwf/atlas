/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE atlas_test_trans
#include "ecbuild/boost_test_framework.h"

#include <algorithm>    // std::min_element, std::max_element, std::copy

#include "eckit/mpi/mpi.h"
#include "eckit/config/ResourceMgr.h"
#include "atlas/atlas.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/TransPartitioner.h"
#include "atlas/grids/grids.h"
#include "atlas/meshgen/EqualAreaPartitioner.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/LogFormat.h"
#include "atlas/GridDistribution.h"
#include "atlas/io/Gmsh.h"

#include "transi/trans.h"

using namespace eckit;
using namespace atlas;
//using namespace atlas::trans;
using namespace atlas::grids;

struct Fixture   {
       Fixture() {

         atlas_init(boost::unit_test::framework::master_test_suite().argc,
                    boost::unit_test::framework::master_test_suite().argv);
         trans_init();

       }
      ~Fixture() {
         trans_finalize();
         atlas_finalize();
       }
};


void read_rspecg(trans::Trans& trans, std::vector<double>& rspecg, std::vector<int>& nfrom, int &nfld )
{
  eckit::Log::info() << "read_rspecg ...\n";
  nfld = 2;
  if( trans.myproc() == 1 )
  {
    rspecg.resize(nfld*trans.nspec2g());
    for( int i=0; i<trans.nspec2g(); ++i )
    {
      rspecg[i*nfld + 0] = (i==0 ? 1. : 0.); // scalar field 1
      rspecg[i*nfld + 1] = (i==0 ? 2. : 0.); // scalar field 2
    }
  }
  nfrom.resize(nfld);
  for (int jfld=0; jfld<nfld; ++jfld)
    nfrom[jfld] = 1;

  eckit::Log::info() << "read_rspecg ... done" << std::endl;
}

//void write_rgpg(Trans* trans, double* rgpg[], int nfld )
//{
//  int jfld;
//  if( trans->myproc == 1 ) printf("write_rgpg ...\n");
//  for( jfld=0; jfld<nfld; ++jfld )
//  {
//    // output global field rgpg[jfld]
//  }
//  if( trans->myproc == 1 ) printf("write_rgpg ... done\n");
//}

namespace atlas {
namespace trans {

class SpectralGlobal;
class SpectralDistributed;
class SpectralScalar;
class GridPointScalar;
class GridPointDistributed;
class GridPointGlobal;

void distspec( trans::Trans&, const SpectralGlobal&, const SpectralDistributed& );
void invtrans( trans::Trans&, const SpectralScalar&, const GridPointScalar& );
void gathgrid( trans::Trans&, const GridPointDistributed&, const GridPointGlobal& );


class SpectralGlobal
{
  friend void distspec( trans::Trans& , const SpectralGlobal& , const SpectralDistributed& );
private:
  virtual int nfld() const = 0;
  virtual double* rspecg() const = 0;
  virtual int* procs() const = 0;
  virtual void transfer() const = 0;
};

class SpectralDistributed
{
  friend void distspec( trans::Trans& , const SpectralGlobal& , const SpectralDistributed& );

private:
  virtual int nfld() const = 0;
  virtual double* rspec() const = 0;
  virtual void transfer() const = 0;
};

class SpectralScalar
{
  friend void invtrans( trans::Trans&, const SpectralScalar&, const GridPointScalar& );
private:
  virtual double* rspscalar() const = 0;
  virtual void transfer() const = 0;
};

class SpectralVorDiv
{
private:
  virtual double* rspvor() const = 0;
  virtual double* rspdiv() const = 0;
  virtual void transfer() const = 0;
};

class GridPointGlobal
{
  friend void gathgrid( trans::Trans&, const GridPointDistributed&, const GridPointGlobal& );
private:
  virtual int     nfld() const = 0;
  virtual double* rgpg() const = 0;
  virtual int*    procs() const = 0;
  virtual void transfer() const = 0;
};

class GridPointDistributed
{
  friend void gathgrid( trans::Trans&, const GridPointDistributed&, const GridPointGlobal& );
private:
  virtual int     nfld() const = 0;
  virtual double* rgp() const = 0;
  virtual void transfer() const = 0;
  virtual int nproma() const = 0;
  virtual int ngpblks() const = 0;
};

class GridPointScalar
{
  friend void invtrans( trans::Trans&, const SpectralScalar&, const GridPointScalar& );
private:
  virtual double* rgp() const = 0;
  virtual void transfer() const = 0;
  virtual int nproma() const = 0;
  virtual int ngpblks() const = 0;
};

class GridPointWind
{
private:
  virtual double* rgp() const = 0;
  virtual void transfer() const = 0;
};



// --------------------------------------------------------------------

class SpectralGlobalData: public SpectralGlobal
{
public:
  SpectralGlobalData( int _nfld, double* _rspecg, int* _procs )
  {
    nfld_ = _nfld;
    rspecg_ = _rspecg;
    procs_ = _procs;
    procs_owned_ = false;
  }
  SpectralGlobalData( int _nfld, double* _rspecg )
  {
    nfld_ = _nfld;
    rspecg_ = _rspecg;
    procs_ = new int(nfld_);
    for( int j=0; j<nfld_; ++j )
      procs_[j] = 1;
    procs_owned_ = true;
  }
  SpectralGlobalData( int _nfld, const std::vector<double>& _rspecg, const std::vector<int>& _procs )
  {
    ASSERT(_nfld == _procs.size());
    nfld_ = _nfld;
    rspecg_ = const_cast<std::vector<double>&>(_rspecg).data();
    procs_ = const_cast<std::vector<int>&>(_procs).data();
    procs_owned_ = false;
  }
  SpectralGlobalData( int nfld, const std::vector<double>& _rspecg )
  {
    nfld_ = nfld;
    rspecg_ = const_cast<std::vector<double>&>(_rspecg).data();
    procs_ = new int(nfld_);
    for( int j=0; j<nfld_; ++j )
      procs_[j] = 1;
    procs_owned_ = true;
  }
  virtual ~SpectralGlobalData()
  {
    if( procs_owned_ )
      delete[] procs_;
  }

private: // methods
  virtual int nfld() const { return nfld_; }
  virtual double* rspecg() const { return rspecg_; }
  virtual int* procs() const { return procs_; }
  virtual void transfer() const {} // nothing to do

private: // data
  int nfld_;
  double* rspecg_;
  int* procs_;
  bool procs_owned_;
};

// --------------------------------------------------------------------

class SpectralDistributedData: public SpectralDistributed
{
public:
  SpectralDistributedData( int _nfld, double* _rspec )
  {
    nfld_ = _nfld;
    rspec_ = _rspec;
  }
  SpectralDistributedData( int _nfld, const std::vector<double>& _rspec )
  {
    nfld_ = _nfld;
    rspec_ = const_cast<std::vector<double>&>(_rspec).data();
  }


private: // methods
  virtual int nfld() const { return nfld_; }
  virtual double* rspec() const { return rspec_; }
  virtual void transfer() const {} // nothing to do

private: // data
  int nfld_;
  double* rspec_;
};

// --------------------------------------------------------------------


class SpectralScalarData: public SpectralScalar
{
public:
  SpectralScalarData( double* _rspscalar )
  {
    rspscalar_ = _rspscalar;
  }
  SpectralScalarData( const std::vector<double>& _rspscalar )
  {
    rspscalar_ = const_cast<std::vector<double>&>(_rspscalar).data();
  }


private: // methods
  virtual double* rspscalar() const { return rspscalar_; }
  virtual void transfer() const {} // nothing to do

private: // data
  double* rspscalar_;
};


// --------------------------------------------------------------------


class GridPointScalarData: public GridPointScalar
{
public:
  GridPointScalarData( double* _rgp )
  {
    rgp_ = _rgp;
  }
  GridPointScalarData( const std::vector<double>& _rgp )
  {
    rgp_ = const_cast<std::vector<double>&>(_rgp).data();
  }


private: // methods
  virtual double* rgp() const { return rgp_; }
  virtual void transfer() const {} // nothing to do
  virtual int nproma() const  { return 0; }
  virtual int ngpblks() const { return 1; }

private: // data
  double* rgp_;
};


// --------------------------------------------------------------------

class GridPointGlobalData: public GridPointGlobal
{
public:
  GridPointGlobalData( int _nfld, double* _rgpg, int* _procs )
  {
    nfld_ = _nfld;
    rgpg_ = _rgpg;
    procs_ = _procs;
    procs_owned_ = false;
  }
  GridPointGlobalData( int _nfld, double* _rgpg )
  {
    nfld_ = _nfld;
    rgpg_ = _rgpg;
    procs_ = new int(nfld_);
    for( int j=0; j<nfld_; ++j )
      procs_[j] = 1;
    procs_owned_ = true;
  }
  GridPointGlobalData( int _nfld, const std::vector<double>& _rgpg, const std::vector<int>& _procs )
  {
    ASSERT(_nfld == _procs.size());
    nfld_ = _nfld;
    rgpg_ = const_cast<std::vector<double>&>(_rgpg).data();
    procs_ = const_cast<std::vector<int>&>(_procs).data();
    procs_owned_ = false;
  }
  GridPointGlobalData( int nfld, const std::vector<double>& _rgpg )
  {
    nfld_ = nfld;
    rgpg_ = const_cast<std::vector<double>&>(_rgpg).data();
    procs_ = new int(nfld_);
    for( int j=0; j<nfld_; ++j )
      procs_[j] = 1;
    procs_owned_ = true;
  }
  virtual ~GridPointGlobalData()
  {
    if( procs_owned_ )
      delete[] procs_;
  }

private: // methods
  virtual int nfld() const { return nfld_; }
  virtual double* rgpg() const { return rgpg_; }
  virtual int* procs() const { return procs_; }
  virtual void transfer() const {} // nothing to do

private: // data
  int nfld_;
  double* rgpg_;
  int* procs_;
  bool procs_owned_;
};

// --------------------------------------------------------------------

class GridPointDistributedData: public GridPointDistributed
{
public:
  GridPointDistributedData( int _nfld, double* _rgp )
  {
    nfld_ = _nfld;
    rgp_ = _rgp;
  }
  GridPointDistributedData( int _nfld, const std::vector<double>& _rgp )
  {
    nfld_ = _nfld;
    rgp_ = const_cast<std::vector<double>&>(_rgp).data();
  }


private: // methods
  virtual int nfld() const { return nfld_; }
  virtual double* rgp() const { return rgp_; }
  virtual void transfer() const {} // nothing to do
  virtual int nproma() const  { return 0; }
  virtual int ngpblks() const { return 1; }

private: // data
  int nfld_;
  double* rgp_;
};

// --------------------------------------------------------------------

class GridPointDistributedField: public GridPointDistributed
{
public:
  GridPointDistributedField( Field& field )
  {
    field_ = &field;
    nfld_ = field_->nb_vars();
    ngptot_ = field_->size()/nfld_;
  }

private: // methods
  virtual int nfld() const { return nfld_; }
  virtual double* rgp() const { Log::warning() << "ghost nodes are inside!!!" << std::endl; return field_->data<double>(); }
  virtual void transfer() const {}
  virtual int nproma() const  { return 1; }
  virtual int ngpblks() const { return ngptot_; }

private: // data
  int nfld_;
  int ngptot_;
  mutable std::vector<double> rgp_;
  Field* field_;
};

// --------------------------------------------------------------------

class GridPointScalarField: public GridPointScalar
{
public:
  GridPointScalarField( Field& field )
  {
    field_ = &field;
    nfld_ = field_->nb_vars();
    ngptot_ = field_->size()/nfld_;
  }

private: // methods
  virtual double* rgp() const { Log::warning() << "ghost nodes are inside!!!" << std::endl; return field_->data<double>();  }
  virtual void transfer() const {}
  virtual int nproma() const  { return 1; }
  virtual int ngpblks() const { return ngptot_; }

private: // data
  int nfld_;
  int ngptot_;
  mutable std::vector<double> rgp_;
  Field* field_;
};

// --------------------------------------------------------------------


void distspec( trans::Trans& t, const SpectralGlobal& glb, const SpectralDistributed& dist)
{
  ASSERT( glb.nfld() == dist.nfld() );
  struct ::DistSpec_t args = new_distspec(t);
    args.nfld = glb.nfld();
    args.rspecg = glb.rspecg();
    args.nfrom = glb.procs();
    args.rspec = dist.rspec();
  ::trans_distspec(&args);
  dist.transfer();
}

void invtrans( trans::Trans& t, const SpectralScalar& sp, const GridPointScalar& gp )
{
  struct ::InvTrans_t args = new_invtrans(t);
    args.rspscalar = sp.rspscalar();
    args.rgp = gp.rgp();
    args.ngpblks = gp.ngpblks() != 1 ? t.ngptot() : 1;
    args.nproma = gp.nproma() != 0 ? gp.nproma() : t.ngptot();
  ::trans_invtrans(&args);
  gp.transfer();
}

void gathgrid( trans::Trans& t, const GridPointDistributed& dist, const GridPointGlobal& glb )
{
  ASSERT( glb.nfld() == dist.nfld() );
  struct ::GathGrid_t args = new_gathgrid(t);
    args.nfld = glb.nfld();
    args.rgpg = glb.rgpg();
    args.nto  = glb.procs();
    args.rgp  = dist.rgp();
    args.ngpblks = dist.ngpblks() != 1 ? t.ngptot() : 1;
    args.nproma = dist.nproma() != 0 ? 1 : t.ngptot();
  ::trans_gathgrid(&args);
  glb.transfer();
}

//gathspec( trans, SpectralDistributed(rspec),    SpectralGlobal(rspecg,nto) );

//distgrid( trans, GridPointGlobal(rgpg,nfrom),   GridPointDistributed(rgp) );

//gathgrid( trans, GridPointDistributed(rgp),     GridPointGlobal(rgpg,nto) );

//invtrans( trans, SpectralScalar(rspscalar),     GridPointScalar(rgp) );

//invtrans( trans, SpectralVorDiv(rspvor,rspdiv), GridPointWind(rgp) );

//dirtrans( trans, GridPointScalar(rgp),          SpectralScalar(rspscalar) );

//dirtrans( trans, GridPointWind(rgp),            SpectralVorDiv(rspvor,rspdiv) );


/*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void invtrans( trans::Trans& t, const int nb_fields, const double vorticity_spectra[], const double divergence_spectra[], double wind_fields[] )
{
  struct ::InvTrans_t args = new_invtrans(t);
    args.nvordiv = nb_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp = wind_fields;
  if( int errcode = ::trans_invtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("invtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

/*!
 * @brief Direct transform of scalar fields
 */
void dirtrans( trans::Trans& t, const int nb_fields, const double scalar_fields[], double scalar_spectra[] )
{
  struct ::DirTrans_t args = new_dirtrans(t);
    args.nscalar = nb_fields;
    args.rgp = scalar_fields;
    args.rspscalar = scalar_spectra;
  if( int errcode = ::trans_dirtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("dirtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

/*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void dirtrans( trans::Trans& t, const int nb_fields, const double wind_fields[], double vorticity_spectra[], double divergence_spectra[] )
{
  struct ::DirTrans_t args = new_dirtrans(t);
    args.nvordiv = nb_fields;
    args.rspvor = vorticity_spectra;
    args.rspdiv = divergence_spectra;
    args.rgp    = wind_fields;
  if( int errcode = ::trans_dirtrans(&args) != TRANS_SUCCESS )
    throw eckit::Exception("dirtrans failed: "+std::string(::trans_error_msg(errcode)),Here());
}

class DistSpec
{
public:
  DistSpec(Trans& trans)
  {
    distspec_ = new_distspec(trans);
  }

  DistSpec& nfld(int _nfld) { distspec_.nfld = _nfld; return *this; }

  DistSpec& nfrom (const int*    _nfrom)  { distspec_.nfrom  = const_cast<int*>(_nfrom);     return *this; }
  DistSpec& rspecg(const double* _rspecg) { distspec_.rspecg = const_cast<double*>(_rspecg); return *this; }
  DistSpec& rspec (      double* _rspec)  { distspec_.rspec  = _rspec;                        return *this; }

  DistSpec& nfrom (const std::vector<int>&    _nfrom)  { return nfrom(_nfrom.size()?_nfrom.data():NULL); }
  DistSpec& rspecg(const std::vector<double>& _rspecg) { return rspecg(_rspecg.size()?_rspecg.data():NULL); }
  DistSpec& rspec (      std::vector<double>& _rspec)
  { ASSERT(_rspec.size()); return rspec(_rspec.data()); }

  void execute() { ::trans_distspec(&distspec_); }

private:
  struct DistSpec_t distspec_;
};
}
}

BOOST_GLOBAL_FIXTURE( Fixture )

BOOST_AUTO_TEST_CASE( test_trans_distribution_matches_atlas )
{
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg.N80" ) );

  BOOST_CHECK_EQUAL( g->nlat() , 160 );

  trans::Trans trans( *g );

  BOOST_CHECK_EQUAL( trans.nsmax() , 0 );

  trans::TransPartitioner partitioner(*g,trans);
  GridDistribution distribution( partitioner );

  // -------------- do checks -------------- //
  BOOST_CHECK_EQUAL( trans.nproc(),  eckit::mpi::size() );
  BOOST_CHECK_EQUAL( trans.myproc(), eckit::mpi::rank()+1 );

  if( eckit::mpi::rank() == 0 ) // all tasks do the same, so only one needs to check
  {
    int max_nb_regions_EW(0);
    for( int j=0; j<partitioner.nb_bands(); ++j )
      max_nb_regions_EW = std::max(max_nb_regions_EW, partitioner.nb_regions(j));

    BOOST_CHECK_EQUAL( trans.n_regions_NS(), partitioner.nb_bands() );
    BOOST_CHECK_EQUAL( trans.n_regions_EW(), max_nb_regions_EW );

    BOOST_CHECK_EQUAL( distribution.nb_partitions(), eckit::mpi::size() );
    BOOST_CHECK_EQUAL( distribution.partition().size(), g->npts() );

    std::vector<int> npts(distribution.nb_partitions(),0);

    for( int j=0; j<g->npts(); ++j )
      ++npts[distribution.partition(j)];

    BOOST_CHECK_EQUAL( trans.ngptotg(), g->npts() );
    BOOST_CHECK_EQUAL( trans.ngptot(),  npts[eckit::mpi::rank()] );
    BOOST_CHECK_EQUAL( trans.ngptotmx(), *std::max_element(npts.begin(),npts.end()) );

    for( int j=0; j<partitioner.nb_bands(); ++j )
      BOOST_CHECK_EQUAL( trans.n_regions()[j] , partitioner.nb_regions(j) );
  }

}

BOOST_AUTO_TEST_CASE( test_trans_partitioner )
{
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg.N80" ) );

  trans::Trans trans( *g, 0 );

  BOOST_CHECK_EQUAL( trans.nsmax() , 0 );
  BOOST_CHECK_EQUAL( trans.ngptotg() , g->npts() );

}

BOOST_AUTO_TEST_CASE( test_trans_options )
{
  trans::Trans::Options opts;
  opts.set_fft(trans::FFTW);
  opts.set_split_latitudes(false);
  opts.set_read("readfile");

  eckit::Log::info() << "trans_opts = " << opts << std::endl;
}


BOOST_AUTO_TEST_CASE( test_distspec )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "oct.N80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  meshgen::ReducedGridMeshGenerator generate;
  BOOST_CHECKPOINT("mesh generator created");
  //trans::Trans trans(*g, 159 );

  trans::Trans::Options p;
  if( eckit::mpi::size() == 1 )
    p.set_write("cached_legendre_coeffs");
  p.set_flt(false);
  trans::Trans trans(400, 159, p);
  BOOST_CHECKPOINT("Trans initialized");
  std::vector<double> rspecg;
  std::vector<int   > nfrom;
  int nfld;
  BOOST_CHECKPOINT("Read rspecg");
  read_rspecg(trans,rspecg,nfrom,nfld);


  std::vector<double> rspec(nfld*trans.nspec2());
  std::vector<int> nto(nfld,1);
  std::vector<double> rgp(nfld*trans.ngptot());
  std::vector<double> rgpg(nfld*trans.ngptotg());

  BOOST_CHECKPOINT("distspec");

  BOOST_CHECK_NO_THROW(
    trans::DistSpec(trans).
        nfld (nfld).
        nfrom (nfrom).
        rspecg (rspecg).
        rspec (rspec).
    execute();
  )

  trans.distspec( nfld, nfrom.data(), rspecg.data(), rspec.data() );
  trans.invtrans( nfld, rspec.data(), rgp.data() );
  trans.gathgrid( nfld, nto.data(),   rgp.data(),    rgpg.data() );

  trans::distspec(trans, trans::SpectralGlobalData(nfld,rspecg,nfrom), trans::SpectralDistributedData(nfld,rspec) );

  trans::invtrans(trans, trans::SpectralScalarData(rspec), trans::GridPointScalarData(rgp) );

  trans::gathgrid(trans, trans::GridPointDistributedData(nfld,rgp), trans::GridPointGlobalData(nfld,rgpg,nto) );

  atlas::Mesh* m = generate(*g);
  Field& f = m->function_space("nodes").create_field<double>("transform_me",nfld);

  trans::invtrans(trans, trans::SpectralScalarData(rspec), trans::GridPointScalarField(f) );

  //trans::gathgrid(trans, trans::GridPointDistributedField(f), trans::GridPointGlobalData(nfld,rgpg,nto) );

  BOOST_CHECKPOINT("end test_distspec");

}

BOOST_AUTO_TEST_CASE( test_distribution )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "oct.N80" ) );

  BOOST_CHECKPOINT("test_distribution");

  GridDistribution::Ptr d_trans( trans::TransPartitioner(*g).distribution() );
  BOOST_CHECKPOINT("trans distribution created");

  GridDistribution::Ptr d_eqreg( meshgen::EqualAreaPartitioner(*g).distribution() );

  BOOST_CHECKPOINT("eqregions distribution created");

  if( eckit::mpi::rank() == 0 )
  {
    BOOST_CHECK_EQUAL( d_trans->nb_partitions(), d_eqreg->nb_partitions() );
    BOOST_CHECK_EQUAL( d_trans->max_pts(), d_eqreg->max_pts() );
    BOOST_CHECK_EQUAL( d_trans->min_pts(), d_eqreg->min_pts() );

    BOOST_CHECK_EQUAL_COLLECTIONS( d_trans->nb_pts().begin(), d_trans->nb_pts().end(),
                                   d_eqreg->nb_pts().begin(), d_eqreg->nb_pts().end() );
  }

}


BOOST_AUTO_TEST_CASE( test_generate_mesh )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "oct.N80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  eckit::ResourceMgr::instance().set("atlas.meshgen.triangulate","true");

  meshgen::ReducedGridMeshGenerator generate;
  trans::Trans trans(*g);

  Mesh::Ptr m_default( generate( *g ) );
  Mesh::Ptr m_trans( generate( *g, trans::TransPartitioner(*g).distribution() ) );
  Mesh::Ptr m_eqreg( generate( *g, meshgen::EqualAreaPartitioner(*g).distribution() ) );

  ArrayView<int,1> p_default( m_default->function_space("nodes").field("partition") );
  ArrayView<int,1> p_trans  ( m_trans  ->function_space("nodes").field("partition") );
  ArrayView<int,1> p_eqreg  ( m_eqreg  ->function_space("nodes").field("partition") );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_trans  .begin(), p_trans  .end() );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_eqreg  .begin(), p_eqreg  .end() );

  //Mesh::Ptr mesh ( generate(*g, meshgen::EqualAreaPartitioner(*g).distribution() ) );

  io::Gmsh().write(*m_trans,"N16_trans.msh");
}

