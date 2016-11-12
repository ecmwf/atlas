/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestFunctionSpace
#include "ecbuild/boost_test_framework.h"

#include "eckit/types/Types.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/atlas.h"
#include "atlas/internals/Debug.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/grid/Grid.h"
#include "atlas/field/Field.h"
#include "atlas/grid/gaussian/ReducedGaussian.h"
#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif
using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

struct AtlasFixture {
    AtlasFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_functionspace_NodeColumns )
{
  //ScopedPtr<grid::Grid> grid( Grid::create("O2") );

  size_t nlat = 2;
  long nlon[] = {4,8};
  ScopedPtr<grid::Grid> grid( new grid::gaussian::ReducedGaussian( nlat, nlon ) );

  mesh::Mesh mesh;
  mesh::generators::Structured generator;
  //generator.options.set("3d",true);
  generator.generate(*grid,mesh);

  //grid.reset();

  DEBUG();
  SharedPtr<functionspace::NodeColumns> nodes_fs( new functionspace::NodeColumns(mesh,mesh::Halo(1)) );
  DEBUG();
  size_t nb_levels = 10;
  //NodesColumnFunctionSpace columns_fs("columns",mesh,nb_levels,Halo(1));

  //BOOST_CHECK_EQUAL( nodes_fs.nb_nodes() , columns_fs.nb_nodes() );
  //BOOST_CHECK_EQUAL( columns_fs.nb_levels() , 10 );


  SharedPtr<field::Field> surface_scalar_field( nodes_fs->createField<double>("scalar") );
  SharedPtr<field::Field> surface_vector_field( nodes_fs->createField<double>("vector",array::make_shape(2)) );
  SharedPtr<field::Field> surface_tensor_field( nodes_fs->createField<double>("tensor",array::make_shape(2,2)) );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( surface_vector_field->name() , std::string("vector") );
  BOOST_CHECK_EQUAL( surface_tensor_field->name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector_field->size() , nodes_fs->nb_nodes()*2 );
  BOOST_CHECK_EQUAL( surface_tensor_field->size() , nodes_fs->nb_nodes()*2*2 );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );
  BOOST_CHECK_EQUAL( surface_vector_field->rank() , 2 );
  BOOST_CHECK_EQUAL( surface_tensor_field->rank() , 3 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( *surface_scalar_field );
  array::ArrayView<double,2> surface_vector = array::make_view<double,2>( *surface_vector_field );
  array::ArrayView<double,3> surface_tensor = array::make_view<double,3>( *surface_tensor_field );

  BOOST_CHECK_EQUAL( surface_scalar.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( surface_tensor.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector.shape(1), 2 );
  BOOST_CHECK_EQUAL( surface_tensor.shape(1), 2 );
  BOOST_CHECK_EQUAL( surface_tensor.shape(2), 2 );

  SharedPtr<field::Field> columns_scalar_field( nodes_fs->createField<double>("scalar",nb_levels) );
  SharedPtr<field::Field> columns_vector_field( nodes_fs->createField<double>("vector",nb_levels,array::make_shape(2)) );
  SharedPtr<field::Field> columns_tensor_field( nodes_fs->createField<double>("tensor",nb_levels,array::make_shape(2,2)) );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( columns_vector_field->name() , std::string("vector") );
  BOOST_CHECK_EQUAL( columns_tensor_field->name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nodes_fs->nb_nodes()*nb_levels );
  BOOST_CHECK_EQUAL( columns_vector_field->size() , nodes_fs->nb_nodes()*nb_levels*2 );
  BOOST_CHECK_EQUAL( columns_tensor_field->size() , nodes_fs->nb_nodes()*nb_levels*2*2 );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );
  BOOST_CHECK_EQUAL( columns_vector_field->rank() , 3 );
  BOOST_CHECK_EQUAL( columns_tensor_field->rank() , 4 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( *columns_scalar_field );
  array::ArrayView<double,3> columns_vector = array::make_view<double,3>( *columns_vector_field );
  array::ArrayView<double,4> columns_tensor = array::make_view<double,4>( *columns_tensor_field );

  BOOST_CHECK_EQUAL( columns_scalar.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( columns_vector.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( columns_tensor.shape(0), nodes_fs->nb_nodes() );
  BOOST_CHECK_EQUAL( columns_scalar.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_vector.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_tensor.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_vector.shape(2), 2 );
  BOOST_CHECK_EQUAL( columns_tensor.shape(2), 2 );
  BOOST_CHECK_EQUAL( columns_tensor.shape(3), 2 );

  field::Field::Ptr field( nodes_fs->createField<int>("partition",nb_levels) );
  array::ArrayView<int,2> arr = array::make_view<int,2>(*field);
  arr.assign(eckit::mpi::rank());
  //field->dump( Log::info() );
  nodes_fs->haloExchange(*field);
  //field->dump( Log::info() );

  field::Field::Ptr field2( nodes_fs->createField<int>("partition2",nb_levels,array::make_shape(2)) );
  Log::info() << "field2.rank() = " << field2->rank() << std::endl;
  array::ArrayView<int,3> arr2 = array::make_view<int,3>(*field2);

  arr2.assign(eckit::mpi::rank());
  
  BOOST_CHECKPOINT(__LINE__);
  
  //field2->dump( Log::info() );
  nodes_fs->haloExchange(*field2);
  //field2->dump( Log::info() );

  BOOST_CHECKPOINT(__LINE__);

  Log::info() << nodes_fs->checksum(*field) << std::endl;

  BOOST_CHECKPOINT(__LINE__);

  size_t root = eckit::mpi::size()-1;
  field::Field::Ptr glb_field( nodes_fs->createField("partition",*field,field::global(root)) );
  nodes_fs->gather(*field,*glb_field);

  Log::info() << "local points = " << nodes_fs->nb_nodes() << std::endl;
  Log::info() << "grid points = " << grid->npts() << std::endl;
  Log::info() << "glb_field.shape(0) = " << glb_field->shape(0) << std::endl;

  BOOST_CHECK_EQUAL( glb_field->metadata().get<bool>("global"), true );
  BOOST_CHECK_EQUAL( glb_field->metadata().get<int>("owner"),  root );

  //glb_field->dump( Log::info() );

  if( eckit::mpi::rank() == root )
    glb_field->metadata().set("test_broadcast",123);

  arr.assign(-1);
  nodes_fs->scatter(*glb_field,*field);
  BOOST_CHECK_EQUAL( field->metadata().get<int>("test_broadcast"), 123 );
  nodes_fs->haloExchange(*field);
  //field->dump( Log::info() );

  Log::info() << nodes_fs->checksum(*field) << std::endl;

  field::FieldSet fields;
  fields.add(*field);
  fields.add(*field2);
  Log::info() << nodes_fs->checksum(fields) << std::endl;



  Log::info() << "Testing collectives for nodes scalar field" << std::endl;
  BOOST_TEST_CHECKPOINT("Testing collectives for nodes scalar field");
  {
    const field::Field& field = *surface_scalar_field;
    const functionspace::NodeColumns& fs = *nodes_fs;

  double max;
  double min;
  double sum;
  double mean;
  double stddev;
  size_t N;
  gidx_t gidx_max;
  gidx_t gidx_min;

  array::ArrayView<double,1> sfc_arr = array::make_view<double,1>( field );
  sfc_arr.assign( eckit::mpi::rank()+1 );
  fs.maximum(*surface_scalar_field,max);
  BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );

  fs.minimum(*surface_scalar_field,min);
  BOOST_CHECK_EQUAL( min, 1 );

  fs.maximumAndLocation(field,max,gidx_max);
  BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );
  Log::info() << "global index for maximum: " << gidx_max << std::endl;

  fs.minimumAndLocation(field,min,gidx_min);
  BOOST_CHECK_EQUAL( min, 1 );
  Log::info() << "global index for minimum: " << gidx_min << std::endl;

  fs.orderIndependentSum(field,sum,N);
  Log::info() << "oisum: " << sum << std::endl;
  Log::info() << "oiN: " << N << std::endl;

  fs.sum(field,sum,N);
  Log::info() << "sum: " << sum << std::endl;
  Log::info() << "N: " << N << std::endl;

  fs.mean(field,mean,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "N: " << N << std::endl;

  fs.meanAndStandardDeviation(field,mean,stddev,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "standard deviation: " << stddev << std::endl;
  Log::info() << "N: " << N << std::endl;

  int sumint;
  fs.orderIndependentSum(field,sumint,N);
  Log::info() << "sumint: " << sumint << std::endl;

  fs.sum(field,sumint,N);
  Log::info() << "sumint: " << sumint << std::endl;

  }


  Log::info() << "Testing collectives for nodes vector field" << std::endl;
  BOOST_TEST_CHECKPOINT("Testing collectives for nodes vector field");
  {
    const field::Field& field = *surface_vector_field;
    const functionspace::NodeColumns& fs = *nodes_fs;

    std::vector<double> max;
    std::vector<double> min;
    std::vector<double> sum;
    std::vector<double> mean;
    std::vector<double> stddev;
    size_t N;
    std::vector<gidx_t> gidx_max;
    std::vector<gidx_t> gidx_min;

    array::ArrayView<double,2> vec_arr = array::make_view<double,2>( field );
    vec_arr.assign( eckit::mpi::rank()+1 );
    fs.maximum(field,max);
    std::vector<double> check_max(field.stride(0),eckit::mpi::size());
    BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );

    fs.minimum(field,min);
    std::vector<double> check_min(field.stride(0),1);
    BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );

    fs.maximumAndLocation(field,max,gidx_max);
    BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );
    Log::info() << "global index for maximum: " << gidx_max << std::endl;

    fs.minimumAndLocation(field,min,gidx_min);
    BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );
    Log::info() << "global index for minimum: " << gidx_min << std::endl;

    fs.orderIndependentSum(field,sum,N);
    Log::info() << "oisum: " << sum << std::endl;
    Log::info() << "oiN: " << N << std::endl;

    fs.mean(field,mean,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.meanAndStandardDeviation(field,mean,stddev,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "standard deviation: " << stddev << std::endl;
    Log::info() << "N: " << N << std::endl;

    std::vector<int> sumint;
    fs.orderIndependentSum(field,sumint,N);
    Log::info() << "sumint: " << sumint << std::endl;

    fs.sum(field,sumint,N);
    Log::info() << "sumint: " << sumint << std::endl;

  }

  Log::info() << "Testing collectives for columns scalar field" << std::endl;
  BOOST_TEST_CHECKPOINT("Testing collectives for columns scalar field");
  if(1){
    const field::Field& field = *columns_scalar_field;
    const functionspace::NodeColumns& fs = *nodes_fs;
    double max;
    double min;
    double sum;
    double mean;
    double stddev;
    size_t N;
    gidx_t gidx_max;
    gidx_t gidx_min;
    size_t level;

    BOOST_CHECK_EQUAL(field.has_levels(),true);

    array::ArrayView<double,2> arr = array::make_view<double,2>( field );
    arr.assign( eckit::mpi::rank()+1 );
    fs.maximum(field,max);
    BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );

    fs.minimum(field,min);
    BOOST_CHECK_EQUAL( min, 1 );

    fs.maximumAndLocation(field,max,gidx_max,level);
    BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );
    Log::info() << "global index for maximum: " << gidx_max << std::endl;
    Log::info() << "level for maximum: " << level << std::endl;

    fs.minimumAndLocation(field,min,gidx_min,level);
    BOOST_CHECK_EQUAL( min, 1 );
    Log::info() << "global index for minimum: " << gidx_min << std::endl;
    Log::info() << "level for minimum: " << level << std::endl;

    fs.orderIndependentSum(field,sum,N);
    Log::info() << "order independent sum: " << sum << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.sum(field,sum,N);
    Log::info() << "sum: " << sum << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.mean(field,mean,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.meanAndStandardDeviation(field,mean,stddev,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "standard deviation: " << stddev << std::endl;
    Log::info() << "N: " << N << std::endl;

    int sumint;
    fs.orderIndependentSum(field,sumint,N);
    Log::info() << "order independent sum in int: " << sumint << std::endl;

    fs.sum(field,sumint,N);
    Log::info() << "sum in int: " << sumint << std::endl;

    field::Field::Ptr max_per_level    ( field::Field::create<double>("max",    array::make_shape(nb_levels)) );
    field::Field::Ptr min_per_level    ( field::Field::create<double>("min",    array::make_shape(nb_levels)) );
    field::Field::Ptr sum_per_level    ( field::Field::create<double>("sum",    array::make_shape(nb_levels)) );
    field::Field::Ptr mean_per_level   ( field::Field::create<double>("mean",   array::make_shape(nb_levels)) );
    field::Field::Ptr stddev_per_level ( field::Field::create<double>("stddev", array::make_shape(nb_levels)) );
    field::Field::Ptr gidx_per_level   ( field::Field::create<gidx_t>("gidx",   array::make_shape(nb_levels)) );

    fs.maximumPerLevel(field,*max_per_level);
    max_per_level->dump(Log::info());
    fs.minimumPerLevel(field,*min_per_level);
    min_per_level->dump(Log::info());
    fs.sumPerLevel(field,*sum_per_level,N);
    sum_per_level->dump(Log::info());
    fs.meanPerLevel(field,*mean_per_level,N);
    mean_per_level->dump(Log::info());
    fs.meanAndStandardDeviationPerLevel(field,*mean_per_level,*stddev_per_level,N);
    mean_per_level->dump(Log::info());
    stddev_per_level->dump(Log::info());
    fs.orderIndependentSumPerLevel(field,*sum_per_level,N);
    sum_per_level->dump(Log::info());

  }

  BOOST_TEST_CHECKPOINT("Testing collectives for columns vector field");
  if(1){
    const field::Field& field = *columns_vector_field;
    const functionspace::NodeColumns& fs = *nodes_fs;
    size_t nvar = field.stride(1);
    std::vector<double> max;
    std::vector<double> min;
    std::vector<double> sum;
    std::vector<double> mean;
    std::vector<double> stddev;
    size_t N;
    std::vector<gidx_t> gidx_max;
    std::vector<gidx_t> gidx_min;
    std::vector<size_t> levels;

    array::ArrayView<double,3> vec_arr = array::make_view<double,3>( field );
    vec_arr.assign( eckit::mpi::rank()+1 );
    fs.maximum(field,max);
    std::vector<double> check_max(nvar,eckit::mpi::size());
    BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );

    fs.minimum(field,min);
    std::vector<double> check_min(nvar,1);
    BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );

    fs.maximumAndLocation(field,max,gidx_max,levels);
    BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );
    Log::info() << "global index for maximum: " << gidx_max << std::endl;

    fs.minimumAndLocation(field,min,gidx_min,levels);
    BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );
    Log::info() << "global index for minimum: " << gidx_min << std::endl;

    fs.orderIndependentSum(field,sum,N);
    Log::info() << "sum: " << sum << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.mean(field,mean,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "N: " << N << std::endl;

    fs.meanAndStandardDeviation(field,mean,stddev,N);
    Log::info() << "mean: " << mean << std::endl;
    Log::info() << "standard deviation: " << stddev << std::endl;
    Log::info() << "N: " << N << std::endl;

    std::vector<int> sumint;
    fs.orderIndependentSum(field,sumint,N);
    Log::info() << "sumint: " << sumint << std::endl;

    fs.sum(field,sumint,N);
    Log::info() << "sumint: " << sumint << std::endl;

    field::Field::Ptr max_per_level    ( field::Field::create<double>("max",    array::make_shape(nb_levels,nvar)) );
    field::Field::Ptr min_per_level    ( field::Field::create<double>("min",    array::make_shape(nb_levels,nvar)) );
    field::Field::Ptr sum_per_level    ( field::Field::create<double>("sum",    array::make_shape(nb_levels,nvar)) );
    field::Field::Ptr mean_per_level   ( field::Field::create<double>("mean",   array::make_shape(nb_levels,nvar)) );
    field::Field::Ptr stddev_per_level ( field::Field::create<double>("stddev", array::make_shape(nb_levels,nvar)) );
    field::Field::Ptr gidx_per_level   ( field::Field::create<gidx_t>("gidx",   array::make_shape(nb_levels,nvar)) );

    fs.maximumPerLevel(field,*max_per_level);
    max_per_level->dump(Log::info());

    fs.minimumPerLevel(field,*min_per_level);
    min_per_level->dump(Log::info());

    fs.sumPerLevel(field,*sum_per_level,N);
    sum_per_level->dump(Log::info());

    fs.meanPerLevel(field,*mean_per_level,N);
    mean_per_level->dump(Log::info());

    fs.meanAndStandardDeviationPerLevel(field,*mean_per_level,*stddev_per_level,N);
    mean_per_level->dump(Log::info());
    stddev_per_level->dump(Log::info());

    fs.orderIndependentSumPerLevel(field,*sum_per_level,N);
    sum_per_level->dump(Log::info());
  }



  SharedPtr<field::Field> tmp (
        nodes_fs->createField( field::datatypeT<double>() | field::global(0) | field::levels(10) | field::name("tmp") ) );

}

BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace )
{
  size_t truncation = 159;
  size_t nb_levels = 10;
  size_t nspec2g = (truncation+1)*(truncation+2);

  SharedPtr<Spectral> spectral_fs( new Spectral(truncation) );

  SharedPtr<field::Field> surface_scalar_field( spectral_fs->createField<double>("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( *surface_scalar_field );

  BOOST_CHECK_EQUAL( surface_scalar.shape(0), nspec2g );

  SharedPtr<field::Field> columns_scalar_field( spectral_fs->createField<double>("scalar",nb_levels) );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2g*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( *columns_scalar_field );

  BOOST_CHECK_EQUAL( columns_scalar.shape(0), nspec2g );
  BOOST_CHECK_EQUAL( columns_scalar.shape(1), nb_levels );

}

#ifdef ATLAS_HAVE_TRANS

BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_dist )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2 = trans.nspec2();

  SharedPtr<Spectral> spectral_fs( new Spectral(trans) );

  SharedPtr<field::Field> surface_scalar_field( spectral_fs->createField<double>("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2 );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  array::ArrayView<double,1> surface_scalar( *surface_scalar_field );

  size_t surface_scalar_shape[] = { nspec2 };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  SharedPtr<field::Field> columns_scalar_field( spectral_fs->createField<double>("scalar",nb_levels) );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  array::ArrayView<double,2> columns_scalar( *columns_scalar_field );

  size_t columns_scalar_shape[] = { nspec2, nb_levels };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}
BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_global )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2g = trans.nspec2g();

  SharedPtr<Spectral> spectral_fs( new Spectral(trans) );

  SharedPtr<field::Field> surface_scalar_field( spectral_fs->createField<double>("scalar",field::global()) );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  BOOST_CHECK_EQUAL( surface_scalar_field->metadata().get<bool>("global"), true );

  BOOST_CHECK_EQUAL( surface_scalar_field->metadata().get<size_t>("owner"), 0 );

  array::ArrayView<double,1> surface_scalar( *surface_scalar_field );

  size_t surface_scalar_shape[] = { nspec2g };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  SharedPtr<field::Field> columns_scalar_field( spectral_fs->createField<double>("scalar",nb_levels,field::global()) );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2g*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  array::ArrayView<double,2> columns_scalar( *columns_scalar_field );

  size_t columns_scalar_shape[] = { nspec2g, nb_levels };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}
#endif


} // namespace test
} // namespace atlas
