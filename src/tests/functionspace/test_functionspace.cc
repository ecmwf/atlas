/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "atlas/library/Library.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/grid/Grid.h"
#include "atlas/field/Field.h"
#include "atlas/parallel/mpi/mpi.h"
#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif

#include "tests/AtlasFixture.h"


using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_functionspace_NodeColumns_no_halo )
{
  Grid grid("O8");
  Mesh mesh = meshgenerator::StructuredMeshGenerator().generate(grid);
  functionspace::NodeColumns nodes_fs(mesh);
  Field field( nodes_fs.createField<int>("field") );
  array::ArrayView<int,1> value = array::make_view<int,1>( field );
  array::ArrayView<int,1> ghost = array::make_view<int,1>( mesh.nodes().ghost() );
  const size_t nb_nodes = mesh.nodes().size();
  for( size_t j=0; j<nb_nodes; ++j )
  {
    if( ghost(j) )
    {
      value(j) = -1;
    }
    else
    {
      value(j) = 1;
    }
  }
  nodes_fs.haloExchange(field);
  for( size_t j=0; j<nb_nodes; ++j )
  {
    BOOST_CHECK_EQUAL( value(j), 1 );
  }
}

BOOST_AUTO_TEST_CASE( test_functionspace_NodeColumns )
{
  //ScopedPtr<grid::Grid> grid( Grid::create("O2") );

  grid::ReducedGaussianGrid grid( {4,8,8,4} );

  meshgenerator::StructuredMeshGenerator generator;
  //generator.options.set("3d",true);
  Mesh mesh = generator.generate(grid);

  //grid.reset();

  functionspace::NodeColumns nodes_fs(mesh,mesh::Halo(1));
  size_t nb_levels = 10;
  //NodesColumnFunctionSpace columns_fs("columns",mesh,nb_levels,Halo(1));

  //BOOST_CHECK_EQUAL( nodes_fs.nb_nodes() , columns_fs.nb_nodes() );
  //BOOST_CHECK_EQUAL( columns_fs.nb_levels() , 10 );


  Field surface_scalar_field = nodes_fs.createField<double>("scalar");
  Field surface_vector_field = nodes_fs.createField<double>("vector",array::make_shape(2));
  Field surface_tensor_field = nodes_fs.createField<double>("tensor",array::make_shape(2,2));

  BOOST_CHECK_EQUAL( surface_scalar_field.name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( surface_vector_field.name() , std::string("vector") );
  BOOST_CHECK_EQUAL( surface_tensor_field.name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( surface_scalar_field.size() , nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector_field.size() , nodes_fs.nb_nodes()*2 );
  BOOST_CHECK_EQUAL( surface_tensor_field.size() , nodes_fs.nb_nodes()*2*2 );

  BOOST_CHECK_EQUAL( surface_scalar_field.rank() , 1 );
  BOOST_CHECK_EQUAL( surface_vector_field.rank() , 2 );
  BOOST_CHECK_EQUAL( surface_tensor_field.rank() , 3 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( surface_scalar_field );
  array::ArrayView<double,2> surface_vector = array::make_view<double,2>( surface_vector_field );
  array::ArrayView<double,3> surface_tensor = array::make_view<double,3>( surface_tensor_field );

  BOOST_CHECK_EQUAL( surface_scalar.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( surface_tensor.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector.shape(1), 2 );
  BOOST_CHECK_EQUAL( surface_tensor.shape(1), 2 );
  BOOST_CHECK_EQUAL( surface_tensor.shape(2), 2 );

  Field columns_scalar_field = nodes_fs.createField<double>("scalar",nb_levels);
  Field columns_vector_field = nodes_fs.createField<double>("vector",nb_levels,array::make_shape(2));
  Field columns_tensor_field = nodes_fs.createField<double>("tensor",nb_levels,array::make_shape(2,2));

  BOOST_CHECK_EQUAL( columns_scalar_field.name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( columns_vector_field.name() , std::string("vector") );
  BOOST_CHECK_EQUAL( columns_tensor_field.name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( columns_scalar_field.size() , nodes_fs.nb_nodes()*nb_levels );
  BOOST_CHECK_EQUAL( columns_vector_field.size() , nodes_fs.nb_nodes()*nb_levels*2 );
  BOOST_CHECK_EQUAL( columns_tensor_field.size() , nodes_fs.nb_nodes()*nb_levels*2*2 );

  BOOST_CHECK_EQUAL( columns_scalar_field.rank() , 2 );
  BOOST_CHECK_EQUAL( columns_vector_field.rank() , 3 );
  BOOST_CHECK_EQUAL( columns_tensor_field.rank() , 4 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( columns_scalar_field );
  array::ArrayView<double,3> columns_vector = array::make_view<double,3>( columns_vector_field );
  array::ArrayView<double,4> columns_tensor = array::make_view<double,4>( columns_tensor_field );

  BOOST_CHECK_EQUAL( columns_scalar.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( columns_vector.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( columns_tensor.shape(0), nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( columns_scalar.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_vector.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_tensor.shape(1), nb_levels );
  BOOST_CHECK_EQUAL( columns_vector.shape(2), 2 );
  BOOST_CHECK_EQUAL( columns_tensor.shape(2), 2 );
  BOOST_CHECK_EQUAL( columns_tensor.shape(3), 2 );

  Field field = nodes_fs.createField<int>("partition",nb_levels);
  array::ArrayView<int,2> arr = array::make_view<int,2>(field);
  arr.assign(parallel::mpi::comm().rank());
  //field->dump( Log::info() );
  nodes_fs.haloExchange(field);
  //field->dump( Log::info() );

  Field field2 = nodes_fs.createField<int>("partition2",nb_levels,array::make_shape(2));
  Log::info() << "field2.rank() = " << field2.rank() << std::endl;
  array::ArrayView<int,3> arr2 = array::make_view<int,3>(field2);
  arr2.assign(parallel::mpi::comm().rank());

  //field2->dump( Log::info() );
  nodes_fs.haloExchange(field2);
  //field2->dump( Log::info() );

  Log::info() << nodes_fs.checksum(field) << std::endl;

  size_t root = parallel::mpi::comm().size()-1;
  Field glb_field = nodes_fs.createField("partition",field,field::global(root));
  nodes_fs.gather(field,glb_field);

  Log::info() << "local points = " << nodes_fs.nb_nodes() << std::endl;
  Log::info() << "grid points = " << grid.size() << std::endl;
  Log::info() << "glb_field.shape(0) = " << glb_field.shape(0) << std::endl;

  BOOST_CHECK_EQUAL( glb_field.metadata().get<bool>("global"), true );
  BOOST_CHECK_EQUAL( glb_field.metadata().get<int>("owner"),  root );

  //glb_field->dump( Log::info() );

  if( parallel::mpi::comm().rank() == root )
    glb_field.metadata().set("test_broadcast",123);

  arr.assign(-1);
  nodes_fs.scatter(glb_field,field);
  BOOST_CHECK_EQUAL( field.metadata().get<int>("test_broadcast"), 123 );
  nodes_fs.haloExchange(field);
  //field->dump( Log::info() );

  Log::info() << nodes_fs.checksum(field) << std::endl;

  FieldSet fields;
  fields.add(field);
  fields.add(field2);
  Log::info() << nodes_fs.checksum(fields) << std::endl;



  Log::info() << "Testing collectives for nodes scalar field" << std::endl;
  BOOST_TEST_CHECKPOINT("Testing collectives for nodes scalar field");
  {
    const Field& field = surface_scalar_field;
    const functionspace::NodeColumns fs = nodes_fs;

  double max;
  double min;
  double sum;
  double mean;
  double stddev;
  size_t N;
  gidx_t gidx_max;
  gidx_t gidx_min;

  array::ArrayView<double,1> sfc_arr = array::make_view<double,1>( field );
  sfc_arr.assign( parallel::mpi::comm().rank()+1 );
  fs.maximum(surface_scalar_field,max);
  BOOST_CHECK_EQUAL( max, double(parallel::mpi::comm().size()) );

  fs.minimum(surface_scalar_field,min);
  BOOST_CHECK_EQUAL( min, 1 );

  fs.maximumAndLocation(field,max,gidx_max);
  BOOST_CHECK_EQUAL( max, double(parallel::mpi::comm().size()) );
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
    const Field& field = surface_vector_field;
    const functionspace::NodeColumns fs = nodes_fs;

    std::vector<double> max;
    std::vector<double> min;
    std::vector<double> sum;
    std::vector<double> mean;
    std::vector<double> stddev;
    size_t N;
    std::vector<gidx_t> gidx_max;
    std::vector<gidx_t> gidx_min;

    array::ArrayView<double,2> vec_arr = array::make_view<double,2>( field );
    vec_arr.assign( parallel::mpi::comm().rank()+1 );
    fs.maximum(field,max);
    std::vector<double> check_max(field.stride(0),parallel::mpi::comm().size());
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
    const Field& field = columns_scalar_field;
    const functionspace::NodeColumns fs = nodes_fs;
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
    arr.assign( parallel::mpi::comm().rank()+1 );
    fs.maximum(field,max);
    BOOST_CHECK_EQUAL( max, double(parallel::mpi::comm().size()) );

    fs.minimum(field,min);
    BOOST_CHECK_EQUAL( min, 1 );

    fs.maximumAndLocation(field,max,gidx_max,level);
    BOOST_CHECK_EQUAL( max, double(parallel::mpi::comm().size()) );
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

    Field max_per_level    ( "max",    array::make_datatype<double>(), array::make_shape(nb_levels) );
    Field min_per_level    ( "min",    array::make_datatype<double>(), array::make_shape(nb_levels) );
    Field sum_per_level    ( "sum",    array::make_datatype<double>(), array::make_shape(nb_levels) );
    Field mean_per_level   ( "mean",   array::make_datatype<double>(), array::make_shape(nb_levels) );
    Field stddev_per_level ( "stddev", array::make_datatype<double>(), array::make_shape(nb_levels) );
    Field gidx_per_level   ( "gidx",   array::make_datatype<gidx_t>(), array::make_shape(nb_levels) );

    fs.maximumPerLevel(field,max_per_level);
    max_per_level.dump(Log::info());
    fs.minimumPerLevel(field,min_per_level);
    min_per_level.dump(Log::info());
    fs.sumPerLevel(field,sum_per_level,N);
    sum_per_level.dump(Log::info());
    fs.meanPerLevel(field,mean_per_level,N);
    mean_per_level.dump(Log::info());
    fs.meanAndStandardDeviationPerLevel(field,mean_per_level,stddev_per_level,N);
    mean_per_level.dump(Log::info());
    stddev_per_level.dump(Log::info());
    fs.orderIndependentSumPerLevel(field,sum_per_level,N);
    sum_per_level.dump(Log::info());

  }

  BOOST_TEST_CHECKPOINT("Testing collectives for columns vector field");
  if(1){
    const Field& field = columns_vector_field;
    const functionspace::NodeColumns fs = nodes_fs;
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
    vec_arr.assign( parallel::mpi::comm().rank()+1 );
    fs.maximum(field,max);
    std::vector<double> check_max(nvar,parallel::mpi::comm().size());
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

    Field max_per_level    ( "max",    array::make_datatype<double>(), array::make_shape(nb_levels,nvar) );
    Field min_per_level    ( "min",    array::make_datatype<double>(), array::make_shape(nb_levels,nvar) );
    Field sum_per_level    ( "sum",    array::make_datatype<double>(), array::make_shape(nb_levels,nvar) );
    Field mean_per_level   ( "mean",   array::make_datatype<double>(), array::make_shape(nb_levels,nvar) );
    Field stddev_per_level ( "stddev", array::make_datatype<double>(), array::make_shape(nb_levels,nvar) );
    Field gidx_per_level   ( "gidx",   array::make_datatype<gidx_t>(), array::make_shape(nb_levels,nvar) );

    fs.maximumPerLevel(field,max_per_level);
    max_per_level.dump(Log::info());

    fs.minimumPerLevel(field,min_per_level);
    min_per_level.dump(Log::info());

    fs.sumPerLevel(field,sum_per_level,N);
    sum_per_level.dump(Log::info());

    fs.meanPerLevel(field,mean_per_level,N);
    mean_per_level.dump(Log::info());

    fs.meanAndStandardDeviationPerLevel(field,mean_per_level,stddev_per_level,N);
    mean_per_level.dump(Log::info());
    stddev_per_level.dump(Log::info());

    fs.orderIndependentSumPerLevel(field,sum_per_level,N);
    sum_per_level.dump(Log::info());
  }



  Field tmp =
        nodes_fs.createField( field::datatypeT<double>() | field::global(0) | field::levels(10) | field::name("tmp") );

}


BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace )
{
  size_t truncation = 159;
  size_t nb_levels = 10;
  size_t nspec2g = (truncation+1)*(truncation+2);

  Spectral spectral_fs(truncation);

  Field surface_scalar_field = spectral_fs.createField<double>("scalar");

  BOOST_CHECK_EQUAL( surface_scalar_field.name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field.size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field.rank() , 1 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( surface_scalar_field );

  BOOST_CHECK_EQUAL( surface_scalar.shape(0), nspec2g );

  Field columns_scalar_field = spectral_fs.createField<double>("scalar",nb_levels);

  BOOST_CHECK_EQUAL( columns_scalar_field.name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field.size() , nspec2g*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field.rank() , 2 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( columns_scalar_field );

  BOOST_CHECK_EQUAL( columns_scalar.shape(0), nspec2g );
  BOOST_CHECK_EQUAL( columns_scalar.shape(1), nb_levels );

}


#ifdef ATLAS_HAVE_TRANS

BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_dist )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2 = trans.nspec2();

  Spectral spectral_fs( trans );

  Field surface_scalar_field = spectral_fs.createField<double>("scalar");

  BOOST_CHECK_EQUAL( surface_scalar_field.name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field.size() , nspec2 );

  BOOST_CHECK_EQUAL( surface_scalar_field.rank() , 1 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( surface_scalar_field );

  BOOST_CHECK_EQUAL( surface_scalar.shape(0), nspec2 );
  // size_t surface_scalar_shape[] = { nspec2 };
  // BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  Field columns_scalar_field = spectral_fs.createField<double>("scalar",nb_levels);

  BOOST_CHECK_EQUAL( columns_scalar_field.name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field.size() , nspec2*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field.rank() , 2 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( columns_scalar_field );

  BOOST_CHECK_EQUAL(columns_scalar.shape(0), nspec2);
  BOOST_CHECK_EQUAL(columns_scalar.shape(1), nb_levels);
  // size_t columns_scalar_shape[] = { nspec2, nb_levels };
  // BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}
BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_global )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2g = trans.nspec2g();

  Spectral spectral_fs( trans );

  Field surface_scalar_field = spectral_fs.createField<double>("scalar",field::global());

  BOOST_CHECK_EQUAL( surface_scalar_field.name() , std::string("scalar") );

  if( eckit::mpi::comm().rank() == 0 )
    BOOST_CHECK_EQUAL( surface_scalar_field.size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field.rank() , 1 );

  BOOST_CHECK_EQUAL( surface_scalar_field.metadata().get<bool>("global"), true );

  BOOST_CHECK_EQUAL( surface_scalar_field.metadata().get<size_t>("owner"), 0 );

  array::ArrayView<double,1> surface_scalar = array::make_view<double,1>( surface_scalar_field );

  if( eckit::mpi::comm().rank() == 0 ) {
    BOOST_CHECK_EQUAL( surface_scalar.shape(0), nspec2g );
  }
  Field columns_scalar_field = spectral_fs.createField<double>("scalar",nb_levels,field::global());

  BOOST_CHECK_EQUAL( columns_scalar_field.name() , std::string("scalar") );

  if( eckit::mpi::comm().rank() == 0 ) {
    BOOST_CHECK_EQUAL( columns_scalar_field.size() , nspec2g*nb_levels );
  } else {
    BOOST_CHECK_EQUAL( columns_scalar_field.size() , 0 );
  }

  BOOST_CHECK_EQUAL( columns_scalar_field.rank() , 2 );

  array::ArrayView<double,2> columns_scalar = array::make_view<double,2>( columns_scalar_field );

  if( eckit::mpi::comm().rank() == 0 ) {
    BOOST_CHECK_EQUAL( columns_scalar.shape(0), nspec2g );
    BOOST_CHECK_EQUAL( columns_scalar.shape(1), nb_levels );
  }
}
BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_norm )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  Spectral spectral_fs( trans );

  Field twoD_field   = spectral_fs.createField<double>("2d");
  Field threeD_field = spectral_fs.createField<double>("3d",nb_levels);

  // Set first wave number
  {
    array::ArrayView<double,1> twoD = array::make_view<double,1>( twoD_field );
    twoD.assign(0.);
    if( parallel::mpi::comm().rank() == 0 ) twoD(0) = 1.;

    array::ArrayView<double,2> threeD = array::make_view<double,2>( threeD_field );
    threeD.assign(0.);
    for( size_t jlev=0; jlev<nb_levels; ++jlev) {
      if( parallel::mpi::comm().rank() == 0 ) threeD(0,jlev) = jlev;
    }
  }

  double twoD_norm(0.);
  std::vector<double> threeD_norms(threeD_field.levels(),0.);

  spectral_fs.norm(twoD_field,twoD_norm);
  spectral_fs.norm(threeD_field,threeD_norms);

  if( eckit::mpi::comm().rank() == 0 ) {
    BOOST_CHECK_CLOSE( twoD_norm, 1.0, 1.e-10 );
    for( size_t jlev=0; jlev<nb_levels; ++jlev) {
      BOOST_CHECK_CLOSE( threeD_norms[jlev], double(jlev), 1.e-10 );
    }
  }
}

#endif

} // namespace test
} // namespace atlas
