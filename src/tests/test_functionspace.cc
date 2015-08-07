/*
 * (C) Copyright 1996-2014 ECMWF.
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
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"
#include "atlas/functionspace/NodesFunctionSpace.h"
#include "atlas/functionspace/SpectralFunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/Grid.h"
#include "atlas/Field.h"
#include "atlas/grids/ReducedGaussianGrid.h"
#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif
using namespace eckit;
using namespace atlas::functionspace;

namespace atlas {
namespace test {

template<class T>
inline std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
    return eckit::__print_list(s, v);
}

struct AtlasFixture {
    AtlasFixture()  { atlas_init(); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture )

BOOST_AUTO_TEST_CASE( test_NodesFunctionSpace )
{
  //ScopedPtr<Grid> grid( Grid::create("oct.N2") );

  size_t nlat = 2;
  long nlon[] = {4,8};
  ScopedPtr<Grid> grid( new grids::ReducedGaussianGrid( nlat, nlon ) );

  Mesh mesh;
  meshgen::ReducedGridMeshGenerator generator;
  //generator.options.set("three_dimensional",true);
  generator.generate(*grid,mesh);

  //grid.reset();

  NodesFunctionSpace nodes_fs("nodes",mesh,Halo(1));

  size_t nb_levels = 10;
  NodesColumnFunctionSpace columns_fs("columns",mesh,nb_levels,Halo(1));

  BOOST_CHECK_EQUAL( nodes_fs.nb_nodes() , columns_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( columns_fs.nb_levels() , 10 );


  size_t nb_nodes = nodes_fs.nb_nodes();

  ScopedPtr<Field> surface_scalar_field( nodes_fs.createField<double>("scalar") );
  ScopedPtr<Field> surface_vector_field( nodes_fs.createField<double>("vector",make_shape(2)) );
  ScopedPtr<Field> surface_tensor_field( nodes_fs.createField<double>("tensor",make_shape(2,2)) );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( surface_vector_field->name() , std::string("vector") );
  BOOST_CHECK_EQUAL( surface_tensor_field->name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nodes_fs.nb_nodes() );
  BOOST_CHECK_EQUAL( surface_vector_field->size() , nodes_fs.nb_nodes()*2 );
  BOOST_CHECK_EQUAL( surface_tensor_field->size() , nodes_fs.nb_nodes()*2*2 );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );
  BOOST_CHECK_EQUAL( surface_vector_field->rank() , 2 );
  BOOST_CHECK_EQUAL( surface_tensor_field->rank() , 3 );

  ArrayView<double,1> surface_scalar( *surface_scalar_field );
  ArrayView<double,2> surface_vector( *surface_vector_field );
  ArrayView<double,3> surface_tensor( *surface_tensor_field );

  size_t surface_scalar_shape[] = { nodes_fs.nb_nodes() };
  size_t surface_vector_shape[] = { nodes_fs.nb_nodes(), 2 };
  size_t surface_tensor_shape[] = { nodes_fs.nb_nodes(), 2, 2 };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_vector.shape(),surface_vector.shape()+2, surface_vector_shape,surface_vector_shape+2 );
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_tensor.shape(),surface_tensor.shape()+3, surface_tensor_shape,surface_tensor_shape+3 );


  ScopedPtr<Field> columns_scalar_field( columns_fs.createField<double>("scalar") );
  ScopedPtr<Field> columns_vector_field( columns_fs.createField<double>("vector",make_shape(2)) );
  ScopedPtr<Field> columns_tensor_field( columns_fs.createField<double>("tensor",make_shape(2,2)) );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );
  BOOST_CHECK_EQUAL( columns_vector_field->name() , std::string("vector") );
  BOOST_CHECK_EQUAL( columns_tensor_field->name() , std::string("tensor") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , columns_fs.nb_nodes()*columns_fs.nb_levels() );
  BOOST_CHECK_EQUAL( columns_vector_field->size() , columns_fs.nb_nodes()*columns_fs.nb_levels()*2 );
  BOOST_CHECK_EQUAL( columns_tensor_field->size() , columns_fs.nb_nodes()*columns_fs.nb_levels()*2*2 );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );
  BOOST_CHECK_EQUAL( columns_vector_field->rank() , 3 );
  BOOST_CHECK_EQUAL( columns_tensor_field->rank() , 4 );

  ArrayView<double,2> columns_scalar( *columns_scalar_field );
  ArrayView<double,3> columns_vector( *columns_vector_field );
  ArrayView<double,4> columns_tensor( *columns_tensor_field );

  size_t columns_scalar_shape[] = { columns_fs.nb_nodes(), columns_fs.nb_levels() };
  size_t columns_vector_shape[] = { columns_fs.nb_nodes(), columns_fs.nb_levels(), 2 };
  size_t columns_tensor_shape[] = { columns_fs.nb_nodes(), columns_fs.nb_levels(), 2, 2 };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_vector.shape(),columns_vector.shape()+3, columns_vector_shape,columns_vector_shape+3);
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_tensor.shape(),columns_tensor.shape()+4, columns_tensor_shape,columns_tensor_shape+4);

  Field::Ptr field( columns_fs.createField<int>("partition") );
  ArrayView<int,2> arr(*field);
  arr = eckit::mpi::rank();
  //field->dump( eckit::Log::info() );
  columns_fs.haloExchange(*field);
  //field->dump( eckit::Log::info() );

  Field::Ptr field2( columns_fs.createField<int>("partition2", make_shape(2) ) );
  ArrayView<int,3> arr2(*field2);
  arr2 = eckit::mpi::rank();
  //field2->dump( eckit::Log::info() );
  columns_fs.haloExchange(*field2);
  //field2->dump( eckit::Log::info() );

  Log::info() << nodes_fs.checksum(*field) << std::endl;

  Field::Ptr glb_field( nodes_fs.createGlobalField("partition",*field) );
  columns_fs.gather(*field,*glb_field);

  Log::info() << "local points = " << nodes_fs.nb_nodes() << std::endl;
  Log::info() << "grid points = " << grid->npts() << std::endl;
  Log::info() << "glb_field.shape(0) = " << glb_field->shape(0) << std::endl;

  //glb_field->dump( eckit::Log::info() );

  arr = -1;
  columns_fs.scatter(*glb_field,*field);
  columns_fs.haloExchange(*field);
  //field->dump( eckit::Log::info() );

  Log::info() << columns_fs.checksum(*field) << std::endl;

  FieldSet fields;
  fields.add(*field);
  fields.add(*field2);
  Log::info() << columns_fs.checksum(fields) << std::endl;




  BOOST_CHECKPOINT("Testing collectives for nodes scalar field");
  {
  double max;
  double min;
  double sum;
  double mean;
  double stddev;
  size_t N;
  gidx_t gidx_max;
  gidx_t gidx_min;

  ArrayView<double,1> sfc_arr( *surface_scalar_field );
  sfc_arr = eckit::mpi::rank()+1;
  nodes_fs.maximum(*surface_scalar_field,max);
  BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );

  nodes_fs.minimum(*surface_scalar_field,min);
  BOOST_CHECK_EQUAL( min, 1 );

  nodes_fs.maximum_and_global_index(*surface_scalar_field,max,gidx_max);
  BOOST_CHECK_EQUAL( max, double(eckit::mpi::size()) );
  Log::info() << "global index for maximum: " << gidx_max << std::endl;

  nodes_fs.minimum_and_global_index(*surface_scalar_field,min,gidx_min);
  BOOST_CHECK_EQUAL( min, 1 );
  Log::info() << "global index for minimum: " << gidx_min << std::endl;

  nodes_fs.sum(*surface_scalar_field,sum,N);
  Log::info() << "sum: " << sum << std::endl;
  Log::info() << "N: " << N << std::endl;

  nodes_fs.mean(*surface_scalar_field,mean,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "N: " << N << std::endl;

  nodes_fs.mean_and_standard_deviation(*surface_scalar_field,mean,stddev,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "standard deviation: " << stddev << std::endl;
  Log::info() << "N: " << N << std::endl;
  }


  BOOST_CHECKPOINT("Testing collectives for nodes vector field");
  {
  std::vector<double> max;
  std::vector<double> min;
  std::vector<double> sum;
  std::vector<double> mean;
  std::vector<double> stddev;
  size_t N;
  std::vector<gidx_t> gidx_max;
  std::vector<gidx_t> gidx_min;

  ArrayView<double,2> vec_arr( *surface_vector_field );
  vec_arr = eckit::mpi::rank()+1;
  nodes_fs.maximum(*surface_vector_field,max);
  std::vector<double> check_max(surface_vector_field->stride(0),eckit::mpi::size());
  BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );

  nodes_fs.minimum(*surface_vector_field,min);
  std::vector<double> check_min(surface_vector_field->stride(0),1);
  BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );

  nodes_fs.maximum_and_global_index(*surface_vector_field,max,gidx_max);
  BOOST_CHECK_EQUAL_COLLECTIONS( max.begin(),max.end(), check_max.begin(), check_max.end() );
  Log::info() << "global index for maximum: " << gidx_max << std::endl;

  nodes_fs.minimum_and_global_index(*surface_vector_field,min,gidx_min);
  BOOST_CHECK_EQUAL_COLLECTIONS( min.begin(),min.end(), check_min.begin(), check_min.end() );
  Log::info() << "global index for minimum: " << gidx_min << std::endl;

  nodes_fs.sum(*surface_vector_field,sum,N);
  Log::info() << "sum: " << sum << std::endl;
  Log::info() << "N: " << N << std::endl;

  nodes_fs.mean(*surface_vector_field,mean,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "N: " << N << std::endl;

  nodes_fs.mean_and_standard_deviation(*surface_vector_field,mean,stddev,N);
  Log::info() << "mean: " << mean << std::endl;
  Log::info() << "standard deviation: " << stddev << std::endl;
  Log::info() << "N: " << N << std::endl;
  }


}

BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace )
{
  size_t truncation = 159;
  size_t nb_levels = 10;
  size_t nspec2g = (truncation+1)*(truncation+2);

  SpectralFunctionSpace spectral_fs("nodes",truncation);
  SpectralColumnFunctionSpace columns_fs("columns",truncation,nb_levels);

  ScopedPtr<Field> surface_scalar_field( spectral_fs.createField("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  ArrayView<double,1> surface_scalar( *surface_scalar_field );

  size_t surface_scalar_shape[] = { nspec2g };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  ScopedPtr<Field> columns_scalar_field( columns_fs.createField("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2g*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  ArrayView<double,2> columns_scalar( *columns_scalar_field );

  size_t columns_scalar_shape[] = { nspec2g, nb_levels };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}

#ifdef ATLAS_HAVE_TRANS

BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_dist )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2 = trans.nspec2();
  size_t nspec2g = trans.nspec2g();

  SpectralFunctionSpace spectral_fs("nodes",trans);
  SpectralColumnFunctionSpace columns_fs("columns",trans,nb_levels);

  ScopedPtr<Field> surface_scalar_field( spectral_fs.createField("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2 );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  ArrayView<double,1> surface_scalar( *surface_scalar_field );

  size_t surface_scalar_shape[] = { nspec2 };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  ScopedPtr<Field> columns_scalar_field( columns_fs.createField("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  ArrayView<double,2> columns_scalar( *columns_scalar_field );

  size_t columns_scalar_shape[] = { nspec2, nb_levels };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}
BOOST_AUTO_TEST_CASE( test_SpectralFunctionSpace_trans_global )
{
  trans::Trans trans(80,159);
  size_t nb_levels(10);

  size_t nspec2 = trans.nspec2();
  size_t nspec2g = trans.nspec2g();

  SpectralFunctionSpace spectral_fs("nodes",trans);
  SpectralColumnFunctionSpace columns_fs("columns",trans,nb_levels);

  ScopedPtr<Field> surface_scalar_field( spectral_fs.createGlobalField("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( surface_scalar_field->size() , nspec2g );

  BOOST_CHECK_EQUAL( surface_scalar_field->rank() , 1 );

  ArrayView<double,1> surface_scalar( *surface_scalar_field );

  size_t surface_scalar_shape[] = { nspec2g };
  BOOST_CHECK_EQUAL_COLLECTIONS( surface_scalar.shape(),surface_scalar.shape()+1, surface_scalar_shape,surface_scalar_shape+1 );

  ScopedPtr<Field> columns_scalar_field( columns_fs.createGlobalField("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->name() , std::string("scalar") );

  BOOST_CHECK_EQUAL( columns_scalar_field->size() , nspec2g*nb_levels );

  BOOST_CHECK_EQUAL( columns_scalar_field->rank() , 2 );

  ArrayView<double,2> columns_scalar( *columns_scalar_field );

  size_t columns_scalar_shape[] = { nspec2g, nb_levels };
  BOOST_CHECK_EQUAL_COLLECTIONS(columns_scalar.shape(),columns_scalar.shape()+2, columns_scalar_shape,columns_scalar_shape+2);

}
#endif


} // namespace test
} // namespace atlas
