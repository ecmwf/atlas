/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>

#include "eckit/utils/MD5.h"

#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/MakeView.h"

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {
namespace functionspace {
namespace detail {

namespace {

void set_field_metadata(const eckit::Parametrisation& config, Field& field)
{
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      field.metadata().set("owner",owner);
    }
  }
  field.metadata().set("global",global);
}

}

size_t EdgeColumns::config_size(const eckit::Parametrisation& config) const
{
  size_t size = nb_edges();
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      size = (parallel::mpi::comm().rank() == owner ? nb_edges_global() : 0);
    }
  }
  return size;
}

EdgeColumns::EdgeColumns( const mesh::Mesh& mesh )
  : mesh_(mesh),
    edges_(mesh_.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  constructor();
}

EdgeColumns::EdgeColumns( const mesh::Mesh& mesh, const mesh::Halo &halo, const eckit::Parametrisation &params ) :
    mesh_(mesh),
    edges_(mesh_.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  size_t mesh_halo_size_;
  mesh.metadata().get("halo",mesh_halo_size_);
  ASSERT( mesh_halo_size_ == halo.size() );
  constructor();
}

EdgeColumns::EdgeColumns( const mesh::Mesh& mesh, const mesh::Halo &halo) :
    mesh_(mesh),
    edges_(mesh_.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  size_t mesh_halo_size_;
  mesh.metadata().get("halo",mesh_halo_size_);
  ASSERT( mesh_halo_size_ == halo.size() );
  constructor();
}


void EdgeColumns::constructor()
{
  nb_edges_ = mesh().edges().size();

  gather_scatter_.reset(new parallel::GatherScatter());
  halo_exchange_.reset(new parallel::HaloExchange());
  checksum_.reset(new parallel::Checksum());

  const Field& partition    = edges().partition();
  const Field& remote_index = edges().remote_index();
  const Field& global_index = edges().global_index();

  halo_exchange_->setup(
        array::make_view<int,1>(partition).data(),
        array::make_view<int,1>(remote_index).data(),REMOTE_IDX_BASE,
        nb_edges_);

  gather_scatter_->setup(
        array::make_view<int,1>(partition).data(),
        array::make_view<int,1>(remote_index).data(),REMOTE_IDX_BASE,
        array::make_view<gidx_t,1>(global_index).data(),
        nb_edges_);

  checksum_->setup(
        array::make_view<int,1>(partition).data(),
        array::make_view<int,1>(remote_index).data(),REMOTE_IDX_BASE,
        array::make_view<gidx_t,1>(global_index).data(),
        nb_edges_);

  nb_edges_global_ =  gather_scatter_->glb_dof();
}

EdgeColumns::~EdgeColumns() {}

size_t EdgeColumns::footprint() const {
  size_t size = sizeof(*this);
  // TODO
  return size;
}

size_t EdgeColumns::nb_edges() const
{
  return nb_edges_;
}

size_t EdgeColumns::nb_edges_global() const
{
  return nb_edges_global_;
}

Field EdgeColumns::createField(const std::string& name,array::DataType datatype,const eckit::Parametrisation& options) const
{
  size_t nb_edges = config_size(options);
  Field field = Field(name,datatype,array::make_shape(nb_edges));
  field.set_functionspace(this);
  set_field_metadata(options,field);
  return field;
}

Field EdgeColumns::createField(const std::string& name,array::DataType datatype, size_t levels,const eckit::Parametrisation& options) const
{
  size_t nb_edges = config_size(options);
  Field field = Field(name,datatype,array::make_shape(nb_edges,levels));
  field.set_levels(levels);
  field.set_functionspace(this);
  set_field_metadata(options,field);
  return field;
}

Field EdgeColumns::createField(const std::string& name,array::DataType datatype, const std::vector<size_t>& variables,const eckit::Parametrisation& options) const
{
  size_t nb_edges = config_size(options);
  std::vector<size_t> shape(1,nb_edges);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field field = Field(name,datatype,shape);
  field.set_functionspace(this);
  set_field_metadata(options,field);
  return field;
}

Field EdgeColumns::createField(const std::string& name, array::DataType datatype, size_t levels, const std::vector<size_t>& variables,const eckit::Parametrisation& options) const
{
  size_t nb_edges = config_size(options);
  std::vector<size_t> shape(1,nb_edges); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field field = Field(name,datatype,shape);
  field.set_levels(levels);
  field.set_functionspace(this);
  set_field_metadata(options,field);
  return field;
}

Field EdgeColumns::createField(const std::string& name, const Field& other,const eckit::Parametrisation& options) const {
  size_t nb_edges = config_size(options);
  array::ArrayShape shape = other.shape();
  shape[0] = nb_edges;
  Field field = Field(name,other.datatype(),shape);
  if( other.has_levels() )
    field.set_levels(field.shape(1));
  field.set_functionspace(this);
  set_field_metadata(options,field);
  return field;
}

void EdgeColumns::haloExchange( FieldSet& fieldset ) const
{
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field = fieldset[f];
    if     ( field.datatype() == array::DataType::kind<int>() ) {
      array::ArrayView<int,2> view = array::make_view<int,2>(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<long>() ) {
      array::ArrayView<long,2> view = array::make_view<long,2>(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<float>() ) {
      array::ArrayView<float,2> view = array::make_view<float,2>(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<double>() ) {
      array::ArrayView<double,2> view = array::make_view<double,2>(field);
      halo_exchange().execute( view );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void EdgeColumns::haloExchange( Field& field ) const
{
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
}
const parallel::HaloExchange& EdgeColumns::halo_exchange() const
{
  return *halo_exchange_;
}


void EdgeColumns::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];
    const size_t nb_fields = 1;
    size_t root(0);
    glb.metadata().get("owner",root);
    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
      parallel::Field<int      > glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
      parallel::Field<long      > glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
      parallel::Field<float      > glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
      parallel::Field<double      > glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void EdgeColumns::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}
const parallel::GatherScatter& EdgeColumns::gather() const
{
  return *gather_scatter_;
}
const parallel::GatherScatter& EdgeColumns::scatter() const
{
  return *gather_scatter_;
}


void EdgeColumns::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    const size_t nb_fields = 1;
    size_t root(0);
    glb.metadata().get("owner",root);
    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
      parallel::Field<int      > loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
      parallel::Field<long      > loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
      parallel::Field<float      > loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
      parallel::Field<double      > loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else throw eckit::Exception("datatype not supported",Here());
    glb.metadata().broadcast(loc.metadata(),root);
    loc.metadata().set("global",false);
  }
}
void EdgeColumns::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const parallel::Checksum& checksum, const Field& field )
{
  array::ArrayView<T,3> values = array::make_view<T,3>(field);
  array::ArrayT<T> surface_field( field.shape(0),field.shape(2) );
  array::ArrayView<T,2> surface = array::make_view<T,2>(surface_field);
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<surface.shape(1); ++j )
    {
      surface(n,j) = 0.;
      for( size_t l=0; l<values.shape(1);++l )
        surface(n,j) += values(n,l,j);
    }
  }
  return checksum.execute( surface.data(), surface_field.stride(0) );
}
template <typename T>
std::string checksum_2d_field(const parallel::Checksum& checksum, const Field& field )
{
  array::ArrayView<T,2> values = array::make_view<T,2>(field);
  return checksum.execute( values.data(), field.stride(0) );
}

}

std::string EdgeColumns::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field=fieldset[f];
    if     ( field.datatype() == array::DataType::kind<int>() ) {
      if( field.has_levels() )
        md5 << checksum_3d_field<int>(checksum(),field);
      else
        md5 << checksum_2d_field<int>(checksum(),field);
    }
    else if( field.datatype() == array::DataType::kind<long>() ) {
      if( field.has_levels() )
        md5 << checksum_3d_field<long>(checksum(),field);
      else
        md5 << checksum_2d_field<long>(checksum(),field);
    }
    else if( field.datatype() == array::DataType::kind<float>() ) {
      if( field.has_levels() )
        md5 << checksum_3d_field<float>(checksum(),field);
      else
        md5 << checksum_2d_field<float>(checksum(),field);
    }
    else if( field.datatype() == array::DataType::kind<double>() ) {
      if( field.has_levels() )
        md5 << checksum_3d_field<double>(checksum(),field);
      else
        md5 << checksum_2d_field<double>(checksum(),field);
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string EdgeColumns::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

const parallel::Checksum& EdgeColumns::checksum() const
{
  return *checksum_;
}

namespace {
void reverse_copy(const int variables[], const int size, std::vector<size_t> &reverse)
{
  int r=size;
  for( int i=0; i<size; ++i)
  {
    reverse[--r] = static_cast<size_t>(variables[i]);
  }
}

void copy(const int variables[], const int size, std::vector<size_t> &cpy)
{
  for( int i=0; i<size; ++i)
  {
    cpy[i] = static_cast<size_t>(variables[i]);
  }
}

std::vector<size_t> variables_to_vector( const int variables[], const int size, bool fortran_ordering )
{
  std::vector<size_t> vec(size);
  if (fortran_ordering)
    reverse_copy(variables,size,vec);
  else
    copy(variables,size,vec);
  return vec;
}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
extern "C" {
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


EdgeColumns* atlas__functionspace__Edges__new ( mesh::Mesh::Implementation* mesh, int halo )
{
  EdgeColumns* edges(0);
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      mesh::Mesh m(mesh);
      edges = new EdgeColumns(m,mesh::Halo(halo));
  );
  return edges;
}

//------------------------------------------------------------------------------

EdgeColumns* atlas__functionspace__Edges__new_mesh ( mesh::Mesh::Implementation* mesh )
{
  EdgeColumns* edges(0);
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      mesh::Mesh m(mesh);
      edges = new EdgeColumns(m);
  );
  return edges;
}

//------------------------------------------------------------------------------

void atlas__functionspace__Edges__delete (EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete(This);
  );
}

//------------------------------------------------------------------------------

int atlas__functionspace__Edges__nb_edges(const EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->nb_edges();
  );
  return 0;
}

//------------------------------------------------------------------------------

mesh::Mesh::Implementation* atlas__functionspace__Edges__mesh(EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return This->mesh().get();
  );
  return 0;
}

//------------------------------------------------------------------------------

mesh::Edges* atlas__functionspace__Edges__edges(EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
      ASSERT(This);
      return &This->edges();
  );
  return 0;
}

//------------------------------------------------------------------------------

using field::FieldImpl;
using field::FieldSetImpl;

field::FieldImpl* atlas__functionspace__Edges__create_field (
    const EdgeColumns* This,
    const char* name,
    int kind,
    const eckit::Parametrisation* options )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(options);
    FieldImpl* field;
    {
      Field f = This->createField(std::string(name),array::DataType(kind));
      field = f.get();
      field->attach();
    }
    field->detach();
    return field
  );
  return 0;
}

//------------------------------------------------------------------------------

field::FieldImpl* atlas__functionspace__Edges__create_field_vars (
    const EdgeColumns* This,
    const char* name,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind,
    const eckit::Parametrisation* options )
{

  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    ASSERT(options);
    FieldImpl* field;
    {
      Field f = This->createField(
        std::string(name),
        array::DataType(kind),
        variables_to_vector(variables,variables_size,fortran_ordering) );
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::FieldImpl* atlas__functionspace__Edges__create_field_lev (
    const EdgeColumns* This,
    const char* name,
    int levels,
    int kind,
    const eckit::Parametrisation* options )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(options);
    FieldImpl* field;
    {
      Field f = This->createField(std::string(name),array::DataType(kind),size_t(levels));
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::FieldImpl* atlas__functionspace__Edges__create_field_lev_vars (
    const EdgeColumns* This,
    const char* name,
    int levels,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind,
    const eckit::Parametrisation* options )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    ASSERT(options);
    FieldImpl* field;
    {
      Field f = This->createField(
        std::string(name),
        array::DataType(kind),
        size_t(levels),
        variables_to_vector(variables,variables_size,fortran_ordering) );
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::FieldImpl* atlas__functionspace__Edges__create_field_template (
    const EdgeColumns* This,
    const char* name,
    const field::FieldImpl* field_template,
    const eckit::Parametrisation* options )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(options);
    FieldImpl* field;
    {
      Field f = This->createField(std::string(name),field_template);
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_fieldset(
    const EdgeColumns* This,
    field::FieldSetImpl* fieldset)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    FieldSet f(fieldset);
    This->haloExchange(f);
  );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_field(const EdgeColumns* This, field::FieldImpl* field)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(field);
        Field f(field);
        This->haloExchange(f);
   );
}

// -----------------------------------------------------------------------------------

const parallel::HaloExchange* atlas__functionspace__Edges__get_halo_exchange(const EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->halo_exchange();
  );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__gather_fieldset(
    const EdgeColumns* This,
    const field::FieldSetImpl* local,
    field::FieldSetImpl* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        const FieldSet l(local);
        FieldSet g(global);
        This->gather(l,g); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__gather_field(
    const EdgeColumns* This,
    const field::FieldImpl* local,
    field::FieldImpl* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        const Field l(local);
        Field g(global);
        This->gather(l,g); );
}

// -----------------------------------------------------------------------------------

const parallel::GatherScatter* atlas__functionspace__Edges__get_gather(const EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->gather(); );
  return 0;
}

// -----------------------------------------------------------------------------------

const parallel::GatherScatter* atlas__functionspace__Edges__get_scatter(const EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->scatter(); );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__scatter_fieldset(const EdgeColumns* This, const field::FieldSetImpl* global, field::FieldSetImpl* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        const FieldSet g(global);
        FieldSet l(local);
        This->scatter(g,l); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__scatter_field(const EdgeColumns* This, const field::FieldImpl* global, field::FieldImpl* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(global);
        ASSERT(local);
        const Field g(global);
        Field l(local);
        This->scatter(g,l); );
}

// -----------------------------------------------------------------------------------

const parallel::Checksum* atlas__functionspace__Edges__get_checksum(const EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->checksum();
  );
  return 0;
}

// -----------------------------------------------------------------------------------


void atlas__functionspace__Edges__checksum_fieldset(
    const EdgeColumns* This,
    const field::FieldSetImpl* fieldset,
    char* &checksum,
    int &size,
    int &allocated)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    std::string checksum_str (This->checksum(fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__checksum_field(
    const EdgeColumns* This,
    const field::FieldImpl* field,
    char* &checksum,
    int &size,
    int &allocated)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    std::string checksum_str (This->checksum(field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

}

// -----------------------------------------------------------------------------------

} // namespace detail

// -----------------------------------------------------------------------------------

EdgeColumns::EdgeColumns() :
    FunctionSpace(),
    functionspace_(nullptr) {
}

EdgeColumns::EdgeColumns( const FunctionSpace& functionspace ) :
    FunctionSpace(functionspace),
    functionspace_( dynamic_cast< const detail::EdgeColumns* >( get() ) ) {
}

EdgeColumns::EdgeColumns( 
  const mesh::Mesh& mesh,
  const mesh::Halo& halo,
  const eckit::Parametrisation& config ) :
  FunctionSpace( new detail::EdgeColumns(mesh,halo,config) ),
  functionspace_( dynamic_cast< const detail::EdgeColumns* >( get() ) ) {
}

EdgeColumns::EdgeColumns( 
  const mesh::Mesh& mesh,
  const mesh::Halo& halo) :
  FunctionSpace( new detail::EdgeColumns(mesh,halo) ),
  functionspace_( dynamic_cast< const detail::EdgeColumns* >( get() ) ) {
}

EdgeColumns::EdgeColumns( 
  const mesh::Mesh& mesh) :
  FunctionSpace( new detail::EdgeColumns(mesh) ),
  functionspace_( dynamic_cast< const detail::EdgeColumns* >( get() ) ) {
}
  
size_t EdgeColumns::nb_edges() const {
  return functionspace_->nb_edges();
}

size_t EdgeColumns::nb_edges_global() const { // Only on MPI rank 0, will this be different from 0
  return functionspace_->nb_edges_global();
}

const mesh::Mesh& EdgeColumns::mesh() const {
  return functionspace_->mesh();
}

const mesh::HybridElements& EdgeColumns::edges() const {
  return functionspace_->edges();
}

Field EdgeColumns::createField(
        const std::string& name,
        array::DataType datatype,
        const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,datatype,config);
}

Field EdgeColumns::createField(
        const std::string& name,
        array::DataType datatype,
        size_t levels,
        const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,datatype,levels,config);
}

Field EdgeColumns::createField(const std::string& name,
                          array::DataType datatype,
                          const std::vector<size_t>& variables,
                          const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,datatype,variables,config);
}

Field EdgeColumns::createField(const std::string& name,
                          array::DataType datatype,
                          size_t levels,
                          const std::vector<size_t>& variables,
                          const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,datatype,levels,variables,config);
}

Field EdgeColumns::createField(const std::string& name,
                          const Field& field,
                          const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,field,config);
}



void EdgeColumns::haloExchange( FieldSet& fieldset ) const {
  functionspace_->haloExchange(fieldset);
}

void EdgeColumns::haloExchange( Field& field) const {
  functionspace_->haloExchange(field);
}

const parallel::HaloExchange& EdgeColumns::halo_exchange() const {
  return functionspace_->halo_exchange();
}

void EdgeColumns::gather( const FieldSet& local, FieldSet& global ) const {
  functionspace_->gather(local,global);
}

void EdgeColumns::gather( const Field& local, Field& global ) const {
  functionspace_->gather(local,global);
}

const parallel::GatherScatter& EdgeColumns::gather() const {
  return functionspace_->gather();
}

void EdgeColumns::scatter( const FieldSet& global, FieldSet& local ) const {
  functionspace_->scatter(global,local);
}

void EdgeColumns::scatter( const Field& global, Field& local ) const {
  functionspace_->scatter(global,local);
}

const parallel::GatherScatter& EdgeColumns::scatter() const {
  return functionspace_->scatter();
}

std::string EdgeColumns::checksum( const FieldSet& fieldset ) const {
  return functionspace_->checksum(fieldset);
}

std::string EdgeColumns::checksum( const Field& field ) const {
  return functionspace_->checksum(field);
}

const parallel::Checksum& EdgeColumns::checksum() const {
  return functionspace_->checksum();
}

} // namespace functionspace
} // namespace atlas

