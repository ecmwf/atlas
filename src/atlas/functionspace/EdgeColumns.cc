/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "eckit/os/BackTrace.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/internals/IsGhost.h"
#include "atlas/parallel/mpi/Collectives.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/runtime/Log.h"

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {
namespace functionspace {

namespace {

template <typename T>
array::ArrayView<T,3> leveled_view(const field::Field &field)
{
  if( field.has_levels() )
    return array::ArrayView<T,3> ( field.data<T>(), array::make_shape(field.shape(0),field.shape(1),field.stride(1)) );
  else
    return array::ArrayView<T,3> ( field.data<T>(), array::make_shape(field.shape(0),1,field.stride(0)) );
}

template <typename T>
array::ArrayView<T,2> surface_view(const field::Field &field)
{
  return array::ArrayView<T,2> ( field.data<T>(), array::make_shape(field.shape(0),field.stride(0)) );
}

template <typename T>
array::ArrayView<T,2> leveled_scalar_view(const field::Field &field)
{
  if( field.has_levels() )
    return array::ArrayView<T,2> ( field.data<T>(), array::make_shape(field.shape(0),field.shape(1)) );
  else
    return array::ArrayView<T,2> ( field.data<T>(), array::make_shape(field.shape(0),1) );
}

template <typename T>
array::ArrayView<T,1> surface_scalar_view(const field::Field &field)
{
  return array::ArrayView<T,1> ( field.data<T>(), array::make_shape(field.size()) );
}

void print_warning(const eckit::CodeLocation& here)
{

  Log::warning()
      << "  " << here << '\n'
      << "  Function createGlobalField is deprecated. Please use createField and pass \n"
      << "  the extra argument \n"
      << "      atlas::field::global() \n"
      << "  with owner a mpi rank"
      << std::endl;
  Log::debug()
      << eckit::BackTrace::dump()
      << std::endl;
}

}

EdgeColumns::EdgeColumns( mesh::Mesh& mesh )
  : mesh_(&mesh),
    edges_(mesh.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  constructor();
}

EdgeColumns::EdgeColumns( mesh::Mesh& mesh, const mesh::Halo &halo, const eckit::Parametrisation &params )
  : FunctionSpace(),
    mesh_(&mesh),
    edges_(mesh.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  size_t mesh_halo_size_;
  mesh.metadata().get("halo",mesh_halo_size_);
  ASSERT( mesh_halo_size_ == halo.size() );
  constructor();
}

EdgeColumns::EdgeColumns(mesh::Mesh& mesh, const mesh::Halo &halo)
  : FunctionSpace(),
    mesh_(&mesh),
    edges_(mesh.edges()),
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

  const field::Field& partition    = edges().partition();
  const field::Field& remote_index = edges().remote_index();
  const field::Field& global_index = edges().global_index();

  halo_exchange_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        nb_edges_);

  gather_scatter_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        edges_.global_index().data<gidx_t>(),
        nb_edges_);

  checksum_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        global_index.data<gidx_t>(),
        nb_edges_);

  size_t root = 0;
  nb_edges_global_ = eckit::mpi::rank() == root ? gather_scatter_->glb_dof() : 0;
}

EdgeColumns::~EdgeColumns() {}

size_t EdgeColumns::nb_edges() const
{
  return nb_edges_;
}

size_t EdgeColumns::nb_edges_global() const
{
  return nb_edges_global_;
}

field::Field* EdgeColumns::createField(const std::string& name,array::DataType datatype) const {
  field::Field* field = field::Field::create(name,datatype,array::make_shape(nb_edges()));
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createField(const std::string& name,array::DataType datatype, size_t levels) const {
  field::Field* field = field::Field::create(name,datatype,array::make_shape(nb_edges(),levels));
  field->set_levels(levels);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createField(const std::string& name,array::DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field* field = field::Field::create(name,datatype,shape);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createField(const std::string& name, array::DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field* field = field::Field::create(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createField(const std::string& name, const field::Field& other) const {
  array::ArrayShape shape = other.shape();
  shape[0] = nb_edges();
  field::Field* field = field::Field::create(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createGlobalField(const std::string& name,array::DataType datatype) const {
  field::Field* field = field::Field::create(name,datatype,array::make_shape(nb_edges_global()));
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createGlobalField(const std::string& name, array::DataType datatype, size_t levels) const {
  field::Field* field = field::Field::create(name,datatype,array::make_shape(nb_edges_global(),levels));
  field->set_levels(levels);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createGlobalField(const std::string& name, array::DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field* field = field::Field::create(name,datatype,shape);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createGlobalField(const std::string& name, array::DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges_global()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field* field = field::Field::create(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(*this);
  return field;
}

field::Field* EdgeColumns::createGlobalField(const std::string& name,const field::Field& other) const {
  array::ArrayShape shape = other.shape();
  shape[0] = nb_edges_global();
  field::Field* field = field::Field::create(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(*this);
  return field;
}

void EdgeColumns::haloExchange( field::FieldSet& fieldset ) const
{
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const field::Field& field = fieldset[f];
    if     ( field.datatype() == array::DataType::kind<int>() ) {
      array::ArrayView<int,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<long>() ) {
      array::ArrayView<long,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<float>() ) {
      array::ArrayView<float,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == array::DataType::kind<double>() ) {
      array::ArrayView<double,2> view(field);
      halo_exchange().execute( view );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void EdgeColumns::haloExchange( field::Field& field ) const
{
    field::FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
}
const parallel::HaloExchange& EdgeColumns::halo_exchange() const
{
  return *halo_exchange_;
}


void EdgeColumns::gather( const field::FieldSet& local_fieldset, field::FieldSet& global_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const field::Field& loc = local_fieldset[f];
    field::Field& glb = global_fieldset[f];
    const size_t nb_fields = 1;
    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> loc_field(loc.data<int>(),loc.stride(0));
      parallel::Field<int      > glb_field(glb.data<int>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> loc_field(loc.data<long>(),loc.stride(0));
      parallel::Field<long      > glb_field(glb.data<long>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> loc_field(loc.data<float>(),loc.stride(0));
      parallel::Field<float      > glb_field(glb.data<float>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> loc_field(loc.data<double>(),loc.stride(0));
      parallel::Field<double      > glb_field(glb.data<double>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void EdgeColumns::gather( const field::Field& local, field::Field& global ) const
{
  field::FieldSet local_fields;
  field::FieldSet global_fields;
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


void EdgeColumns::scatter( const field::FieldSet& global_fieldset, field::FieldSet& local_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const field::Field& glb = global_fieldset[f];
    field::Field& loc = local_fieldset[f];
    const size_t nb_fields = 1;

    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> glb_field(glb.data<int>(),glb.stride(0));
      parallel::Field<int      > loc_field(loc.data<int>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> glb_field(glb.data<long>(),glb.stride(0));
      parallel::Field<long      > loc_field(loc.data<long>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> glb_field(glb.data<float>(),glb.stride(0));
      parallel::Field<float      > loc_field(loc.data<float>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> glb_field(glb.data<double>(),glb.stride(0));
      parallel::Field<double      > loc_field(loc.data<double>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void EdgeColumns::scatter( const field::Field& global, field::Field& local ) const
{
  field::FieldSet global_fields;
  field::FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const parallel::Checksum& checksum, const field::Field& field )
{
  array::ArrayView<T,3> values = leveled_view<T>(field);
  array::ArrayT<T> surface_field ( array::make_shape(values.shape(0),values.shape(2) ) );
  array::ArrayView<T,2> surface(surface_field);
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<surface.shape(1); ++j )
    {
      surface(n,j) = 0.;
      for( size_t l=0; l<values.shape(1);++l )
        surface(n,j) += values(n,l,j);
    }
  }
  return checksum.execute( surface.data(), surface.stride(0) );
}
}

std::string EdgeColumns::checksum( const field::FieldSet& fieldset ) const {
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const field::Field& field=fieldset[f];
    if     ( field.datatype() == array::DataType::kind<int>() )
      md5 << checksum_3d_field<int>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<long>() )
      md5 << checksum_3d_field<long>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<float>() )
      md5 << checksum_3d_field<float>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<double>() )
      md5 << checksum_3d_field<double>(checksum(),field);
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string EdgeColumns::checksum( const field::Field& field ) const {
  field::FieldSet fieldset;
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


EdgeColumns* atlas__functionspace__Edges__new ( mesh::Mesh* mesh, int halo )
{
  EdgeColumns* edges;
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      edges = new EdgeColumns(*mesh,mesh::Halo(halo));
  );
  return edges;
}

//------------------------------------------------------------------------------

EdgeColumns* atlas__functionspace__Edges__new_mesh ( mesh::Mesh* mesh )
{
  EdgeColumns* edges;
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      edges = new EdgeColumns(*mesh);
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

mesh::Mesh* atlas__functionspace__Edges__mesh(EdgeColumns* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->mesh();
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

field::Field* atlas__functionspace__Edges__create_field (const EdgeColumns* This, const char* name, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),array::DataType(kind));
  );
  return 0;
}

//------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_field_vars (
    const EdgeColumns* This,
    const char* name,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{

  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    return This->createField(
      std::string(name),
      array::DataType(kind),
      variables_to_vector(variables,variables_size,fortran_ordering) );
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_field_lev (const EdgeColumns* This, const char* name, int levels, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),array::DataType(kind),size_t(levels));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_field_lev_vars (
    const EdgeColumns* This,
    const char* name,
    int levels,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    return This->createField(
      std::string(name),
      array::DataType(kind),
      size_t(levels),
      variables_to_vector(variables,variables_size,fortran_ordering) );
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_field_template (const EdgeColumns* This, const char* name, const field::Field* field_template )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),*field_template);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_global_field (const EdgeColumns* This, const char* name, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),array::DataType(kind));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_global_field_vars (
    const EdgeColumns* This,
    const char* name,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    return This->createGlobalField(
      std::string(name),
      array::DataType(kind),
      variables_to_vector(variables,variables_size,fortran_ordering) );
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_global_field_lev (
    const EdgeColumns* This,
    const char* name,
    int levels,
    int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),array::DataType(kind),size_t(levels));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_global_field_lev_vars (
    const EdgeColumns* This,
    const char* name,
    int levels,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    return This->createGlobalField(
      std::string(name),
      array::DataType(kind),
      size_t(levels),
      variables_to_vector(variables,variables_size,fortran_ordering) );
  );
  return 0;
}

// -----------------------------------------------------------------------------------

field::Field* atlas__functionspace__Edges__create_global_field_template (
    const EdgeColumns* This,
    const char* name,
    const field::Field* field_template )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),*field_template);
  );
  return 0;
}


// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_fieldset(
    const EdgeColumns* This,
    field::FieldSet* fieldset)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    This->haloExchange(*fieldset);
  );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_field(const EdgeColumns* This, field::Field* field)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(field);
        This->haloExchange(*field);
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
    const field::FieldSet* local,
    field::FieldSet* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->gather(*local,*global); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__gather_field(
    const EdgeColumns* This,
    const field::Field* local,
    field::Field* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->gather(*local,*global); );
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

void atlas__functionspace__Edges__scatter_fieldset(const EdgeColumns* This, const field::FieldSet* global, field::FieldSet* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->scatter(*global,*local); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__scatter_field(const EdgeColumns* This, const field::Field* global, field::Field* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(global);
        ASSERT(local);
        This->scatter(*global,*local); );
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
    const field::FieldSet* fieldset,
    char* &checksum,
    int &size,
    int &allocated)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__checksum_field(
    const EdgeColumns* This,
    const field::Field* field,
    char* &checksum,
    int &size,
    int &allocated)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

}

// -----------------------------------------------------------------------------------



} // namespace functionspace
} // namespace atlas

