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
#include "atlas/atlas_config.h"
#include "atlas/mpi/Collectives.h"
#include "atlas/Mesh.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mpl/HaloExchange.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/util/IsGhost.h"
#include "atlas/functionspace/Edges.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/atlas_omp.h"
#include "atlas/runtime/ErrorHandling.h"

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {
namespace functionspace {

namespace {

template <typename T>
ArrayView<T,3> leveled_view(const Field &field)
{
  if( field.has_levels() )
    return ArrayView<T,3> ( field.data<T>(), make_shape(field.shape(0),field.shape(1),field.stride(1)) );
  else
    return ArrayView<T,3> ( field.data<T>(), make_shape(field.shape(0),1,field.stride(0)) );
}

template <typename T>
ArrayView<T,2> surface_view(const Field &field)
{
  return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),field.stride(0)) );
}

template <typename T>
ArrayView<T,2> leveled_scalar_view(const Field &field)
{
  if( field.has_levels() )
    return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),field.shape(1)) );
  else
    return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),1) );
}

template <typename T>
ArrayView<T,1> surface_scalar_view(const Field &field)
{
  return ArrayView<T,1> ( field.data<T>(), make_shape(field.size()) );
}


}

Edges::Edges( Mesh& mesh )
  : mesh_(&mesh),
    edges_(mesh.edges()),
    nb_edges_(0),
    nb_edges_global_(0)
{
  constructor();
}

Edges::Edges( Mesh& mesh, const Halo &halo, const eckit::Parametrisation &params )
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

Edges::Edges(Mesh& mesh, const Halo &halo)
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


void Edges::constructor()
{
  nb_edges_ = mesh().edges().size();

  gather_scatter_.reset(new mpl::GatherScatter());
  halo_exchange_.reset(new mpl::HaloExchange());
  checksum_.reset(new mpl::Checksum());

  const Field& partition    = edges().partition();
  const Field& remote_index = edges().remote_index();
  const Field& global_index = edges().global_index();

  halo_exchange_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        nb_edges_);

  gather_scatter_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        edges_.global_index().data<gidx_t>(),
        -1,
        nb_edges_);

  checksum_->setup(
        partition.data<int>(),
        remote_index.data<int>(),REMOTE_IDX_BASE,
        global_index.data<gidx_t>(),
        -1,
        nb_edges_);

  nb_edges_global_ = gather_scatter_->glb_dof();
}

Edges::~Edges() {}

size_t Edges::nb_edges() const
{
  return nb_edges_;
}

size_t Edges::nb_edges_global() const
{
  return nb_edges_global_;
}

Field* Edges::createField(const std::string& name,DataType datatype) const {
  Field* field = Field::create(name,datatype,make_shape(nb_edges()));
  field->set_functionspace(this);
  return field;
}

Field* Edges::createField(const std::string& name,DataType datatype, size_t levels) const {
  Field* field = Field::create(name,datatype,make_shape(nb_edges(),levels));
  field->set_levels(levels);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createField(const std::string& name,DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = Field::create(name,datatype,shape);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createField(const std::string& name, DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = Field::create(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_edges();
  Field* field = Field::create(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(this);
  return field;
}

Field* Edges::createGlobalField(const std::string& name,DataType datatype) const {
  Field* field = Field::create(name,datatype,make_shape(nb_edges_global()));
  field->set_functionspace(this);
  return field;
}

Field* Edges::createGlobalField(const std::string& name, DataType datatype, size_t levels) const {
  Field* field = Field::create(name,datatype,make_shape(nb_edges_global(),levels));
  field->set_levels(levels);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createGlobalField(const std::string& name, DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = Field::create(name,datatype,shape);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createGlobalField(const std::string& name, DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_edges_global()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = Field::create(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(this);
  return field;
}

Field* Edges::createGlobalField(const std::string& name,const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_edges_global();
  Field* field = Field::create(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(this);
  return field;
}

void Edges::haloExchange( FieldSet& fieldset ) const
{
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field = fieldset[f];
    if     ( field.datatype() == DataType::kind<int>() ) {
      ArrayView<int,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == DataType::kind<long>() ) {
      ArrayView<long,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == DataType::kind<float>() ) {
      ArrayView<float,2> view(field);
      halo_exchange().execute( view );
    }
    else if( field.datatype() == DataType::kind<double>() ) {
      ArrayView<double,2> view(field);
      halo_exchange().execute( view );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void Edges::haloExchange( Field& field ) const
{
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
}
const mpl::HaloExchange& Edges::halo_exchange() const
{
  return *halo_exchange_;
}


void Edges::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];
    const size_t nb_fields = 1;
    if     ( loc.datatype() == DataType::kind<int>() ) {
      mpl::Field<int const> loc_field(loc.data<int>(),loc.stride(0));
      mpl::Field<int      > glb_field(glb.data<int>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<long>() ) {
      mpl::Field<long const> loc_field(loc.data<long>(),loc.stride(0));
      mpl::Field<long      > glb_field(glb.data<long>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<float>() ) {
      mpl::Field<float const> loc_field(loc.data<float>(),loc.stride(0));
      mpl::Field<float      > glb_field(glb.data<float>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<double>() ) {
      mpl::Field<double const> loc_field(loc.data<double>(),loc.stride(0));
      mpl::Field<double      > glb_field(glb.data<double>(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void Edges::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}
const mpl::GatherScatter& Edges::gather() const
{
  return *gather_scatter_;
}
const mpl::GatherScatter& Edges::scatter() const
{
  return *gather_scatter_;
}


void Edges::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    const size_t nb_fields = 1;

    if     ( loc.datatype() == DataType::kind<int>() ) {
      mpl::Field<int const> glb_field(glb.data<int>(),glb.stride(0));
      mpl::Field<int      > loc_field(loc.data<int>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<long>() ) {
      mpl::Field<long const> glb_field(glb.data<long>(),glb.stride(0));
      mpl::Field<long      > loc_field(loc.data<long>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<float>() ) {
      mpl::Field<float const> glb_field(glb.data<float>(),glb.stride(0));
      mpl::Field<float      > loc_field(loc.data<float>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else if( loc.datatype() == DataType::kind<double>() ) {
      mpl::Field<double const> glb_field(glb.data<double>(),glb.stride(0));
      mpl::Field<double      > loc_field(loc.data<double>(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void Edges::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const mpl::Checksum& checksum, const Field& field )
{
  ArrayView<T,3> values = leveled_view<T>(field);
  ArrayT<T> surface_field ( make_shape(values.shape(0),values.shape(2) ) );
  ArrayView<T,2> surface(surface_field);
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

std::string Edges::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field=fieldset[f];
    if     ( field.datatype() == DataType::kind<int>() )
      md5 << checksum_3d_field<int>(checksum(),field);
    else if( field.datatype() == DataType::kind<long>() )
      md5 << checksum_3d_field<long>(checksum(),field);
    else if( field.datatype() == DataType::kind<float>() )
      md5 << checksum_3d_field<float>(checksum(),field);
    else if( field.datatype() == DataType::kind<double>() )
      md5 << checksum_3d_field<double>(checksum(),field);
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string Edges::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

const mpl::Checksum& Edges::checksum() const
{
  return *checksum_;
}

// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
extern "C" {
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------


Edges* atlas__functionspace__Edges__new ( Mesh* mesh, int halo )
{
  Edges* edges;
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      edges = new Edges(*mesh,Halo(halo));
  );
  return edges;
}

// -----------------------------------------------------------------------------------


Edges* atlas__functionspace__Edges__new_mesh ( Mesh* mesh )
{
  Edges* edges;
  ATLAS_ERROR_HANDLING(
      ASSERT(mesh);
      edges = new Edges(*mesh);
  );
  return edges;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__delete (Edges* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete(This);
  );
}

// -----------------------------------------------------------------------------------

int atlas__functionspace__Edges__nb_edges(const Edges* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->nb_edges();
  );
  return 0;
}

// -----------------------------------------------------------------------------------

atlas::Mesh* atlas__functionspace__Edges__mesh(Edges* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->mesh();
  );
  return 0;
}

// -----------------------------------------------------------------------------------

mesh::Edges* atlas__functionspace__Edges__edges(Edges* This)
{
  ATLAS_ERROR_HANDLING(
      ASSERT(This);
      return &This->edges();
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_field (const Edges* This, const char* name, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),DataType(kind));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_field_vars (
    const Edges* This,
    const char* name,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    std::vector<size_t> variables_(variables_size);
    if( fortran_ordering )
      std::reverse_copy( variables, variables+variables_size,variables_.begin() );
    else
      variables_.assign(variables,variables+variables_size);
    return This->createField(std::string(name),DataType(kind),variables_);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_field_lev (const Edges* This, const char* name, int levels, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),DataType(kind),size_t(levels));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_field_lev_vars (
    const Edges* This,
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
    std::vector<size_t> variables_(variables_size);
    if( fortran_ordering )
      std::reverse_copy( variables, variables+variables_size,variables_.begin() );
    else
      variables_.assign(variables,variables+variables_size);
    return This->createField(std::string(name),DataType(kind),size_t(levels),variables_);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_field_template (const Edges* This, const char* name, const Field* field_template )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField(std::string(name),*field_template);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_global_field (const Edges* This, const char* name, int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),DataType(kind));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_global_field_vars (
    const Edges* This,
    const char* name,
    int variables[],
    int variables_size,
    int fortran_ordering,
    int kind)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(variables_size);
    std::vector<size_t> variables_(variables_size);
    if( fortran_ordering )
      std::reverse_copy( variables, variables+variables_size, variables_.begin() );
    else
      variables_.assign(variables,variables+variables_size);
    return This->createGlobalField(std::string(name),DataType(kind),variables_);
  );
    return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_global_field_lev (
    const Edges* This,
    const char* name,
    int levels,
    int kind )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),DataType(kind),size_t(levels));
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_global_field_lev_vars (
    const Edges* This,
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
    std::vector<size_t> variables_(variables_size);
    if( fortran_ordering )
      std::reverse_copy( variables, variables+variables_size, variables_.begin() );
    else
      variables_.assign(variables,variables+variables_size);
    return This->createGlobalField(std::string(name),DataType(kind),size_t(levels),variables_);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

Field* atlas__functionspace__Edges__create_global_field_template (
    const Edges* This,
    const char* name,
    const Field* field_template )
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField(std::string(name),*field_template);
  );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_fieldset(
    const Edges* This,
    FieldSet* fieldset)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    This->haloExchange(*fieldset);
  );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__halo_exchange_field(const Edges* This, Field* field)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(field);
        This->haloExchange(*field);
   );
}

// -----------------------------------------------------------------------------------

const mpl::HaloExchange* atlas__functionspace__Edges__get_halo_exchange(const Edges* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->halo_exchange();
  );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__gather_fieldset(
    const Edges* This,
    const FieldSet* local,
    FieldSet* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->gather(*local,*global); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__gather_field(
    const Edges* This,
    const Field* local,
    Field* global)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->gather(*local,*global); );
}

// -----------------------------------------------------------------------------------

const mpl::GatherScatter* atlas__functionspace__Edges__get_gather(const Edges* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->gather(); );
  return 0;
}

// -----------------------------------------------------------------------------------

const mpl::GatherScatter* atlas__functionspace__Edges__get_scatter(const Edges* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->scatter(); );
  return 0;
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__scatter_fieldset(const Edges* This, const FieldSet* global, FieldSet* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(local);
        ASSERT(global);
        This->scatter(*global,*local); );
}

// -----------------------------------------------------------------------------------

void atlas__functionspace__Edges__scatter_field(const Edges* This, const Field* global, Field* local)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        ASSERT(global);
        ASSERT(local);
        This->scatter(*global,*local); );
}

// -----------------------------------------------------------------------------------

const mpl::Checksum* atlas__functionspace__Edges__get_checksum(const Edges* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->checksum();
  );
  return 0;
}

// -----------------------------------------------------------------------------------


void atlas__functionspace__Edges__checksum_fieldset(
    const Edges* This,
    const FieldSet* fieldset,
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
    const Edges* This,
    const Field* field,
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

