/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/functionspace/NodesFunctionSpaceInterface.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

extern "C" {
NodesFunctionSpace* atlas__NodesFunctionSpace__new (const char* name, Mesh* mesh, int halo)
{
  ASSERT(mesh);
  return new NodesFunctionSpace(std::string(name),*mesh,Halo(halo));
}

void atlas__NodesFunctionSpace__delete (NodesFunctionSpace* This)
{
  ASSERT(This);
  delete(This);
}

int atlas__NodesFunctionSpace__nb_nodes(const NodesFunctionSpace* This)
{
  ASSERT(This);
  return This->nb_nodes();
}

Mesh* atlas__NodesFunctionSpace__mesh(NodesFunctionSpace* This)
{
  ASSERT(This);
  return &This->mesh();
}

Nodes* atlas__NodesFunctionSpace__nodes(NodesFunctionSpace* This)
{
  ASSERT(This);
  return &This->nodes();
}

Field* atlas__NodesFunctionSpace__create_field (const NodesFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createField(std::string(name),DataType(kind));
}

Field* atlas__NodesFunctionSpace__create_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size,variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createField(std::string(name),DataType(kind),variables_);
}

Field* atlas__NodesFunctionSpace__create_field_lev (const NodesFunctionSpace* This, const char* name, int levels, int kind )
{
  ASSERT(This);
  return This->createField(std::string(name),DataType(kind),size_t(levels));
}

Field* atlas__NodesFunctionSpace__create_field_lev_vars (const NodesFunctionSpace* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size,variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createField(std::string(name),DataType(kind),size_t(levels),variables_);
}

Field* atlas__NodesFunctionSpace__create_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createField(std::string(name),*field_template);
}

Field* atlas__NodesFunctionSpace__create_global_field (const NodesFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),DataType(kind));
}

Field* atlas__NodesFunctionSpace__create_global_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size, variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createGlobalField(std::string(name),DataType(kind),variables_);
}

Field* atlas__NodesFunctionSpace__create_global_field_lev (const NodesFunctionSpace* This, const char* name, int levels, int kind )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),DataType(kind),size_t(levels));
}

Field* atlas__NodesFunctionSpace__create_global_field_lev_vars (const NodesFunctionSpace* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size, variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createGlobalField(std::string(name),DataType(kind),size_t(levels),variables_);
}

Field* atlas__NodesFunctionSpace__create_global_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),*field_template);
}


void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodesFunctionSpace* This, FieldSet* fieldset)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING( This->haloExchange(*fieldset); );
}

void atlas__NodesFunctionSpace__halo_exchange_field(const NodesFunctionSpace* This, Field* field)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->haloExchange(*field); );
}

const mpl::HaloExchange* atlas__NodesFunctionSpace__get_halo_exchange(const NodesFunctionSpace* This)
{
  ASSERT(This);
  ATLAS_ERROR_HANDLING( return &This->halo_exchange(); );
  return 0;
}

void atlas__NodesFunctionSpace__gather_fieldset(const NodesFunctionSpace* This, const FieldSet* local, FieldSet* global)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->gather(*local,*global); );
}

void atlas__NodesFunctionSpace__gather_field(const NodesFunctionSpace* This, const Field* local, Field* global)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->gather(*local,*global); );
}

const mpl::GatherScatter* atlas__NodesFunctionSpace__get_gather(const NodesFunctionSpace* This)
{
  ASSERT(This);
  ATLAS_ERROR_HANDLING( return &This->gather(); );
  return 0;
}

const mpl::GatherScatter* atlas__NodesFunctionSpace__get_scatter(const NodesFunctionSpace* This)
{
  ASSERT(This);
  ATLAS_ERROR_HANDLING( return &This->scatter(); );
  return 0;
}

void atlas__NodesFunctionSpace__scatter_fieldset(const NodesFunctionSpace* This, const FieldSet* global, FieldSet* local)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->scatter(*global,*local); );
}

void atlas__NodesFunctionSpace__scatter_field(const NodesFunctionSpace* This, const Field* global, Field* local)
{
  ASSERT(This);
  ASSERT(global);
  ASSERT(local);
  ATLAS_ERROR_HANDLING( This->scatter(*global,*local); );
}

const mpl::Checksum* atlas__NodesFunctionSpace__get_checksum(const NodesFunctionSpace* This)
{
  ASSERT(This);
  ATLAS_ERROR_HANDLING( return &This->checksum(); );
  return 0;
}

void atlas__NodesFunctionSpace__checksum_fieldset(const NodesFunctionSpace* This, const FieldSet* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__NodesFunctionSpace__checksum_field(const NodesFunctionSpace* This, const Field* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__NodesFunctionSpace__sum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_long(const NodesFunctionSpace* This, const Field* field, long &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_int(const NodesFunctionSpace* This, const Field* field, int &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}


void atlas__NodesFunctionSpace__sum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new double[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<float> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new float[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_long(const NodesFunctionSpace* This, const Field* field, long* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<long> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new long[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_int(const NodesFunctionSpace* This, const Field* field, int* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<int> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new int[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}


void atlas__NodesFunctionSpace__oisum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new double[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<float> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new float[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__min_double(const NodesFunctionSpace* This, const Field* field, double &minimum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->minimum(*field,minimum) );
}

void atlas__NodesFunctionSpace__min_float(const NodesFunctionSpace* This, const Field* field, float &minimum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->minimum(*field,minimum) );
}

void atlas__NodesFunctionSpace__min_long(const NodesFunctionSpace* This, const Field* field, long &minimum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->minimum(*field,minimum) );
}

void atlas__NodesFunctionSpace__min_int(const NodesFunctionSpace* This, const Field* field, int &minimum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->minimum(*field,minimum) );
}

void atlas__NodesFunctionSpace__max_double(const NodesFunctionSpace* This, const Field* field, double &maximum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->maximum(*field,maximum) );
}

void atlas__NodesFunctionSpace__max_float(const NodesFunctionSpace* This, const Field* field, float &maximum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->maximum(*field,maximum) );
}

void atlas__NodesFunctionSpace__max_long(const NodesFunctionSpace* This, const Field* field, long &maximum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->maximum(*field,maximum) );
}

void atlas__NodesFunctionSpace__max_int(const NodesFunctionSpace* This, const Field* field, int &maximum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->maximum(*field,maximum) );
}

void atlas__NodesFunctionSpace__min_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> minvec;
        This->minimum(*field,minvec);
        size = minvec.size();
        minimum = new double[size];
        for( size_t j=0; j<size; ++j ) minimum[j] = minvec[j];
  );
}

void atlas__NodesFunctionSpace__min_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> minvec;
        This->minimum(*field,minvec);
        size = minvec.size();
        minimum = new float[size];
        for( size_t j=0; j<size; ++j ) minimum[j] = minvec[j];
  );
}

void atlas__NodesFunctionSpace__min_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> minvec;
        This->minimum(*field,minvec);
        size = minvec.size();
        minimum = new long[size];
        for( size_t j=0; j<size; ++j ) minimum[j] = minvec[j];
  );
}

void atlas__NodesFunctionSpace__min_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> minvec;
        This->minimum(*field,minvec);
        size = minvec.size();
        minimum = new int[size];
        for( size_t j=0; j<size; ++j ) minimum[j] = minvec[j];
  );
}


void atlas__NodesFunctionSpace__max_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> maxvec;
        This->maximum(*field,maxvec);
        size = maxvec.size();
        maximum = new double[size];
        for( size_t j=0; j<size; ++j ) maximum[j] = maxvec[j];
  );
}

void atlas__NodesFunctionSpace__max_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> maxvec;
        This->maximum(*field,maxvec);
        size = maxvec.size();
        maximum = new float[size];
        for( size_t j=0; j<size; ++j ) maximum[j] = maxvec[j];
  );
}

void atlas__NodesFunctionSpace__max_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> maxvec;
        This->maximum(*field,maxvec);
        size = maxvec.size();
        maximum = new long[size];
        for( size_t j=0; j<size; ++j ) maximum[j] = maxvec[j];
  );
}

void atlas__NodesFunctionSpace__max_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> maxvec;
        This->maximum(*field,maxvec);
        size = maxvec.size();
        maximum = new int[size];
        for( size_t j=0; j<size; ++j ) maximum[j] = maxvec[j];
  );
}

void atlas__NodesFunctionSpace__minloc_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_float(const NodesFunctionSpace* This, const Field* field, float &minimum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_long(const NodesFunctionSpace* This, const Field* field, long &minimum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_int(const NodesFunctionSpace* This, const Field* field, int &minimum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_float(const NodesFunctionSpace* This, const Field* field, float &maximum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_long(const NodesFunctionSpace* This, const Field* field, long &maximum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_int(const NodesFunctionSpace* This, const Field* field, int &maximum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> minvec;
        std::vector<gidx_t> gidxvec;
        This->minimumAndLocation(*field,minvec,gidxvec);
        size = minvec.size();
        minimum = new double[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> minvec;
        std::vector<gidx_t> gidxvec;
        This->minimumAndLocation(*field,minvec,gidxvec);
        size = minvec.size();
        minimum = new float[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> minvec;
        std::vector<gidx_t> gidxvec;
        This->minimumAndLocation(*field,minvec,gidxvec);
        size = minvec.size();
        minimum = new long[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> minvec;
        std::vector<gidx_t> gidxvec;
        This->minimumAndLocation(*field,minvec,gidxvec);
        size = minvec.size();
        minimum = new int[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> maxvec;
        std::vector<gidx_t> gidxvec;
        This->maximumAndLocation(*field,maxvec,gidxvec);
        size = maxvec.size();
        maximum = new double[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> maxvec;
        std::vector<gidx_t> gidxvec;
        This->maximumAndLocation(*field,maxvec,gidxvec);
        size = maxvec.size();
        maximum = new float[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> maxvec;
        std::vector<gidx_t> gidxvec;
        This->maximumAndLocation(*field,maxvec,gidxvec);
        size = maxvec.size();
        maximum = new long[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> maxvec;
        std::vector<gidx_t> gidxvec;
        This->maximumAndLocation(*field,maxvec,gidxvec);
        size = maxvec.size();
        maximum = new int[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__mean_double(const NodesFunctionSpace* This, const Field* field, double &mean, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->mean(*field,mean,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_float(const NodesFunctionSpace* This, const Field* field, float &mean, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->mean(*field,mean,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_long(const NodesFunctionSpace* This, const Field* field, long &mean, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> meanvec;
        This->mean(*field,meanvec,size_t_N);
        size = meanvec.size();
        mean = new double[size];
        for( size_t j=0; j<size; ++j ) mean[j] = meanvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<float> meanvec;
        This->mean(*field,meanvec,size_t_N);
        size = meanvec.size();
        mean = new float[size];
        for( size_t j=0; j<size; ++j ) mean[j] = meanvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, int &size, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int &size, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodesFunctionSpace* This, const Field* field, double &mean, double &stddev, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->meanAndStandardDeviation(*field,mean,stddev,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodesFunctionSpace* This, const Field* field, float &mean, float &stddev, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->meanAndStandardDeviation(*field,mean,stddev,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodesFunctionSpace* This, const Field* field, long &mean, long &stddev, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &stddev, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, double* &stddev, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> meanvec;
        std::vector<double> stddevvec;
        This->meanAndStandardDeviation(*field,meanvec,stddevvec,size_t_N);
        size = meanvec.size();
        mean = new double[size];
        stddev = new double[size];
        for( size_t j=0; j<size; ++j ) {
          mean[j] = meanvec[j];
          stddev[j] = stddevvec[j];
        }
  );
  N = size_t_N;
}


void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, float* &stddev, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<float> meanvec;
        std::vector<float> stddevvec;
        This->meanAndStandardDeviation(*field,meanvec,stddevvec,size_t_N);
        size = meanvec.size();
        mean = new float[size];
        stddev = new float[size];
        for( size_t j=0; j<size; ++j ) {
          mean[j] = meanvec[j];
          stddev[j] = stddevvec[j];
        }
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, long* &stddev, int &size, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int* &stddev, int &size, int &N)
{
NOTIMP;
}

void atlas__NodesFunctionSpace__minloclev_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__minloclev_float(const NodesFunctionSpace* This, const Field* field, float &minimum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__minloclev_long(const NodesFunctionSpace* This, const Field* field, long &minimum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__minloclev_int(const NodesFunctionSpace* This, const Field* field, int &minimum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__maxloclev_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__maxloclev_float(const NodesFunctionSpace* This, const Field* field, float &maximum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__maxloclev_long(const NodesFunctionSpace* This, const Field* field, long &maximum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__maxloclev_int(const NodesFunctionSpace* This, const Field* field, int &maximum, long &glb_idx, int &level)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  size_t lev;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx,lev) );
  glb_idx = gidx;
  level = lev;
}

void atlas__NodesFunctionSpace__minloclev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> minvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->minimumAndLocation(*field,minvec,gidxvec,levvec);
        size = minvec.size();
        minimum = new double[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloclev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> minvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->minimumAndLocation(*field,minvec,gidxvec,levvec);
        size = minvec.size();
        minimum = new float[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloclev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> minvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->minimumAndLocation(*field,minvec,gidxvec,levvec);
        size = minvec.size();
        minimum = new long[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__minloclev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> minvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->minimumAndLocation(*field,minvec,gidxvec,levvec);
        size = minvec.size();
        minimum = new int[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloclev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> maxvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->maximumAndLocation(*field,maxvec,gidxvec,levvec);
        size = maxvec.size();
        maximum = new double[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloclev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<float> maxvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->maximumAndLocation(*field,maxvec,gidxvec,levvec);
        size = maxvec.size();
        maximum = new float[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloclev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<long> maxvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->maximumAndLocation(*field,maxvec,gidxvec,levvec);
        size = maxvec.size();
        maximum = new long[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloclev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, long* &glb_idx, int* &level, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<int> maxvec;
        std::vector<gidx_t> gidxvec;
        std::vector<size_t>    levvec;
        This->maximumAndLocation(*field,maxvec,gidxvec,levvec);
        size = maxvec.size();
        maximum = new int[size];
        glb_idx = new long[size];
        level   = new int[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
          level[j] = levvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__sum_per_level(const NodesFunctionSpace* This, const Field* field, Field* sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(sum);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sumPerLevel(*field,*sum,size_t_N); );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_per_level(const NodesFunctionSpace* This, const Field* field, Field* sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(sum);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->orderIndependentSumPerLevel(*field,*sum,size_t_N); );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__min_per_level(const NodesFunctionSpace* This, const Field* field, Field* min)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(min);
  ATLAS_ERROR_HANDLING( This->minimumPerLevel(*field,*min); );
}

void atlas__NodesFunctionSpace__max_per_level(const NodesFunctionSpace* This, const Field* field, Field* max)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(max);
  ATLAS_ERROR_HANDLING( This->maximumPerLevel(*field,*max); );
}

void atlas__NodesFunctionSpace__minloc_per_level(const NodesFunctionSpace* This, const Field* field, Field* min, Field* glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(min);
  ASSERT(glb_idx);
  ATLAS_ERROR_HANDLING( This->minimumAndLocationPerLevel(*field,*min,*glb_idx); );
}

void atlas__NodesFunctionSpace__maxloc_per_level(const NodesFunctionSpace* This, const Field* field, Field* max, Field* glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(max);
  ASSERT(glb_idx);
  ATLAS_ERROR_HANDLING( This->maximumAndLocationPerLevel(*field,*max,*glb_idx); );
}

void atlas__NodesFunctionSpace__mean_per_level(const NodesFunctionSpace* This, const Field* field, Field* mean, int &N)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(mean);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->meanPerLevel(*field,*mean,size_t_N); );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_per_level(const NodesFunctionSpace* This, const Field* field, Field* mean, Field* stddev, int &N)
{
  ASSERT(This);
  ASSERT(field);
  ASSERT(mean);
  ASSERT(stddev);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->meanAndStandardDeviationPerLevel(*field,*mean,*stddev,size_t_N); );
  N = size_t_N;
}


}


} // namespace functionspace
} // namespace atlas

