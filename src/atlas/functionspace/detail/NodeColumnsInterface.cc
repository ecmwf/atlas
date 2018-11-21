/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "NodeColumnsInterface.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {
namespace detail {

using atlas::FieldSet;
using atlas::field::FieldImpl;
using atlas::field::FieldSetImpl;

// ----------------------------------------------------------------------

extern "C" {
const NodeColumns* atlas__NodesFunctionSpace__new( Mesh::Implementation* mesh, const eckit::Configuration* config ) {
    ASSERT( mesh );
    Mesh m( mesh );
    return new NodeColumns( m, *config );
}

void atlas__NodesFunctionSpace__delete( NodeColumns* This ) {
    ASSERT( This );
    delete ( This );
}

int atlas__NodesFunctionSpace__nb_nodes( const NodeColumns* This ) {
    ASSERT( This );
    return This->nb_nodes();
}

const Mesh::Implementation* atlas__NodesFunctionSpace__mesh( const NodeColumns* This ) {
    ASSERT( This );
    return This->mesh().get();
}

mesh::Nodes* atlas__NodesFunctionSpace__nodes( const NodeColumns* This ) {
    ASSERT( This );
    return &This->nodes();
}

field::FieldImpl* atlas__NodesFunctionSpace__create_field( const NodeColumns* This,
                                                           const eckit::Configuration* options ) {
    ASSERT( This );
    ASSERT( options );
    FieldImpl* field;
    {
        Field f = This->createField( util::Config( *options ) );
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

field::FieldImpl* atlas__NodesFunctionSpace__create_field_template( const NodeColumns* This,
                                                                    const field::FieldImpl* field_template,
                                                                    const eckit::Configuration* options ) {
    ASSERT( This );
    ASSERT( options );
    FieldImpl* field;
    {
        Field f = This->createField( Field( field_template ), *options );
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

void atlas__NodesFunctionSpace__halo_exchange_fieldset( const NodeColumns* This, field::FieldSetImpl* fieldset ) {
    ASSERT( This );
    ASSERT( fieldset );
    FieldSet f( fieldset );
    ATLAS_ERROR_HANDLING( This->haloExchange( f ); );
}

void atlas__NodesFunctionSpace__halo_exchange_field( const NodeColumns* This, field::FieldImpl* field ) {
    ASSERT( This );
    ASSERT( field );
    Field f( field );
    ATLAS_ERROR_HANDLING( This->haloExchange( f ); );
}

const parallel::HaloExchange* atlas__NodesFunctionSpace__get_halo_exchange( const NodeColumns* This ) {
    ASSERT( This );
    ATLAS_ERROR_HANDLING( return &This->halo_exchange(); );
    return 0;
}

void atlas__NodesFunctionSpace__gather_fieldset( const NodeColumns* This, const field::FieldSetImpl* local,
                                                 field::FieldSetImpl* global ) {
    ASSERT( This );
    ASSERT( local );
    ASSERT( global );
    const FieldSet l( local );
    FieldSet g( global );
    ATLAS_ERROR_HANDLING( This->gather( l, g ); );
}

void atlas__NodesFunctionSpace__gather_field( const NodeColumns* This, const field::FieldImpl* local,
                                              field::FieldImpl* global ) {
    ASSERT( This );
    ASSERT( local );
    ASSERT( global );
    const Field l( local );
    Field g( global );
    ATLAS_ERROR_HANDLING( This->gather( l, g ); );
}

const parallel::GatherScatter* atlas__NodesFunctionSpace__get_gather( const NodeColumns* This ) {
    ASSERT( This );
    ATLAS_ERROR_HANDLING( return &This->gather(); );
    return 0;
}

const parallel::GatherScatter* atlas__NodesFunctionSpace__get_scatter( const NodeColumns* This ) {
    ASSERT( This );
    ATLAS_ERROR_HANDLING( return &This->scatter(); );
    return 0;
}

void atlas__NodesFunctionSpace__scatter_fieldset( const NodeColumns* This, const field::FieldSetImpl* global,
                                                  field::FieldSetImpl* local ) {
    ASSERT( This );
    ASSERT( local );
    ASSERT( global );
    const FieldSet g( global );
    FieldSet l( local );
    ATLAS_ERROR_HANDLING( This->scatter( g, l ); );
}

void atlas__NodesFunctionSpace__scatter_field( const NodeColumns* This, const field::FieldImpl* global,
                                               field::FieldImpl* local ) {
    ASSERT( This );
    ASSERT( global );
    ASSERT( local );
    const Field g( global );
    Field l( local );
    ATLAS_ERROR_HANDLING( This->scatter( g, l ); );
}

const parallel::Checksum* atlas__NodesFunctionSpace__get_checksum( const NodeColumns* This ) {
    ASSERT( This );
    ATLAS_ERROR_HANDLING( return &This->checksum(); );
    return 0;
}

void atlas__NodesFunctionSpace__checksum_fieldset( const NodeColumns* This, const field::FieldSetImpl* fieldset,
                                                   char*& checksum, int& size, int& allocated ) {
    ASSERT( This );
    ASSERT( fieldset );
    ATLAS_ERROR_HANDLING( std::string checksum_str( This->checksum( fieldset ) ); size = checksum_str.size();
                          checksum = new char[size + 1]; allocated = true; strcpy( checksum, checksum_str.c_str() ); );
}

void atlas__NodesFunctionSpace__checksum_field( const NodeColumns* This, const field::FieldImpl* field, char*& checksum,
                                                int& size, int& allocated ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::string checksum_str( This->checksum( field ) ); size = checksum_str.size();
                          checksum = new char[size + 1]; allocated = true; strcpy( checksum, checksum_str.c_str() ); );
}

void atlas__NodesFunctionSpace__sum_double( const NodeColumns* This, const field::FieldImpl* field, double& sum,
                                            int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->sum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_float( const NodeColumns* This, const field::FieldImpl* field, float& sum,
                                           int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->sum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_long( const NodeColumns* This, const field::FieldImpl* field, long& sum, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->sum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_int( const NodeColumns* This, const field::FieldImpl* field, int& sum, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->sum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_double( const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                                int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<double> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new double[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_float( const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                               int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<float> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new float[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& sum,
                                              int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<long> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& sum,
                                             int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<int> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_double( const NodeColumns* This, const field::FieldImpl* field, double& sum,
                                              int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->orderIndependentSum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_float( const NodeColumns* This, const field::FieldImpl* field, float& sum,
                                             int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->orderIndependentSum( field, sum, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_arr_double( const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                                  int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<double> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new double[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_arr_float( const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                                 int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<float> sumvec; This->orderIndependentSum( field, sumvec, idx_t_N );
                          size = sumvec.size(); sum                        = new float[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) sum[j] = sumvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__min_double( const NodeColumns* This, const field::FieldImpl* field, double& minimum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->minimum( field, minimum ) );
}

void atlas__NodesFunctionSpace__min_float( const NodeColumns* This, const field::FieldImpl* field, float& minimum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->minimum( field, minimum ) );
}

void atlas__NodesFunctionSpace__min_long( const NodeColumns* This, const field::FieldImpl* field, long& minimum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->minimum( field, minimum ) );
}

void atlas__NodesFunctionSpace__min_int( const NodeColumns* This, const field::FieldImpl* field, int& minimum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->minimum( field, minimum ) );
}

void atlas__NodesFunctionSpace__max_double( const NodeColumns* This, const field::FieldImpl* field, double& maximum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->maximum( field, maximum ) );
}

void atlas__NodesFunctionSpace__max_float( const NodeColumns* This, const field::FieldImpl* field, float& maximum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->maximum( field, maximum ) );
}

void atlas__NodesFunctionSpace__max_long( const NodeColumns* This, const field::FieldImpl* field, long& maximum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->maximum( field, maximum ) );
}

void atlas__NodesFunctionSpace__max_int( const NodeColumns* This, const field::FieldImpl* field, int& maximum ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( This->maximum( field, maximum ) );
}

void atlas__NodesFunctionSpace__min_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                double*& minimum, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> minvec; This->minimum( field, minvec ); size = minvec.size();
                          minimum                                                          = new double[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) minimum[j]             = minvec[j]; );
}

void atlas__NodesFunctionSpace__min_arr_float( const NodeColumns* This, const field::FieldImpl* field, float*& minimum,
                                               int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> minvec; This->minimum( field, minvec ); size = minvec.size();
                          minimum                                                         = new float[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) minimum[j]            = minvec[j]; );
}

void atlas__NodesFunctionSpace__min_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                              int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> minvec; This->minimum( field, minvec ); size = minvec.size();
                          minimum = new long[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) minimum[j] = minvec[j]; );
}

void atlas__NodesFunctionSpace__min_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                             int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> minvec; This->minimum( field, minvec ); size = minvec.size();
                          minimum = new int[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) minimum[j] = minvec[j]; );
}

void atlas__NodesFunctionSpace__max_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                double*& maximum, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> maxvec; This->maximum( field, maxvec ); size = maxvec.size();
                          maximum                                                          = new double[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) maximum[j]             = maxvec[j]; );
}

void atlas__NodesFunctionSpace__max_arr_float( const NodeColumns* This, const field::FieldImpl* field, float*& maximum,
                                               int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> maxvec; This->maximum( field, maxvec ); size = maxvec.size();
                          maximum                                                         = new float[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) maximum[j]            = maxvec[j]; );
}

void atlas__NodesFunctionSpace__max_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                              int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> maxvec; This->maximum( field, maxvec ); size = maxvec.size();
                          maximum = new long[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) maximum[j] = maxvec[j]; );
}

void atlas__NodesFunctionSpace__max_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                             int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> maxvec; This->maximum( field, maxvec ); size = maxvec.size();
                          maximum = new int[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) maximum[j] = maxvec[j]; );
}

void atlas__NodesFunctionSpace__minloc_double( const NodeColumns* This, const field::FieldImpl* field, double& minimum,
                                               long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_float( const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                              long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_long( const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                             long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_int( const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                            long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_double( const NodeColumns* This, const field::FieldImpl* field, double& maximum,
                                               long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_float( const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                              long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_long( const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                             long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_int( const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                            long& glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx ) );
    glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                   double*& minimum, long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> minvec; std::vector<gidx_t> gidxvec;
                          This->minimumAndLocation( field, minvec, gidxvec ); size = minvec.size();
                          minimum = new double[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloc_arr_float( const NodeColumns* This, const field::FieldImpl* field,
                                                  float*& minimum, long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> minvec; std::vector<gidx_t> gidxvec;
                          This->minimumAndLocation( field, minvec, gidxvec ); size = minvec.size();
                          minimum = new float[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloc_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                                 long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> minvec; std::vector<gidx_t> gidxvec;
                          This->minimumAndLocation( field, minvec, gidxvec ); size = minvec.size();
                          minimum = new long[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloc_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                                long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> minvec; std::vector<gidx_t> gidxvec;
                          This->minimumAndLocation( field, minvec, gidxvec ); size = minvec.size();
                          minimum = new int[size]; glb_idx = new long[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloc_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                   double*& maximum, long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> maxvec; std::vector<gidx_t> gidxvec;
                          This->maximumAndLocation( field, maxvec, gidxvec ); size = maxvec.size();
                          maximum = new double[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloc_arr_float( const NodeColumns* This, const field::FieldImpl* field,
                                                  float*& maximum, long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> maxvec; std::vector<gidx_t> gidxvec;
                          This->maximumAndLocation( field, maxvec, gidxvec ); size = maxvec.size();
                          maximum = new float[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloc_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                                 long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> maxvec; std::vector<gidx_t> gidxvec;
                          This->maximumAndLocation( field, maxvec, gidxvec ); size = maxvec.size();
                          maximum = new long[size]; glb_idx = new long[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloc_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                                long*& glb_idx, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> maxvec; std::vector<gidx_t> gidxvec;
                          This->maximumAndLocation( field, maxvec, gidxvec ); size = maxvec.size();
                          maximum = new int[size]; glb_idx = new long[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                          } );
}

void atlas__NodesFunctionSpace__mean_double( const NodeColumns* This, const field::FieldImpl* field, double& mean,
                                             int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->mean( field, mean, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_float( const NodeColumns* This, const field::FieldImpl* field, float& mean,
                                            int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->mean( field, mean, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_long( const NodeColumns* This, const field::FieldImpl* field, long& mean,
                                           int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_int( const NodeColumns* This, const field::FieldImpl* field, int& mean, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_arr_double( const NodeColumns* This, const field::FieldImpl* field, double*& mean,
                                                 int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<double> meanvec; This->mean( field, meanvec, idx_t_N ); size = meanvec.size();
                          mean = new double[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) mean[j] = meanvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_arr_float( const NodeColumns* This, const field::FieldImpl* field, float*& mean,
                                                int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<float> meanvec; This->mean( field, meanvec, idx_t_N ); size = meanvec.size();
                          mean = new float[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) mean[j] = meanvec[j]; );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_arr_long( const NodeColumns* This, const field::FieldImpl* field, long*& mean,
                                               int& size, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_arr_int( const NodeColumns* This, const field::FieldImpl* field, int*& mean,
                                              int& size, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_double( const NodeColumns* This, const field::FieldImpl* field,
                                                        double& mean, double& stddev, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->meanAndStandardDeviation( field, mean, stddev, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_float( const NodeColumns* This, const field::FieldImpl* field,
                                                       float& mean, float& stddev, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( This->meanAndStandardDeviation( field, mean, stddev, idx_t_N ) );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_long( const NodeColumns* This, const field::FieldImpl* field,
                                                      long& mean, long& stddev, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_int( const NodeColumns* This, const field::FieldImpl* field, int& mean,
                                                     int& stddev, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                            double*& mean, double*& stddev, int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<double> meanvec; std::vector<double> stddevvec;
                          This->meanAndStandardDeviation( field, meanvec, stddevvec, idx_t_N ); size = meanvec.size();
                          mean = new double[size]; stddev = new double[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              mean[j]   = meanvec[j];
                              stddev[j] = stddevvec[j];
                          } );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_float( const NodeColumns* This, const field::FieldImpl* field,
                                                           float*& mean, float*& stddev, int& size, int& N ) {
    ASSERT( This );
    ASSERT( field );
    idx_t idx_t_N;
    ATLAS_ERROR_HANDLING( std::vector<float> meanvec; std::vector<float> stddevvec;
                          This->meanAndStandardDeviation( field, meanvec, stddevvec, idx_t_N ); size = meanvec.size();
                          mean = new float[size]; stddev = new float[size]; for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              mean[j]   = meanvec[j];
                              stddev[j] = stddevvec[j];
                          } );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_long( const NodeColumns* This, const field::FieldImpl* field,
                                                          long*& mean, long*& stddev, int& size, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_int( const NodeColumns* This, const field::FieldImpl* field,
                                                         int*& mean, int*& stddev, int& size, int& N ) {
    NOTIMP;
}

void atlas__NodesFunctionSpace__minloclev_double( const NodeColumns* This, const field::FieldImpl* field,
                                                  double& minimum, long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_float( const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                                 long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_long( const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                                long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_int( const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                               long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->minimumAndLocation( field, minimum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_double( const NodeColumns* This, const field::FieldImpl* field,
                                                  double& maximum, long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_float( const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                                 long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_long( const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                                long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__maxloclev_int( const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                               long& glb_idx, int& level ) {
    ASSERT( This );
    ASSERT( field );
    gidx_t gidx;
    idx_t lev;
    ATLAS_ERROR_HANDLING( This->maximumAndLocation( field, maximum, gidx, lev ) );
    glb_idx = gidx;
    level   = lev;
}

void atlas__NodesFunctionSpace__minloclev_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                      double*& minimum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> minvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->minimumAndLocation( field, minvec, gidxvec, levvec ); size = minvec.size();
                          minimum = new double[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloclev_arr_float( const NodeColumns* This, const field::FieldImpl* field,
                                                     float*& minimum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> minvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->minimumAndLocation( field, minvec, gidxvec, levvec ); size = minvec.size();
                          minimum = new float[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloclev_arr_long( const NodeColumns* This, const field::FieldImpl* field,
                                                    long*& minimum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> minvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->minimumAndLocation( field, minvec, gidxvec, levvec ); size = minvec.size();
                          minimum = new long[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__minloclev_arr_int( const NodeColumns* This, const field::FieldImpl* field,
                                                   int*& minimum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> minvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->minimumAndLocation( field, minvec, gidxvec, levvec ); size = minvec.size();
                          minimum = new int[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              minimum[j] = minvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloclev_arr_double( const NodeColumns* This, const field::FieldImpl* field,
                                                      double*& maximum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<double> maxvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->maximumAndLocation( field, maxvec, gidxvec, levvec ); size = maxvec.size();
                          maximum = new double[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloclev_arr_float( const NodeColumns* This, const field::FieldImpl* field,
                                                     float*& maximum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<float> maxvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->maximumAndLocation( field, maxvec, gidxvec, levvec ); size = maxvec.size();
                          maximum = new float[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloclev_arr_long( const NodeColumns* This, const field::FieldImpl* field,
                                                    long*& maximum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<long> maxvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->maximumAndLocation( field, maxvec, gidxvec, levvec ); size = maxvec.size();
                          maximum = new long[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__maxloclev_arr_int( const NodeColumns* This, const field::FieldImpl* field,
                                                   int*& maximum, long*& glb_idx, int*& level, int& size ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::vector<int> maxvec; std::vector<gidx_t> gidxvec; std::vector<idx_t> levvec;
                          This->maximumAndLocation( field, maxvec, gidxvec, levvec ); size = maxvec.size();
                          maximum = new int[size]; glb_idx = new long[size]; level = new int[size];
                          for ( idx_t j = 0; j < (idx_t)size; ++j ) {
                              maximum[j] = maxvec[j];
                              glb_idx[j] = gidxvec[j];
                              level[j]   = levvec[j];
                          } );
}

void atlas__NodesFunctionSpace__sum_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* column, int& N ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( column );
    idx_t idx_t_N;
    Field sum( column );
    ATLAS_ERROR_HANDLING( This->sumPerLevel( field, sum, idx_t_N ); );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__oisum_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* column, int& N ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( column );
    idx_t idx_t_N;
    Field sum( column );
    ATLAS_ERROR_HANDLING( This->orderIndependentSumPerLevel( field, sum, idx_t_N ); );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__min_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* min ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( min );
    Field fmin( min );
    ATLAS_ERROR_HANDLING( This->minimumPerLevel( field, fmin ); );
}

void atlas__NodesFunctionSpace__max_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* max ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( max );
    Field fmax( max );
    ATLAS_ERROR_HANDLING( This->maximumPerLevel( field, fmax ); );
}

void atlas__NodesFunctionSpace__minloc_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                                  field::FieldImpl* min, field::FieldImpl* glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( min );
    ASSERT( glb_idx );
    Field fmin( min );
    Field fglb_idx( glb_idx );
    ATLAS_ERROR_HANDLING( This->minimumAndLocationPerLevel( field, fmin, fglb_idx ); );
}

void atlas__NodesFunctionSpace__maxloc_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                                  field::FieldImpl* max, field::FieldImpl* glb_idx ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( max );
    ASSERT( glb_idx );
    Field fmax( max );
    Field fglb_idx( glb_idx );
    ATLAS_ERROR_HANDLING( This->maximumAndLocationPerLevel( field, fmax, fglb_idx ); );
}

void atlas__NodesFunctionSpace__mean_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                                field::FieldImpl* mean, int& N ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( mean );
    idx_t idx_t_N;
    Field fmean( mean );
    ATLAS_ERROR_HANDLING( This->meanPerLevel( field, fmean, idx_t_N ); );
    N = idx_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_per_level( const NodeColumns* This, const field::FieldImpl* field,
                                                           field::FieldImpl* mean, field::FieldImpl* stddev, int& N ) {
    ASSERT( This );
    ASSERT( field );
    ASSERT( mean );
    ASSERT( stddev );
    idx_t idx_t_N;
    Field fmean( mean );
    Field fstddev( stddev );
    ATLAS_ERROR_HANDLING( This->meanAndStandardDeviationPerLevel( field, fmean, fstddev, idx_t_N ); );
    N = idx_t_N;
}
}

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
