/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace field {
class FieldSetImpl;
class FieldImpl;
}  // namespace field
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {

extern "C" {
const NodeColumns* atlas__NodesFunctionSpace__new(Mesh::Implementation* mesh, const eckit::Configuration* config);
void atlas__NodesFunctionSpace__delete(NodeColumns* This);
int atlas__NodesFunctionSpace__nb_nodes(const NodeColumns* This);
const Mesh::Implementation* atlas__NodesFunctionSpace__mesh(const NodeColumns* This);
mesh::Nodes* atlas__NodesFunctionSpace__nodes(const NodeColumns* This);

void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodeColumns* This, field::FieldSetImpl* fieldset);
void atlas__NodesFunctionSpace__halo_exchange_field(const NodeColumns* This, field::FieldImpl* field);
const parallel::HaloExchange* atlas__NodesFunctionSpace__get_halo_exchange(const NodeColumns* This);

void atlas__NodesFunctionSpace__gather_fieldset(const NodeColumns* This, const field::FieldSetImpl* local,
                                                field::FieldSetImpl* global);
void atlas__NodesFunctionSpace__gather_field(const NodeColumns* This, const field::FieldImpl* local,
                                             field::FieldImpl* global);
const parallel::GatherScatter* atlas__NodesFunctionSpace__get_gather(const NodeColumns* This);

void atlas__NodesFunctionSpace__scatter_fieldset(const NodeColumns* This, const field::FieldSetImpl* global,
                                                 field::FieldSetImpl* local);
void atlas__NodesFunctionSpace__scatter_field(const NodeColumns* This, const field::FieldImpl* global,
                                              field::FieldImpl* local);
const parallel::GatherScatter* atlas__NodesFunctionSpace__get_scatter(const NodeColumns* This);

void atlas__NodesFunctionSpace__checksum_fieldset(const NodeColumns* This, const field::FieldSetImpl* fieldset,
                                                  char*& checksum, int& size, int& allocated);
void atlas__NodesFunctionSpace__checksum_field(const NodeColumns* This, const field::FieldImpl* field, char*& checksum,
                                               int& size, int& allocated);
const parallel::Checksum* atlas__NodesFunctionSpace__get_checksum(const NodeColumns* This);

void atlas__NodesFunctionSpace__sum_double(const NodeColumns* This, const field::FieldImpl* field, double& sum, int& N);
void atlas__NodesFunctionSpace__sum_float(const NodeColumns* This, const field::FieldImpl* field, float& sum, int& N);
void atlas__NodesFunctionSpace__sum_int(const NodeColumns* This, const field::FieldImpl* field, int& sum, int& N);
void atlas__NodesFunctionSpace__sum_long(const NodeColumns* This, const field::FieldImpl* field, long& sum, int& N);
void atlas__NodesFunctionSpace__sum_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                               int& size, int& N);
void atlas__NodesFunctionSpace__sum_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                              int& size, int& N);
void atlas__NodesFunctionSpace__sum_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& sum,
                                            int& size, int& N);
void atlas__NodesFunctionSpace__sum_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& sum,
                                             int& size, int& N);

void atlas__NodesFunctionSpace__oisum_double(const NodeColumns* This, const field::FieldImpl* field, double& sum,
                                             int& N);
void atlas__NodesFunctionSpace__oisum_float(const NodeColumns* This, const field::FieldImpl* field, float& sum, int& N);
void atlas__NodesFunctionSpace__oisum_int(const NodeColumns* This, const field::FieldImpl* field, int& sum, int& N);
void atlas__NodesFunctionSpace__oisum_long(const NodeColumns* This, const field::FieldImpl* field, long& sum, int& N);
void atlas__NodesFunctionSpace__oisum_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& sum,
                                                 int& size, int& N);
void atlas__NodesFunctionSpace__oisum_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& sum,
                                                int& size, int& N);
void atlas__NodesFunctionSpace__oisum_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& sum,
                                              int& size, int& N);
void atlas__NodesFunctionSpace__oisum_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& sum,
                                               int& size, int& N);

void atlas__NodesFunctionSpace__min_double(const NodeColumns* This, const field::FieldImpl* field, double& minimum);
void atlas__NodesFunctionSpace__min_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum);
void atlas__NodesFunctionSpace__min_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum);
void atlas__NodesFunctionSpace__min_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum);

void atlas__NodesFunctionSpace__max_double(const NodeColumns* This, const field::FieldImpl* field, double& maximum);
void atlas__NodesFunctionSpace__max_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum);
void atlas__NodesFunctionSpace__max_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum);
void atlas__NodesFunctionSpace__max_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum);

void atlas__NodesFunctionSpace__min_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& minimum,
                                               int& size);
void atlas__NodesFunctionSpace__min_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& minimum,
                                              int& size);
void atlas__NodesFunctionSpace__min_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                            int& size);
void atlas__NodesFunctionSpace__min_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                             int& size);

void atlas__NodesFunctionSpace__max_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& maximum,
                                               int& size);
void atlas__NodesFunctionSpace__max_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& maximum,
                                              int& size);
void atlas__NodesFunctionSpace__max_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                            int& size);
void atlas__NodesFunctionSpace__max_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                             int& size);

void atlas__NodesFunctionSpace__minloc_double(const NodeColumns* This, const field::FieldImpl* field, double& minimum,
                                              long& glb_idx);
void atlas__NodesFunctionSpace__minloc_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                             long& glb_idx);
void atlas__NodesFunctionSpace__minloc_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                           long& glb_idx);
void atlas__NodesFunctionSpace__minloc_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                            long& glb_idx);

void atlas__NodesFunctionSpace__maxloc_double(const NodeColumns* This, const field::FieldImpl* field, double& maximum,
                                              long& glb_idx);
void atlas__NodesFunctionSpace__maxloc_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                             long& glb_idx);
void atlas__NodesFunctionSpace__maxloc_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                           long& glb_idx);
void atlas__NodesFunctionSpace__maxloc_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                            long& glb_idx);

void atlas__NodesFunctionSpace__minloc_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                  double*& minimum, long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__minloc_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                 float*& minimum, long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__minloc_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                               long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__minloc_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& minimum,
                                                long*& glb_idx, int& size);

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                  double*& maximum, long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__maxloc_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                 float*& maximum, long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__maxloc_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                               long*& glb_idx, int& size);
void atlas__NodesFunctionSpace__maxloc_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& maximum,
                                                long*& glb_idx, int& size);

void atlas__NodesFunctionSpace__mean_double(const NodeColumns* This, const field::FieldImpl* field, double& mean,
                                            int& N);
void atlas__NodesFunctionSpace__mean_float(const NodeColumns* This, const field::FieldImpl* field, float& mean, int& N);
void atlas__NodesFunctionSpace__mean_int(const NodeColumns* This, const field::FieldImpl* field, int& mean, int& N);
void atlas__NodesFunctionSpace__mean_long(const NodeColumns* This, const field::FieldImpl* field, long& mean, int& N);
void atlas__NodesFunctionSpace__mean_arr_double(const NodeColumns* This, const field::FieldImpl* field, double*& mean,
                                                int& size, int& N);
void atlas__NodesFunctionSpace__mean_arr_float(const NodeColumns* This, const field::FieldImpl* field, float*& mean,
                                               int& size, int& N);
void atlas__NodesFunctionSpace__mean_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& mean,
                                             int& size, int& N);
void atlas__NodesFunctionSpace__mean_arr_long(const NodeColumns* This, const field::FieldImpl* field, long*& mean,
                                              int& size, int& N);

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                       double& mean, double& stddev, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodeColumns* This, const field::FieldImpl* field,
                                                      float& mean, float& stddev, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodeColumns* This, const field::FieldImpl* field, int& mean,
                                                    int& stddev, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodeColumns* This, const field::FieldImpl* field, long& mean,
                                                     long& stddev, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                           double*& mean, double*& stddev, int& size, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                          float*& mean, float*& stddev, int& size, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodeColumns* This, const field::FieldImpl* field,
                                                        int*& mean, int*& stddev, int& size, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                         long*& mean, long*& stddev, int& size, int& N);

void atlas__NodesFunctionSpace__minloclev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                 double& minimum, long& glb_idx, int& level);
void atlas__NodesFunctionSpace__minloclev_float(const NodeColumns* This, const field::FieldImpl* field, float& minimum,
                                                long& glb_idx, int& level);
void atlas__NodesFunctionSpace__minloclev_int(const NodeColumns* This, const field::FieldImpl* field, int& minimum,
                                              long& glb_idx, int& level);
void atlas__NodesFunctionSpace__minloclev_long(const NodeColumns* This, const field::FieldImpl* field, long& minimum,
                                               long& glb_idx, int& level);

void atlas__NodesFunctionSpace__maxloclev_double(const NodeColumns* This, const field::FieldImpl* field,
                                                 double& maximum, long& glb_idx, int& level);
void atlas__NodesFunctionSpace__maxloclev_float(const NodeColumns* This, const field::FieldImpl* field, float& maximum,
                                                long& glb_idx, int& level);
void atlas__NodesFunctionSpace__maxloclev_int(const NodeColumns* This, const field::FieldImpl* field, int& maximum,
                                              long& glb_idx, int& level);
void atlas__NodesFunctionSpace__maxloclev_long(const NodeColumns* This, const field::FieldImpl* field, long& maximum,
                                               long& glb_idx, int& level);

void atlas__NodesFunctionSpace__minloclev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                     double*& minimum, long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__minloclev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                    float*& minimum, long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__minloclev_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& minimum,
                                                  long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__minloclev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                   long*& minimum, long*& glb_idx, int*& level, int& size);

void atlas__NodesFunctionSpace__maxloclev_arr_double(const NodeColumns* This, const field::FieldImpl* field,
                                                     double*& maximum, long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__maxloclev_arr_float(const NodeColumns* This, const field::FieldImpl* field,
                                                    float*& maximum, long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__maxloclev_arr_int(const NodeColumns* This, const field::FieldImpl* field, int*& maximum,
                                                  long*& glb_idx, int*& level, int& size);
void atlas__NodesFunctionSpace__maxloclev_arr_long(const NodeColumns* This, const field::FieldImpl* field,
                                                   long*& maximum, long*& glb_idx, int*& level, int& size);

void atlas__NodesFunctionSpace__sum_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* sum, int& N);
void atlas__NodesFunctionSpace__oisum_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                field::FieldImpl* sum, int& N);
void atlas__NodesFunctionSpace__min_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* min);
void atlas__NodesFunctionSpace__max_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                              field::FieldImpl* max);
void atlas__NodesFunctionSpace__minloc_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* min, field::FieldImpl* glb_idx);
void atlas__NodesFunctionSpace__maxloc_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                 field::FieldImpl* max, field::FieldImpl* glb_idx);
void atlas__NodesFunctionSpace__mean_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                               field::FieldImpl* mean, int& N);
void atlas__NodesFunctionSpace__mean_and_stddev_per_level(const NodeColumns* This, const field::FieldImpl* field,
                                                          field::FieldImpl* mean, field::FieldImpl* stddev, int& N);
}

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
