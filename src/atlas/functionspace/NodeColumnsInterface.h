/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_NodeColumnsFunctionSpaceInterface_h
#define atlas_functionspace_NodeColumnsFunctionSpaceInterface_h


#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace functionspace {

extern "C" {
NodeColumns* atlas__NodesFunctionSpace__new (mesh::Mesh::mesh_t* mesh, int halo);
NodeColumns* atlas__NodesFunctionSpace__new_mesh (mesh::Mesh::mesh_t* mesh);
void atlas__NodesFunctionSpace__delete (NodeColumns* This);
int atlas__NodesFunctionSpace__nb_nodes(const NodeColumns* This);
mesh::Mesh::mesh_t* atlas__NodesFunctionSpace__mesh(NodeColumns* This);
mesh::Nodes* atlas__NodesFunctionSpace__nodes(NodeColumns* This);
field::Field* atlas__NodesFunctionSpace__create_field (const NodeColumns* This, const char* name, int kind, const eckit::Parametrisation* options);
field::Field* atlas__NodesFunctionSpace__create_field_vars (const NodeColumns* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind, const eckit::Parametrisation* options);

field::Field* atlas__NodesFunctionSpace__create_field_lev (const NodeColumns* This, const char* name, int levels, int kind, const eckit::Parametrisation* options);
field::Field* atlas__NodesFunctionSpace__create_field_lev_vars (const NodeColumns* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind, const eckit::Parametrisation* options);


field::Field* atlas__NodesFunctionSpace__create_field_template (const NodeColumns* This, const char* name, const field::Field* field_template, const eckit::Parametrisation* options);

void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodeColumns* This, field::FieldSet* fieldset);
void atlas__NodesFunctionSpace__halo_exchange_field(const NodeColumns* This, field::Field* field);
const parallel::HaloExchange* atlas__NodesFunctionSpace__get_halo_exchange(const NodeColumns* This);

void atlas__NodesFunctionSpace__gather_fieldset(const NodeColumns* This, const field::FieldSet* local, field::FieldSet* global);
void atlas__NodesFunctionSpace__gather_field(const NodeColumns* This, const field::Field* local, field::Field* global);
const parallel::GatherScatter* atlas__NodesFunctionSpace__get_gather(const NodeColumns* This);

void atlas__NodesFunctionSpace__scatter_fieldset(const NodeColumns* This, const field::FieldSet* global, field::FieldSet* local);
void atlas__NodesFunctionSpace__scatter_field(const NodeColumns* This, const field::Field* global, field::Field* local);
const parallel::GatherScatter* atlas__NodesFunctionSpace__get_scatter(const NodeColumns* This);

void atlas__NodesFunctionSpace__checksum_fieldset(const NodeColumns* This, const field::FieldSet* fieldset, char* &checksum, int &size, int &allocated);
void atlas__NodesFunctionSpace__checksum_field(const NodeColumns* This, const field::Field* field, char* &checksum, int &size, int &allocated);
const parallel::Checksum* atlas__NodesFunctionSpace__get_checksum(const NodeColumns* This);

void atlas__NodesFunctionSpace__sum_double(const NodeColumns* This, const field::Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__sum_float(const NodeColumns* This, const field::Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__sum_int(const NodeColumns* This, const field::Field* field, int &sum, int &N);
void atlas__NodesFunctionSpace__sum_long(const NodeColumns* This, const field::Field* field, long &sum, int &N);
void atlas__NodesFunctionSpace__sum_arr_double(const NodeColumns* This, const field::Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_float(const NodeColumns* This, const field::Field* field, float* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_int(const NodeColumns* This, const field::Field* field, int* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_long(const NodeColumns* This, const field::Field* field, long* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__oisum_double(const NodeColumns* This, const field::Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__oisum_float(const NodeColumns* This, const field::Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__oisum_arr_double(const NodeColumns* This, const field::Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__oisum_arr_float(const NodeColumns* This, const field::Field* field, float* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__min_double(const NodeColumns* This, const field::Field* field, double &minimum);
void atlas__NodesFunctionSpace__min_float(const NodeColumns* This, const field::Field* field, float &minimum);
void atlas__NodesFunctionSpace__min_int(const NodeColumns* This, const field::Field* field, int &minimum);
void atlas__NodesFunctionSpace__min_long(const NodeColumns* This, const field::Field* field, long &minimum);

void atlas__NodesFunctionSpace__max_double(const NodeColumns* This, const field::Field* field, double &maximum);
void atlas__NodesFunctionSpace__max_float(const NodeColumns* This, const field::Field* field, float &maximum);
void atlas__NodesFunctionSpace__max_int(const NodeColumns* This, const field::Field* field, int &maximum);
void atlas__NodesFunctionSpace__max_long(const NodeColumns* This, const field::Field* field, long &maximum);

void atlas__NodesFunctionSpace__min_arr_double(const NodeColumns* This, const field::Field* field, double* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_float(const NodeColumns* This, const field::Field* field, float* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_int(const NodeColumns* This, const field::Field* field, int* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_long(const NodeColumns* This, const field::Field* field, long* &minimum, int &size);

void atlas__NodesFunctionSpace__max_arr_double(const NodeColumns* This, const field::Field* field, double* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_float(const NodeColumns* This, const field::Field* field, float* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_int(const NodeColumns* This, const field::Field* field, int* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_long(const NodeColumns* This, const field::Field* field, long* &maximum, int &size);

void atlas__NodesFunctionSpace__minloc_double(const NodeColumns* This, const field::Field* field, double &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_float(const NodeColumns* This, const field::Field* field, float &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_int(const NodeColumns* This, const field::Field* field, int &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_long(const NodeColumns* This, const field::Field* field, long &minimum, long &glb_idx);

void atlas__NodesFunctionSpace__maxloc_double(const NodeColumns* This, const field::Field* field, double &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_float(const NodeColumns* This, const field::Field* field, float &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_int(const NodeColumns* This, const field::Field* field, int &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_long(const NodeColumns* This, const field::Field* field, long &maximum, long &glb_idx);

void atlas__NodesFunctionSpace__minloc_arr_double(const NodeColumns* This, const field::Field* field, double* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_float(const NodeColumns* This, const field::Field* field, float* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_int(const NodeColumns* This, const field::Field* field, int* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_long(const NodeColumns* This, const field::Field* field, long* &minimum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodeColumns* This, const field::Field* field, double* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_float(const NodeColumns* This, const field::Field* field, float* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_int(const NodeColumns* This, const field::Field* field, int* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_long(const NodeColumns* This, const field::Field* field, long* &maximum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__mean_double(const NodeColumns* This, const field::Field* field, double &mean, int &N);
void atlas__NodesFunctionSpace__mean_float(const NodeColumns* This, const field::Field* field, float &mean, int &N);
void atlas__NodesFunctionSpace__mean_int(const NodeColumns* This, const field::Field* field, int &mean, int &N);
void atlas__NodesFunctionSpace__mean_long(const NodeColumns* This, const field::Field* field, long &mean, int &N);
void atlas__NodesFunctionSpace__mean_arr_double(const NodeColumns* This, const field::Field* field, double* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_float(const NodeColumns* This, const field::Field* field, float* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_int(const NodeColumns* This, const field::Field* field, int* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_long(const NodeColumns* This, const field::Field* field, long* &mean, int &size, int &N);

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodeColumns* This, const field::Field* field, double &mean, double &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodeColumns* This, const field::Field* field, float &mean, float &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodeColumns* This, const field::Field* field, int &mean, int &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodeColumns* This, const field::Field* field, long &mean, long &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodeColumns* This, const field::Field* field, double* &mean, double* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodeColumns* This, const field::Field* field, float* &mean, float* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodeColumns* This, const field::Field* field, int* &mean, int* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodeColumns* This, const field::Field* field, long* &mean, long* &stddev, int &size, int &N);


void atlas__NodesFunctionSpace__minloclev_double(const NodeColumns* This, const field::Field* field, double &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_float(const NodeColumns* This, const field::Field* field, float &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_int(const NodeColumns* This, const field::Field* field, int &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_long(const NodeColumns* This, const field::Field* field, long &minimum, long &glb_idx, int &level);

void atlas__NodesFunctionSpace__maxloclev_double(const NodeColumns* This, const field::Field* field, double &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_float(const NodeColumns* This, const field::Field* field, float &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_int(const NodeColumns* This, const field::Field* field, int &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_long(const NodeColumns* This, const field::Field* field, long &maximum, long &glb_idx, int &level);

void atlas__NodesFunctionSpace__minloclev_arr_double(const NodeColumns* This, const field::Field* field, double* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_float(const NodeColumns* This, const field::Field* field, float* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_int(const NodeColumns* This, const field::Field* field, int* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_long(const NodeColumns* This, const field::Field* field, long* &minimum, long* &glb_idx, int* &level, int &size);

void atlas__NodesFunctionSpace__maxloclev_arr_double(const NodeColumns* This, const field::Field* field, double* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_float(const NodeColumns* This, const field::Field* field, float* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_int(const NodeColumns* This, const field::Field* field, int* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_long(const NodeColumns* This, const field::Field* field, long* &maximum, long* &glb_idx, int* &level, int &size);

void atlas__NodesFunctionSpace__sum_per_level(const NodeColumns* This, const field::Field* field, field::Field* sum, int &N);
void atlas__NodesFunctionSpace__oisum_per_level(const NodeColumns* This, const field::Field* field, field::Field* sum, int &N);
void atlas__NodesFunctionSpace__min_per_level(const NodeColumns* This, const field::Field* field, field::Field* min);
void atlas__NodesFunctionSpace__max_per_level(const NodeColumns* This, const field::Field* field, field::Field* max);
void atlas__NodesFunctionSpace__minloc_per_level(const NodeColumns* This, const field::Field* field, field::Field* min, field::Field* glb_idx);
void atlas__NodesFunctionSpace__maxloc_per_level(const NodeColumns* This, const field::Field* field, field::Field* max, field::Field* glb_idx);
void atlas__NodesFunctionSpace__mean_per_level(const NodeColumns* This, const field::Field* field, field::Field* mean, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_per_level(const NodeColumns* This, const field::Field* field, field::Field* mean, field::Field* stddev, int &N);



}

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodeColumnsFunctionSpaceInterface_h
