/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_NodesFunctionSpaceInterface_h
#define atlas_functionspace_NodesFunctionSpaceInterface_h


#include "atlas/functionspace/NodesFunctionSpace.h"

namespace atlas {
namespace functionspace {

#define Char char
extern "C" {
NodesFunctionSpace* atlas__NodesFunctionSpace__new (const char* name, Mesh* mesh, int halo);
void atlas__NodesFunctionSpace__delete (NodesFunctionSpace* This);
int atlas__NodesFunctionSpace__nb_nodes(const NodesFunctionSpace* This);
Mesh* atlas__NodesFunctionSpace__mesh(NodesFunctionSpace* This);
Nodes* atlas__NodesFunctionSpace__nodes(NodesFunctionSpace* This);
Field* atlas__NodesFunctionSpace__create_field (const NodesFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesFunctionSpace__create_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__NodesFunctionSpace__create_field_lev (const NodesFunctionSpace* This, const char* name, int levels, int kind);
Field* atlas__NodesFunctionSpace__create_field_lev_vars (const NodesFunctionSpace* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind);


Field* atlas__NodesFunctionSpace__create_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template);
Field* atlas__NodesFunctionSpace__create_global_field (const NodesFunctionSpace* This, const char* name, int kind);
Field* atlas__NodesFunctionSpace__create_global_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__NodesFunctionSpace__create_global_field_lev (const NodesFunctionSpace* This, const char* name, int levels, int kind);
Field* atlas__NodesFunctionSpace__create_global_field_lev_vars (const NodesFunctionSpace* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__NodesFunctionSpace__create_global_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template);

void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodesFunctionSpace* This, FieldSet* fieldset);
void atlas__NodesFunctionSpace__halo_exchange_field(const NodesFunctionSpace* This, Field* field);

void atlas__NodesFunctionSpace__gather_fieldset(const NodesFunctionSpace* This, const FieldSet* local, FieldSet* global);
void atlas__NodesFunctionSpace__gather_field(const NodesFunctionSpace* This, const Field* local, Field* global);

void atlas__NodesFunctionSpace__scatter_fieldset(const NodesFunctionSpace* This, const FieldSet* global, FieldSet* local);
void atlas__NodesFunctionSpace__scatter_field(const NodesFunctionSpace* This, const Field* global, Field* local);
void atlas__NodesFunctionSpace__checksum_fieldset(const NodesFunctionSpace* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
void atlas__NodesFunctionSpace__checksum_field(const NodesFunctionSpace* This, const Field* field, Char* &checksum, int &size, int &allocated);

void atlas__NodesFunctionSpace__sum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__sum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__sum_int(const NodesFunctionSpace* This, const Field* field, int &sum, int &N);
void atlas__NodesFunctionSpace__sum_long(const NodesFunctionSpace* This, const Field* field, long &sum, int &N);
void atlas__NodesFunctionSpace__sum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_int(const NodesFunctionSpace* This, const Field* field, int* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__sum_arr_long(const NodesFunctionSpace* This, const Field* field, long* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__oisum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N);
void atlas__NodesFunctionSpace__oisum_float(const NodesFunctionSpace* This, const Field* field, float &sum, int &N);
void atlas__NodesFunctionSpace__oisum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N);
void atlas__NodesFunctionSpace__oisum_arr_float(const NodesFunctionSpace* This, const Field* field, float* &sum, int &size, int &N);

void atlas__NodesFunctionSpace__min_double(const NodesFunctionSpace* This, const Field* field, double &minimum);
void atlas__NodesFunctionSpace__min_float(const NodesFunctionSpace* This, const Field* field, float &minimum);
void atlas__NodesFunctionSpace__min_int(const NodesFunctionSpace* This, const Field* field, int &minimum);
void atlas__NodesFunctionSpace__min_long(const NodesFunctionSpace* This, const Field* field, long &minimum);

void atlas__NodesFunctionSpace__max_double(const NodesFunctionSpace* This, const Field* field, double &maximum);
void atlas__NodesFunctionSpace__max_float(const NodesFunctionSpace* This, const Field* field, float &maximum);
void atlas__NodesFunctionSpace__max_int(const NodesFunctionSpace* This, const Field* field, int &maximum);
void atlas__NodesFunctionSpace__max_long(const NodesFunctionSpace* This, const Field* field, long &maximum);

void atlas__NodesFunctionSpace__min_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, int &size);
void atlas__NodesFunctionSpace__min_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, int &size);

void atlas__NodesFunctionSpace__max_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, int &size);
void atlas__NodesFunctionSpace__max_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, int &size);

void atlas__NodesFunctionSpace__minloc_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_float(const NodesFunctionSpace* This, const Field* field, float &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_int(const NodesFunctionSpace* This, const Field* field, int &minimum, long &glb_idx);
void atlas__NodesFunctionSpace__minloc_long(const NodesFunctionSpace* This, const Field* field, long &minimum, long &glb_idx);

void atlas__NodesFunctionSpace__maxloc_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_float(const NodesFunctionSpace* This, const Field* field, float &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_int(const NodesFunctionSpace* This, const Field* field, int &maximum, long &glb_idx);
void atlas__NodesFunctionSpace__maxloc_long(const NodesFunctionSpace* This, const Field* field, long &maximum, long &glb_idx);

void atlas__NodesFunctionSpace__minloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__minloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, long* &glb_idx, int &size);
void atlas__NodesFunctionSpace__maxloc_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, long* &glb_idx, int &size);

void atlas__NodesFunctionSpace__mean_double(const NodesFunctionSpace* This, const Field* field, double &mean, int &N);
void atlas__NodesFunctionSpace__mean_float(const NodesFunctionSpace* This, const Field* field, float &mean, int &N);
void atlas__NodesFunctionSpace__mean_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &N);
void atlas__NodesFunctionSpace__mean_long(const NodesFunctionSpace* This, const Field* field, long &mean, int &N);
void atlas__NodesFunctionSpace__mean_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int &size, int &N);
void atlas__NodesFunctionSpace__mean_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, int &size, int &N);

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodesFunctionSpace* This, const Field* field, double &mean, double &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_float(const NodesFunctionSpace* This, const Field* field, float &mean, float &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_int(const NodesFunctionSpace* This, const Field* field, int &mean, int &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_long(const NodesFunctionSpace* This, const Field* field, long &mean, long &stddev, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, double* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &mean, float* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &mean, int* &stddev, int &size, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &mean, long* &stddev, int &size, int &N);


void atlas__NodesFunctionSpace__minloclev_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_float(const NodesFunctionSpace* This, const Field* field, float &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_int(const NodesFunctionSpace* This, const Field* field, int &minimum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__minloclev_long(const NodesFunctionSpace* This, const Field* field, long &minimum, long &glb_idx, int &level);

void atlas__NodesFunctionSpace__maxloclev_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_float(const NodesFunctionSpace* This, const Field* field, float &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_int(const NodesFunctionSpace* This, const Field* field, int &maximum, long &glb_idx, int &level);
void atlas__NodesFunctionSpace__maxloclev_long(const NodesFunctionSpace* This, const Field* field, long &maximum, long &glb_idx, int &level);

void atlas__NodesFunctionSpace__minloclev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &minimum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__minloclev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &minimum, long* &glb_idx, int* &level, int &size);

void atlas__NodesFunctionSpace__maxloclev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_float(const NodesFunctionSpace* This, const Field* field, float* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_int(const NodesFunctionSpace* This, const Field* field, int* &maximum, long* &glb_idx, int* &level, int &size);
void atlas__NodesFunctionSpace__maxloclev_arr_long(const NodesFunctionSpace* This, const Field* field, long* &maximum, long* &glb_idx, int* &level, int &size);

void atlas__NodesFunctionSpace__sum_per_level(const NodesFunctionSpace* This, const Field* field, Field* sum, int &N);
void atlas__NodesFunctionSpace__oisum_per_level(const NodesFunctionSpace* This, const Field* field, Field* sum, int &N);
void atlas__NodesFunctionSpace__min_per_level(const NodesFunctionSpace* This, const Field* field, Field* min);
void atlas__NodesFunctionSpace__max_per_level(const NodesFunctionSpace* This, const Field* field, Field* max);
void atlas__NodesFunctionSpace__minloc_per_level(const NodesFunctionSpace* This, const Field* field, Field* min, Field* glb_idx);
void atlas__NodesFunctionSpace__maxloc_per_level(const NodesFunctionSpace* This, const Field* field, Field* max, Field* glb_idx);
void atlas__NodesFunctionSpace__mean_per_level(const NodesFunctionSpace* This, const Field* field, Field* mean, int &N);
void atlas__NodesFunctionSpace__mean_and_stddev_per_level(const NodesFunctionSpace* This, const Field* field, Field* mean, Field* stddev, int &N);



}
#undef Char

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpaceInterface_h
