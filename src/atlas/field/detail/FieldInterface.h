/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file FieldInterface.h
/// @author Willem Deconinck
/// @date Sep 2014

#pragma once

#include "atlas/field/detail/FieldImpl.h"

namespace atlas {
namespace functionspace {
class FunctionSpaceImpl;
}
}  // namespace atlas

namespace atlas {
namespace field {

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
FieldImpl* atlas__Field__wrap_int_specf(const char* name, int data[], int rank, int shapef[], int stridesf[]);
FieldImpl* atlas__Field__wrap_long_specf(const char* name, long data[], int rank, int shapef[], int stridesf[]);
FieldImpl* atlas__Field__wrap_float_specf(const char* name, float data[], int rank, int shapef[], int stridesf[]);
FieldImpl* atlas__Field__wrap_double_specf(const char* name, double data[], int rank, int shapef[], int stridesf[]);
FieldImpl* atlas__Field__create(eckit::Parametrisation* params);
void atlas__Field__delete(FieldImpl* This);
const char* atlas__Field__name(FieldImpl* This);
void atlas__Field__datatype(FieldImpl* This, char*& datatype, int& size, int& allocated);
int atlas__Field__kind(FieldImpl* This);
int atlas__Field__rank(FieldImpl* This);
int atlas__Field__size(FieldImpl* This);
int atlas__Field__levels(FieldImpl* This);
double atlas__Field__bytes(FieldImpl* This);
void atlas__Field__shapef(FieldImpl* This, int*& shape, int& rank);
void atlas__Field__stridesf(FieldImpl* This, int*& strides, int& rank);
void atlas__Field__data_int_specf(FieldImpl* This, int*& field_data, int& rank, int*& field_shapef,
                                  int*& field_stridesf);
void atlas__Field__data_long_specf(FieldImpl* This, long*& field_data, int& rank, int*& field_shapef,
                                   int*& field_stridesf);
void atlas__Field__data_float_specf(FieldImpl* This, float*& field_data, int& rank, int*& field_shapef,
                                    int*& field_stridesf);
void atlas__Field__data_double_specf(FieldImpl* This, double*& field_data, int& rank, int*& field_shapef,
                                     int*& field_stridesf);
void atlas__Field__device_data_int_specf(FieldImpl* This, int*& field_data, int& rank, int*& field_shapef,
                                  int*& field_stridesf);
void atlas__Field__device_data_long_specf(FieldImpl* This, long*& field_data, int& rank, int*& field_shapef,
                                   int*& field_stridesf);
void atlas__Field__device_data_float_specf(FieldImpl* This, float*& field_data, int& rank, int*& field_shapef,
                                    int*& field_stridesf);
void atlas__Field__device_data_double_specf(FieldImpl* This, double*& field_data, int& rank, int*& field_shapef,
                                     int*& field_stridesf);
util::Metadata* atlas__Field__metadata(FieldImpl* This);
const functionspace::FunctionSpaceImpl* atlas__Field__functionspace(FieldImpl* This);
void atlas__Field__rename(FieldImpl* This, const char* name);
void atlas__Field__set_levels(FieldImpl* This, int levels);
void atlas__Field__set_functionspace(FieldImpl* This, const functionspace::FunctionSpaceImpl* functionspace);
int atlas__Field__host_needs_update(const FieldImpl* This);
int atlas__Field__device_needs_update(const FieldImpl* This);
int atlas__Field__device_allocated(const FieldImpl* This);
void atlas__Field__set_host_needs_update(const FieldImpl* This, int value);
void atlas__Field__set_device_needs_update(const FieldImpl* This, int value);
void atlas__Field__update_device(FieldImpl* This);
void atlas__Field__update_host(FieldImpl* This);
void atlas__Field__sync_host_device(FieldImpl* This);
void atlas__Field__allocate_device(FieldImpl* This);
void atlas__Field__deallocate_device(FieldImpl* This);
void atlas__Field__set_dirty(FieldImpl* This, int value);
void atlas__Field__halo_exchange(FieldImpl* This, int on_device);
void atlas__Field__adjoint_halo_exchange(FieldImpl* This, int on_device);
int atlas__Field__dirty(FieldImpl* This);
int atlas__Field__contiguous(FieldImpl* This);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
