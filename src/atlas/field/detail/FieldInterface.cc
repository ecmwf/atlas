/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "atlas/library/config.h"
#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/field/detail/FieldInterface.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_FUNCTIONSPACE
#include "atlas/functionspace/FunctionSpace.h"
#endif
namespace atlas {
namespace field {

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

namespace {
template <typename Value>
void atlas__Field__data_specf(FieldImpl* This, Value*& data, int& rank, int*& shapef, int*& stridesf) {
    ATLAS_ASSERT(This != nullptr, "Cannot access data of uninitialised atlas_Field");
    if (This->datatype() != array::make_datatype<Value>()) {
        throw_Exception("Datatype mismatch for accessing field data");
    }
    This->array().accMap();
    data     = This->data<Value>();
    shapef   = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank     = This->shapef().size();
}

template <typename Value>
FieldImpl* atlas__Field__wrap_specf(const char* name, Value data[], int rank, int shapef[], int stridesf[]) {
    array::ArrayShape shape;
    shape.resize(rank);
    array::ArrayStrides strides;
    strides.resize(rank);
    idx_t jf = rank - 1;
    for (int j = 0; j < rank; ++j) {
        shape[j]   = shapef[jf];
        strides[j] = stridesf[jf];
        --jf;
    }
    FieldImpl* field;
    {
        Field wrapped(std::string(name), data, array::ArraySpec(shape, strides));
        field = wrapped.get();
        field->attach();
    }
    field->detach();
    ATLAS_ASSERT(field);
    return field;
}


}  // namespace

extern "C" {

FieldImpl* atlas__Field__wrap_int_specf(const char* name, int data[], int rank, int shapef[], int stridesf[]) {
    return atlas__Field__wrap_specf(name, data, rank, shapef, stridesf);
}

FieldImpl* atlas__Field__wrap_long_specf(const char* name, long data[], int rank, int shapef[], int stridesf[]) {
    return atlas__Field__wrap_specf(name, data, rank, shapef, stridesf);
}

FieldImpl* atlas__Field__wrap_float_specf(const char* name, float data[], int rank, int shapef[], int stridesf[]) {
    return atlas__Field__wrap_specf(name, data, rank, shapef, stridesf);
}

FieldImpl* atlas__Field__wrap_double_specf(const char* name, double data[], int rank, int shapef[], int stridesf[]) {
    return atlas__Field__wrap_specf(name, data, rank, shapef, stridesf);
}

FieldImpl* atlas__Field__create(eckit::Parametrisation* params) {
    ATLAS_ASSERT(params != nullptr);
    FieldImpl* field;
    {
        Field f(*params);
        field = f.get();
        field->attach();
    }
    field->detach();
    ATLAS_ASSERT(field);
    return field;
}

void atlas__Field__delete(FieldImpl* This) {
    delete This;
}

const char* atlas__Field__name(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access name of uninitialised atlas_Field");
    return This->name().c_str();
}

void atlas__Field__datatype(FieldImpl* This, char*& datatype, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access datatype of uninitialised atlas_Field");
    std::string s = This->datatype().str();
    size          = static_cast<int>(s.size());
    datatype      = new char[size + 1];
    std::strncpy(datatype, s.c_str(), size + 1);
    allocated = true;
}

int atlas__Field__size(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access size of uninitialised atlas_Field");
    return This->size();
}

int atlas__Field__rank(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access rank of uninitialised atlas_Field");
    return This->rank();
}

int atlas__Field__kind(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access data kind of uninitialised atlas_Field");
    return This->datatype().kind();
}

double atlas__Field__bytes(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access bytes occupied by uninitialised atlas_Field");
    return This->bytes();
}

int atlas__Field__levels(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access levels of uninitialised atlas_Field");
    return This->levels();
}

util::Metadata* atlas__Field__metadata(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access metadata of uninitialised atlas_Field");
    return &This->metadata();
}

int atlas__Field__has_functionspace(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
#if ATLAS_HAVE_FUNCTIONSPACE
    return (This->functionspace() != 0);
#else
    return 0;
#endif
}

const functionspace::FunctionSpaceImpl* atlas__Field__functionspace(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
#if ATLAS_HAVE_FUNCTIONSPACE
    return This->functionspace().get();
#else
    throw_Exception("Atlas is not compiled with FunctionSpace support", Here());
#endif
}

void atlas__Field__shapef(FieldImpl* This, int*& shape, int& rank) {
    ATLAS_ASSERT(This != nullptr, "Cannot access bytes occupied by uninitialised atlas_Field");
    shape = const_cast<int*>(&This->shapef().front());
    rank  = This->shapef().size();
}

void atlas__Field__stridesf(FieldImpl* This, int*& shape, int& rank) {
    ATLAS_ASSERT(This != nullptr, "Cannot access bytes occupied by uninitialised atlas_Field");
    shape = const_cast<int*>(&This->stridesf().front());
    rank  = This->stridesf().size();
}

void atlas__Field__data_int_specf(FieldImpl* This, int*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_specf(This, data, rank, shapef, stridesf);
}

void atlas__Field__data_long_specf(FieldImpl* This, long*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_specf(This, data, rank, shapef, stridesf);
}

void atlas__Field__data_float_specf(FieldImpl* This, float*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_specf(This, data, rank, shapef, stridesf);
}

void atlas__Field__data_double_specf(FieldImpl* This, double*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_specf(This, data, rank, shapef, stridesf);
}

int atlas__Field__host_needs_update(const FieldImpl* This) {
    return This->hostNeedsUpdate();
}

int atlas__Field__device_needs_update(const FieldImpl* This) {
    return This->deviceNeedsUpdate();
}

int atlas__Field__device_allocated(const FieldImpl* This) {
    return This->deviceAllocated();
}

void atlas__Field__rename(FieldImpl* This, const char* name) {
    ATLAS_ASSERT(This, "Cannot rename uninitialised atlas_Field");
    This->rename(std::string(name));
}

void atlas__Field__set_levels(FieldImpl* This, int levels) {
    ATLAS_ASSERT(This != nullptr, "Cannot set levels of uninitialised atlas_Field");
    This->set_levels(levels);
}

#
void atlas__Field__set_functionspace(FieldImpl* This, const functionspace::FunctionSpaceImpl* functionspace) {
    ATLAS_ASSERT(This != nullptr, "Cannot set functionspace in uninitialised atlas_Field");
    ATLAS_ASSERT(functionspace != nullptr, "Cannot set uninitialised atlas_FunctionSpace in atlas_Field");
#if ATLAS_HAVE_FUNCTIONSPACE
    This->set_functionspace(functionspace);
#else
    throw_Exception("Atlas is not compiled with FunctionSpace support", Here());
#endif
}

void atlas__Field__set_host_needs_update(const FieldImpl* This, int value) {
    ATLAS_ASSERT(This != nullptr, "Cannot set value for uninitialised atlas_Field");
    This->setHostNeedsUpdate(value);
}

void atlas__Field__set_device_needs_update(const FieldImpl* This, int value) {
    ATLAS_ASSERT(This != nullptr, "Cannot set value for uninitialised atlas_Field");
    This->setDeviceNeedsUpdate(value);
}

void atlas__Field__update_device(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->updateDevice();
}

void atlas__Field__update_host(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->updateHost();
}

void atlas__Field__sync_host_device(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->syncHostDevice();
}

void atlas__Field__allocate_device(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->allocateDevice();
}

void atlas__Field__deallocate_device(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->deallocateDevice();
}

void atlas__Field__set_dirty(FieldImpl* This, int value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    This->set_dirty(value);
}

int atlas__Field__dirty(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    return This->dirty();
}

int atlas__Field__contiguous(FieldImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Field");
    return This->array().contiguous();
}

void atlas__Field__halo_exchange(FieldImpl* This, int on_device) {
    ATLAS_ASSERT(This != nullptr, "Cannot halo-exchange uninitialised atlas_Field");
    This->haloExchange(on_device);
}

void atlas__Field__adjoint_halo_exchange(FieldImpl* This, int on_device) {
    ATLAS_ASSERT(This != nullptr, "Cannot adjoint-halo-exchange uninitialised atlas_Field");
    This->adjointHaloExchange(on_device);
}
}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
