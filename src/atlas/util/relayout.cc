/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/relayout.h"

/**
 * @file relayout.cc
 * @brief Type dispatch, FieldSet dispatch, and C bindings for Atlas relayout utilities.
 *
 * The public C++ entry points validate datatype compatibility, select the supported value
 * type, create host or device views, and dispatch to the host/device implementation files.
 */

#include <type_traits>

#include "pluto/pluto.h"

#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Log.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/mdspan.h"

namespace atlas {

#if !ATLAS_HAVE_GPU
template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    if (pluto::devices() == 0) {
        return host_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
    ATLAS_NOTIMPLEMENTED;
}

template <class Blocked, class Nonblocked>
void device_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    if (pluto::devices() == 0) {
        return host_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
    ATLAS_NOTIMPLEMENTED;
}

template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    if (pluto::devices() == 0) {
        return host_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
    ATLAS_NOTIMPLEMENTED;
}
#endif

template <class Nonblocked, class Blocked>
void copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked, bool on_device) {
    ATLAS_TRACE("copy_nonblocked_to_blocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device && pluto::devices()) {
        device_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
    else {
        host_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
}

template <class Blocked, class Nonblocked>
void copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked, bool on_device) {
    ATLAS_TRACE("copy_blocked_to_nonblocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device && pluto::devices()) {
        device_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
    else {
        host_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
}

template <class BlockedIn, class BlockedOut>
void copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out, bool on_device) {
    ATLAS_TRACE("copy_blocked_to_blocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device && pluto::devices()) {
        device_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
    else {
        host_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
}

template <class ValueType>
void copy_blocked_to_nonblocked_T(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    if (blocked.rank() != nonblocked.rank()+1) {
        ATLAS_THROW_EXCEPTION("copy_blocked_to_nonblocked_T: blocked rank must be one more than nonblocked rank" << " but got blocked rank " << blocked.rank() << " and nonblocked rank " << nonblocked.rank());
    }
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if (blocked.rank()==4) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 4>(blocked)    : array::make_host_view<ValueType, 4>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 3>(nonblocked) : array::make_host_view<ValueType, 3>(nonblocked);
        copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
    }
    else if (blocked.rank()==3) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 3>(blocked)    : array::make_host_view<ValueType, 3>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 2>(nonblocked) : array::make_host_view<ValueType, 2>(nonblocked);
        copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
    }
    else if (blocked.rank()==2) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 2>(blocked)    : array::make_host_view<ValueType, 2>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 1>(nonblocked) : array::make_host_view<ValueType, 1>(nonblocked);
        copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
    }
    else {
        ATLAS_THROW_EXCEPTION("copy_blocked_to_nonblocked not implemented for blocked.rank " + std::to_string(blocked.rank()));
    }
    nonblocked.setHostNeedsUpdate(on_device);
    nonblocked.setDeviceNeedsUpdate(!on_device);
}

template <class ValueType>
void copy_nonblocked_to_blocked_T(const array::Array& nonblocked, array::Array& blocked, bool on_device) {
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if (blocked.rank()==4) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,4>(blocked)    : array::make_host_view<ValueType, 4>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,3>(nonblocked) : array::make_host_view<ValueType, 3>(nonblocked);
        copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
    }
    else if (blocked.rank()==3) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,3>(blocked)    : array::make_host_view<ValueType, 3>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,2>(nonblocked) : array::make_host_view<ValueType, 2>(nonblocked);
        copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
    }
    else if (blocked.rank()==2) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,2>(blocked)    : array::make_host_view<ValueType, 2>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,1>(nonblocked) : array::make_host_view<ValueType, 1>(nonblocked);
        copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
    }
    else {
        ATLAS_THROW_EXCEPTION("copy_nonblocked_to_blocked not implemented for blocked.rank " + std::to_string(blocked.rank()));
    }
    blocked.setHostNeedsUpdate(on_device);
    blocked.setDeviceNeedsUpdate(!on_device);
}

template <class ValueType>
void copy_blocked_to_blocked_T(const array::Array& blocked_in, array::Array& blocked_out, bool on_device) {
    ATLAS_ASSERT(blocked_in.rank() == blocked_out.rank());
    if (blocked_in.rank()==4) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 4>(blocked_in)  : array::make_host_view<ValueType, 4>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 4>(blocked_out) : array::make_host_view<ValueType, 4>(blocked_out);
        copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
    }
    else if (blocked_in.rank()==3) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 3>(blocked_in)  : array::make_host_view<ValueType, 3>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 3>(blocked_out) : array::make_host_view<ValueType, 3>(blocked_out);
        copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
    }
    else if (blocked_in.rank()==2) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 2>(blocked_in)  : array::make_host_view<ValueType, 2>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 2>(blocked_out) : array::make_host_view<ValueType, 2>(blocked_out);
        copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
    }
    else {
        ATLAS_THROW_EXCEPTION("copy_blocked_to_blocked not implemented for blocked.rank " + std::to_string(blocked_in.rank()));
    }
    blocked_out.setHostNeedsUpdate(on_device);
    blocked_out.setDeviceNeedsUpdate(!on_device);
}

void copy_nonblocked_to_blocked(const array::Array& nonblocked, array::Array& blocked, bool on_device) {
    ATLAS_ASSERT(blocked.datatype() == nonblocked.datatype());
    switch (blocked.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_nonblocked_to_blocked_T<int>(nonblocked, blocked, on_device);
        case array::DataType::kind<long>()   : return copy_nonblocked_to_blocked_T<long>(nonblocked, blocked, on_device);
        case array::DataType::kind<float>()  : return copy_nonblocked_to_blocked_T<float>(nonblocked, blocked, on_device);
        case array::DataType::kind<double>() : return copy_nonblocked_to_blocked_T<double>(nonblocked, blocked, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_nonblocked(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    ATLAS_ASSERT(blocked.datatype() == nonblocked.datatype());
    switch (blocked.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_blocked_to_nonblocked_T<int>(blocked, nonblocked, on_device);
        case array::DataType::kind<long>()   : return copy_blocked_to_nonblocked_T<long>(blocked, nonblocked, on_device);
        case array::DataType::kind<float>()  : return copy_blocked_to_nonblocked_T<float>(blocked, nonblocked, on_device);
        case array::DataType::kind<double>() : return copy_blocked_to_nonblocked_T<double>(blocked, nonblocked, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_blocked(const array::Array& blocked_in, array::Array& blocked_out, bool on_device) {
    ATLAS_ASSERT(blocked_in.datatype() == blocked_out.datatype());
    switch (blocked_in.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_blocked_to_blocked_T<int>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<long>()   : return copy_blocked_to_blocked_T<long>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<float>()  : return copy_blocked_to_blocked_T<float>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<double>() : return copy_blocked_to_blocked_T<double>(blocked_in, blocked_out, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_blocked(const FieldSet& blocked_fields_in, FieldSet& blocked_fields_out, bool on_device) {
    ATLAS_ASSERT(blocked_fields_in.size() == blocked_fields_out.size());
    for( int i=0; i<blocked_fields_in.size(); ++i) {
        copy_blocked_to_blocked(blocked_fields_in[i], blocked_fields_out[i], on_device);
    }
}

void copy_blocked_to_nonblocked(const FieldSet& blocked_fields_in, FieldSet& nonblocked_fields_out, bool on_device) {
    ATLAS_ASSERT(blocked_fields_in.size() == nonblocked_fields_out.size());
    for( int i=0; i<blocked_fields_in.size(); ++i) {
        copy_blocked_to_nonblocked(blocked_fields_in[i], nonblocked_fields_out[i], on_device);
    }
}

void copy_nonblocked_to_blocked(const FieldSet& nonblocked_fields_in, FieldSet& blocked_fields_out, bool on_device) {
    ATLAS_ASSERT(nonblocked_fields_in.size() == blocked_fields_out.size());
    for( int i=0; i<blocked_fields_out.size(); ++i) {
        copy_nonblocked_to_blocked(nonblocked_fields_in[i], blocked_fields_out[i], on_device);
    }
}


extern "C" {
/**
 * @brief C binding for copying one blocked field to another blocked field.
 *
 * @param source Non-null pointer to a source `FieldImpl` with blocked rank 2, 3, or 4.
 * @param target Non-null pointer to a target `FieldImpl` with the same datatype and rank as
 *        `source`.  The target may use a different `nproma`.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_blocked_to_blocked_field(const field::FieldImpl* source, field::FieldImpl* target, int on_device) {
    const Field source_field(source);
    Field target_field(target);
    copy_blocked_to_blocked(static_cast<const array::Array&>(source_field),
                            static_cast<array::Array&>(target_field), on_device);
}

/**
 * @brief C binding for copying a blocked field set to a blocked field set.
 *
 * @param source Non-null pointer to a source field set.  Each field must be blocked.
 * @param target Non-null pointer to a target field set with the same number of fields as
 *        `source`; corresponding fields must satisfy the blocked-to-blocked requirements.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_blocked_to_blocked_fieldset(const field::FieldSetImpl* source, field::FieldSetImpl* target, int on_device) {
    const FieldSet source_fields(source);
    FieldSet target_fields(target);
    copy_blocked_to_blocked(source_fields, target_fields, on_device);
}

/**
 * @brief C binding for copying one blocked field to a nonblocked field.
 *
 * @param source Non-null pointer to a blocked source field with rank 2, 3, or 4.
 * @param target Non-null pointer to a nonblocked target field with matching datatype and rank
 *        one less than `source`.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_blocked_to_nonblocked_field(const field::FieldImpl* source, field::FieldImpl* target, int on_device) {
    const Field source_field(source);
    Field target_field(target);
    copy_blocked_to_nonblocked(static_cast<const array::Array&>(source_field),
                               static_cast<array::Array&>(target_field), on_device);
}

/**
 * @brief C binding for copying a blocked field set to a nonblocked field set.
 *
 * @param source Non-null pointer to a source field set.  Each field must be blocked.
 * @param target Non-null pointer to a target field set with the same number of fields as
 *        `source`; corresponding fields must satisfy the blocked-to-nonblocked requirements.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_blocked_to_nonblocked_fieldset(const field::FieldSetImpl* source, field::FieldSetImpl* target, int on_device) {
    const FieldSet source_fields(source);
    FieldSet target_fields(target);
    copy_blocked_to_nonblocked(source_fields, target_fields, on_device);
}

/**
 * @brief C binding for copying one nonblocked field to a blocked field.
 *
 * @param source Non-null pointer to a nonblocked source field with rank 1, 2, or 3.
 * @param target Non-null pointer to a blocked target field with matching datatype and rank one
 *        greater than `source`.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_nonblocked_to_blocked_field(const field::FieldImpl* source, field::FieldImpl* target, int on_device) {
    const Field source_field(source);
    Field target_field(target);
    copy_nonblocked_to_blocked(static_cast<const array::Array&>(source_field),
                               static_cast<array::Array&>(target_field), on_device);
}

/**
 * @brief C binding for copying a nonblocked field set to a blocked field set.
 *
 * @param source Non-null pointer to a source field set.  Each field must be nonblocked.
 * @param target Non-null pointer to a target field set with the same number of fields as
 *        `source`; corresponding fields must satisfy the nonblocked-to-blocked requirements.
 * @param on_device Non-zero selects device execution; zero selects host execution.
 */
void atlas__copy_nonblocked_to_blocked_fieldset(const field::FieldSetImpl* source, field::FieldSetImpl* target, int on_device) {
    const FieldSet source_fields(source);
    FieldSet target_fields(target);
    copy_nonblocked_to_blocked(source_fields, target_fields, on_device);
}
}


}  // namespace atlas
