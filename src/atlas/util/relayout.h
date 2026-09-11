/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

/**
 * @file relayout.h
 * @brief Utilities for copying fields between blocked and nonblocked Atlas array layouts.
 *
 * The blocked layout stores horizontal points as blocks with the block index first and the
 * `nproma` lane index last.  Supported blocked ranks are 2, 3, and 4:
 * - rank 2: `[nblk, nproma]`
 * - rank 3: `[nblk, nlev, nproma]` or `[nblk, nvar, nproma]`
 * - rank 4: `[nblk, nvar, nlev, nproma]`
 *
 * The corresponding nonblocked layouts have one fewer dimension and store the horizontal
 * point index first:
 * - rank 1: `[npoint]`
 * - rank 2: `[npoint, nlev]` or `[npoint, nvar]`
 * - rank 3: `[npoint, nlev, nvar]`
 */

// Forward declarations to avoid including headers in this file

namespace atlas::array {
    class Array;
}
namespace atlas {
    class FieldSet;
}

namespace atlas {

/**
 * @brief Copy a nonblocked array into a blocked array.
 *
 * @param nonblocked Source array in nonblocked layout.  Its rank must be one less than
 *        `blocked.rank()`.  For rank-3 nonblocked data the dimensions must be
 *        `[npoint, nlev, nvar]`.
 * @param blocked Target array in blocked layout.  Its datatype must match `nonblocked`.
 *        Supported ranks are 2, 3, and 4.  For rank-4 blocked data the dimensions must be
 *        `[nblk, nvar, nlev, nproma]`.
 * @param on_device If true, both arrays must have valid device views and the copy is
 *        dispatched to the device implementation.  If false, host views are used.
 *
 * @pre `blocked.datatype() == nonblocked.datatype()`.
 * @pre `blocked.rank() == nonblocked.rank() + 1`.
 * @pre Supported datatypes are `int`, `long`, `float`, and `double`.
 */
void copy_nonblocked_to_blocked(const array::Array& nonblocked, array::Array& blocked, bool on_device);

/**
 * @brief Copy a blocked array into a nonblocked array.
 *
 * @param blocked Source array in blocked layout.  Supported ranks are 2, 3, and 4.  For
 *        rank-4 blocked data the dimensions must be `[nblk, nvar, nlev, nproma]`.
 * @param nonblocked Target array in nonblocked layout.  Its datatype must match `blocked`,
 *        and its rank must be one less than `blocked.rank()`.  For rank-3 nonblocked data
 *        the dimensions must be `[npoint, nlev, nvar]`.
 * @param on_device If true, both arrays must have valid device views and the copy is
 *        dispatched to the device implementation.  If false, host views are used.
 *
 * @pre `blocked.datatype() == nonblocked.datatype()`.
 * @pre `blocked.rank() == nonblocked.rank() + 1`.
 * @pre Supported datatypes are `int`, `long`, `float`, and `double`.
 */
void copy_blocked_to_nonblocked(const array::Array& blocked, array::Array& nonblocked, bool on_device);

/**
 * @brief Copy between two blocked arrays, allowing different `nproma` values.
 *
 * @param blocked_in Source array in blocked layout.  Supported ranks are 2, 3, and 4.
 * @param blocked_out Target array in blocked layout.  Its rank and datatype must match
 *        `blocked_in`.  The output `nproma` and number of blocks may differ from the input.
 * @param on_device If true, both arrays must have valid device views and the copy is
 *        dispatched to the device implementation.  If false, host views are used.
 *
 * @pre `blocked_in.datatype() == blocked_out.datatype()`.
 * @pre `blocked_in.rank() == blocked_out.rank()`.
 * @pre Supported datatypes are `int`, `long`, `float`, and `double`.
 */
void copy_blocked_to_blocked(const array::Array& blocked_in, array::Array& blocked_out, bool on_device);

/**
 * @brief Copy all fields in a blocked field set to a blocked field set.
 *
 * @param blocked_fields_in Source field set.  Each field must satisfy the blocked-array
 *        requirements of copy_blocked_to_blocked().
 * @param blocked_fields_out Target field set.  It must have the same number of fields as
 *        `blocked_fields_in`; each corresponding field must have matching rank and datatype.
 * @param on_device Selects host or device copy for every field.
 */
void copy_blocked_to_blocked(const FieldSet& blocked_fields_in, FieldSet& blocked_fields_out, bool on_device);

/**
 * @brief Copy all fields in a blocked field set to a nonblocked field set.
 *
 * @param blocked_fields_in Source field set.  Each field must satisfy the blocked-array
 *        requirements of copy_blocked_to_nonblocked().
 * @param nonblocked_fields_out Target field set.  It must have the same number of fields as
 *        `blocked_fields_in`; each corresponding field must have matching datatype and rank
 *        one less than the blocked source field.
 * @param on_device Selects host or device copy for every field.
 */
void copy_blocked_to_nonblocked(const FieldSet& blocked_fields_in, FieldSet& nonblocked_fields_out, bool on_device);

/**
 * @brief Copy all fields in a nonblocked field set to a blocked field set.
 *
 * @param nonblocked_fields_in Source field set.  Each field must satisfy the nonblocked-array
 *        requirements of copy_nonblocked_to_blocked().
 * @param blocked_fields_out Target field set.  It must have the same number of fields as
 *        `nonblocked_fields_in`; each corresponding field must have matching datatype and rank
 *        one greater than the nonblocked source field.
 * @param on_device Selects host or device copy for every field.
 */
void copy_nonblocked_to_blocked(const FieldSet& nonblocked_fields_in, FieldSet& blocked_fields_out, bool on_device);

// Implementation in relayout_on_host.cc
/// @brief Host implementation for copying a nonblocked view to a blocked view.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked);

/// @brief Host implementation for copying a blocked view to a nonblocked view.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked);

/// @brief Host implementation for copying between blocked views.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class BlockedIn, class BlockedOut>
void host_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out);

// Implementation in relayout_on_device.hic
/// @brief Device implementation for copying a nonblocked view to a blocked view.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked);

/// @brief Device implementation for copying a blocked view to a nonblocked view.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class Blocked, class Nonblocked>
void device_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked);

/// @brief Device implementation for copying between blocked views.
///
/// The view arguments may be Atlas views (`atlas::View`/`atlas::ArrayView`) or
/// mdspan-like views that provide `rank()`, `extent(i)`, element access, and a
/// compatible `value_type`.
template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out);

}
