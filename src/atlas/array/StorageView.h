/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#pragma once

#include "atlas/array/GridToolsTraits.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class StorageView
{
public:

  using storage_view_t = data_view_tt<DATA_TYPE, 1>;

public:

    StorageView(storage_view_t storage_view) : gt_storage_view_(storage_view) {}

    DATA_TYPE* data() { return gt_storage_view_.data(); }

private:
    storage_view_t gt_storage_view_;
};

//------------------------------------------------------------------------------------------------------
#else
//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class StorageView
{
  typedef void* storage_view_t;
public:
    StorageView(storage_view_t storage_view) : native_storage_view_(storage_view) {}

    DATA_TYPE* data() { return (DATA_TYPE*) native_storage_view_; }

private:
    storage_view_t native_storage_view_;
};

//------------------------------------------------------------------------------------------------------

#endif

} // namespace array
} // namespace atlas
