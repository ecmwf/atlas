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

#include "atlas/internals/atlas_defines.h"
#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "atlas/array/gridtools/GridToolsTraits.h"
#endif

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

    StorageView(storage_view_t storage_view, size_t size, bool contiguous = true) :
        gt_storage_view_(storage_view),
        size_(size),
        contiguous_(contiguous)
    {}

    DATA_TYPE* data() { return gt_storage_view_.data(); }

    size_t size() { return size_; }

    void assign(const DATA_TYPE& value) {
       ASSERT( contiguous() );
       DATA_TYPE* raw_data = data();
       for( size_t j=0; j<size_; ++j ) {
         raw_data[j] = value;
       }
    }

    bool contiguous() const
    {
      return contiguous_;
    }

private:
    storage_view_t gt_storage_view_;
    size_t size_;
    bool contiguous_;
};

//------------------------------------------------------------------------------------------------------
#else
//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class StorageView
{
  typedef void* storage_view_t;
public:
    StorageView(storage_view_t storage_view, size_t size, bool contiguous = true) :
        native_storage_view_(storage_view),
        size_(size),
        contiguous_(contiguous)
    {}

    DATA_TYPE* data() { return (DATA_TYPE*) native_storage_view_; }

    size_t size() { return size_; }

    void assign(const DATA_TYPE& value) {
       ASSERT( contiguous() );
       DATA_TYPE* raw_data = data();
       for( size_t j=0; j<size_; ++j ) {
         raw_data[j] = value;
       }
    }

    bool contiguous() const
    {
      return contiguous_;
    }

private:
    storage_view_t native_storage_view_;
    size_t size_;
    bool contiguous_;
};

//------------------------------------------------------------------------------------------------------

#endif

} // namespace array
} // namespace atlas
