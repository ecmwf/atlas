/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_DataStoreWrapper_h
#define atlas_DataStoreWrapper_h

#include "atlas/internals/atlas_defines.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
    struct DataStoreInterface
    {
        virtual ~DataStoreInterface() {}
        virtual void clone_to_device() const = 0;
        virtual void clone_from_device() const = 0;
        virtual bool valid() const = 0;
        virtual void sync() const = 0;
        virtual bool is_on_host() const = 0;
        virtual bool is_on_device() const = 0;
        virtual void reactivate_device_write_views() const = 0;
        virtual void reactivate_host_write_views() const = 0;
        virtual void* void_data_store() = 0;
    };

    template<typename DataStore>
    struct DataStoreWrapper : DataStoreInterface
    {
        explicit DataStoreWrapper(DataStore const *ds) : data_store_(ds) {}

        ~DataStoreWrapper() {
            assert(data_store_);
            delete data_store_;
        }

        void clone_to_device() const {
            data_store_->clone_to_device();
        }
        void clone_from_device() const {
            data_store_->clone_from_device();
        }

        bool valid() const {
            return data_store_->valid();
        }
        void sync() const {
            data_store_->sync();
        }
        bool is_on_host() const {
            return data_store_->is_on_host();
        }
        bool is_on_device() const {
            return data_store_->is_on_device();
        }
        void reactivate_device_write_views() const {
            data_store_->reactivate_device_write_views();
        }
        void reactivate_host_write_views() const {
            data_store_->reactivate_host_write_views();
        }

        void* void_data_store() {
            return static_cast<void*>(const_cast<DataStore*>(data_store_));
        }

    private:
        DataStore const* data_store_;
    };

#endif

}
}
#endif
