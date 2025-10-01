/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <numeric>

#include "hdf5.h"

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {

class HDF5_Dataset {
public:
    // Callback function for H5Literate
    static herr_t list_datasets(hid_t loc_id, const char *name, const H5L_info_t *info, void *op_data) {
        H5O_info_t object_info;
        H5Oget_info_by_name(loc_id, name, &object_info, H5O_INFO_ALL, H5P_DEFAULT);

        if (object_info.type == H5O_TYPE_DATASET) {
            Log::error() << " - " << name << '\n';
        } else if (object_info.type == H5O_TYPE_GROUP) {
            // Recurse into subgroups
            hid_t group_id = H5Gopen(loc_id, name, H5P_DEFAULT);
            H5Literate(group_id, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, list_datasets, nullptr);
            H5Gclose(group_id);
        }

        return 0;
    }

    HDF5_Dataset(::hid_t file_id, const std::string& dataset_name) {
        name_ = dataset_name;
        dataset_id = H5Dopen(file_id, dataset_name.c_str(), H5P_DEFAULT);
        if (dataset_id < 0) {
            Log::error() << "ERROR: Failed to open dataset '"<<dataset_name<<"'"<<std::endl;
            Log::error() << " Existing datasets: \n";
            H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, list_datasets, nullptr);
        }
        else {
            dataspace_id = H5Dget_space(dataset_id);
            if (dataspace_id < 0) {
                Log::error() << "ERROR: Failed to extract dataspace from dataset '"<<dataset_name<<"'"<<std::endl;
                H5Dclose(dataset_id);
            }
        }
    }
    bool is_open() const {
        return dataset_id >= 0;
    }
    ~HDF5_Dataset() {
        if (dataspace_id >= 0) {
            H5Sclose(dataspace_id);
        }
        if (is_open()) {
            H5Dclose(dataset_id);
        }
    }

    atlas::DataType datatype() const {
        ATLAS_ASSERT(is_open());
        auto datatype_id = H5Dget_type(dataset_id);
        auto dt = to_datatype(datatype_id);
        H5Tclose(datatype_id);
        return dt;
    }

    int rank() const {
        int ndims = H5Sget_simple_extent_ndims(dataspace_id);
        return ndims;
    }
    std::vector<size_t> extents() const {
        std::vector<hsize_t> dims(rank());
        H5Sget_simple_extent_dims(dataspace_id, dims.data(), nullptr);
        return std::vector<size_t>(dims.begin(), dims.end());
    }

    size_t size() const {
        auto dims = extents();
        return std::accumulate(dims.begin(), dims.end(), 1., std::multiplies{});
    }

    template<class T, class Extents, class Layout = layout_right, class Accessor = default_accessor<T>>
    void read(mdspan<T,Extents,Layout,Accessor> data) {
        ATLAS_ASSERT( datatype() == make_datatype<T>());
        ATLAS_ASSERT( rank() == data.rank() );
        ATLAS_ASSERT( size() == data.size() );
        bool contiguous = [&]() {
            int s = data.stride(data.rank()-1);
            if (s != 1) {
                return false;
            }
            for( int i=int(data.rank())-2; i>= 0; --i) {
                s *= data.extent(i+1);
                if (data.stride(i) != s) {
                    return false;
                }
            }
            return true;
        }();
        ATLAS_ASSERT(contiguous);

        /////////// Experiment to see if there's a problem with rowmajor vs colmajor
        // if constexpr( data.rank() == 3) {
        //     std::vector<T> raw(size());
        //     auto status = H5Dread(dataset_id, to_h5_datatype(make_datatype<T>()), H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
        //     ATLAS_ASSERT(status >= 0);
        //     int Ni = data.extent(0);
        //     int Nj = data.extent(1);
        //     int Nk = data.extent(2);
        //     for( int i=0; i<Ni; ++i ) {
        //         for( int j=0; j<Nj; ++j ) {
        //            for( int k=0; k<Nk; ++k ) {
        //               data(i,j,k) = raw[i + Ni*j + Ni*Nj*k];
        //            }
        //         }
        //     }
        // }
        ////////////

        auto status = H5Dread(dataset_id, to_h5_datatype(make_datatype<T>()), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data_handle());
        ATLAS_ASSERT(status >= 0);
    }


    template <typename T>
    void readT(array::Array& array) {
        switch (array.rank()) {
            case 1: return read(array::make_host_view<T,1>(array).as_mdspan());
            case 2: return read(array::make_host_view<T,2>(array).as_mdspan());
            case 3: return read(array::make_host_view<T,3>(array).as_mdspan());
            case 4: return read(array::make_host_view<T,4>(array).as_mdspan());
            default: ATLAS_THROW_EXCEPTION("Unsupported rank " << array.rank());
        }
    }

    void read(array::Array& array) {
        switch (array.datatype().kind()) {
            case DataType::KIND_REAL64 : return readT<double>(array);
            case DataType::KIND_REAL32 : return readT<float>(array);
            case DataType::KIND_INT64  : return readT<long>(array);
            case DataType::KIND_INT32  : return readT<int>(array);
            default: ATLAS_THROW_EXCEPTION("Unsupported datatype " << array.datatype().str());
        }
    }

    const std::string& name() const { return name_; }

    operator ::hid_t() const {
        return dataset_id;
    }
private:

    DataType to_datatype(hid_t datatype_id) const {
        auto type_class = H5Tget_class(datatype_id);
        auto type_size  = H5Tget_size(datatype_id);
        switch (type_class) {
            case H5T_INTEGER:
                if (type_size == 8) {
                    return make_datatype<std::int64_t>();
                }
                else if (type_size == 4) {
                    return make_datatype<std::int32_t>();
                }
                ATLAS_THROW_EXCEPTION("Unsupported type_size");
            case H5T_FLOAT:
                if (type_size == 8) {
                    return make_datatype<double>();
                }
                else if (type_size == 4) {
                    return make_datatype<float>();
                }
                ATLAS_THROW_EXCEPTION("Unsupported type_size");
            case H5T_STRING:
                ATLAS_THROW_EXCEPTION("Unsupported datatype: H5T_STRING");
            case H5T_COMPOUND:
                ATLAS_THROW_EXCEPTION("Unsupported datatype: H5T_COMPOUND");
            case H5T_ARRAY:
                ATLAS_THROW_EXCEPTION("Unsupported datatype: H5T_ARRAY");
            default:
                ATLAS_THROW_EXCEPTION("Unsupported datatype: Other");
        }
    }

    ::hid_t to_h5_datatype(DataType dt) const {
        switch (dt.kind()) {
            case DataType::KIND_INT32:
                return H5T_NATIVE_INT;
            case DataType::KIND_INT64:
                return H5T_NATIVE_LONG;
            case DataType::KIND_REAL32:
                return H5T_NATIVE_FLOAT;
            case DataType::KIND_REAL64:
                return H5T_NATIVE_FLOAT;
            default:
                ATLAS_THROW_EXCEPTION("Unsupported datatype " << dt.str());
        }
    }

    std::string name_;
    ::hid_t dataspace_id{-1};
    ::hid_t dataset_id{-1};
};

class HDF5_File {
public:
    HDF5_File(const std::string& path) {
        file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            Log::error() << "ERROR: Failed to open file '"<<path<<"'"<<std::endl;
        }
    }
    bool is_open() const {
        return file_id >= 0;
    }
    ~HDF5_File() {
        if (is_open()) {
            H5Fclose(file_id);
        }
    }

    HDF5_Dataset dataset(const std::string& dataset_name) {
        ATLAS_ASSERT(is_open());
        return HDF5_Dataset(file_id, dataset_name);
    }

    operator ::hid_t() const {
        return file_id;
    }
    ::hid_t file_id;
};

//-----------------------------------------------------------------------------

} // namespace
