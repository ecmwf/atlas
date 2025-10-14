/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <numeric>
#include <map>
#include <memory>

#include "hdf5.h"

#include "pluto/pluto.h"

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/field.h"
#include "infer_grid.h"

namespace atlas {

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

    operator ::hid_t() const {
        return file_id;
    }
    ::hid_t file_id;
};

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

    HDF5_Dataset(std::string file_name, const std::string& dataset_name) {
        file_ = std::make_unique<HDF5_File>(file_name);
        if (file_->is_open()) {
            setup(*file_, dataset_name);
        }
    }

    HDF5_Dataset(::hid_t file_id, const std::string& dataset_name) {
        setup(file_id, dataset_name);
    }

private:
    std::unique_ptr<HDF5_File> file_;

    void setup(::hid_t file_id, const std::string& dataset_name) {
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
public:

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
        ATLAS_TRACE("read(mdspan)");
        if (data.size() == 0) {
            return;
        }
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

        if (contiguous) {
            ATLAS_TRACE("H5Dread");
            auto status = H5Dread(dataset_id, to_h5_datatype(make_datatype<T>()), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data_handle());
            ATLAS_ASSERT(status >= 0);
        }
        else {
            pluto::allocator<T> allocator;
            T* raw;
            ATLAS_TRACE_SCOPE("allocate temporary buffer") {
                raw = allocator.allocate(size());
            }
            herr_t status;
            ATLAS_TRACE_SCOPE("H5Dread into temporary buffer") {
                status = H5Dread(dataset_id, to_h5_datatype(make_datatype<T>()), H5S_ALL, H5S_ALL, H5P_DEFAULT, raw);
            }
            ATLAS_TRACE_SCOPE("Copy temporary buffer to mdspan") {
                ATLAS_ASSERT(status >= 0);
                auto extent = extents();
                if constexpr (data.rank() == 1) {
                    ATLAS_ASSERT( extent[0] == data.extent(0));
                    size_t Ni = data.extent(0);
                    for( size_t i=0; i<Ni; ++i ) {
                        data(i) = raw[i];
                    }
                }
                else if constexpr (data.rank() == 2) {
                    ATLAS_ASSERT( extent[0] == data.extent(0));
                    ATLAS_ASSERT( extent[1] == data.extent(1));
                    size_t Ni = data.extent(0);
                    size_t Nj = data.extent(1);
                    size_t stride_i = Nj;
                    for( size_t i=0; i<Ni; ++i ) {
                        for( size_t j=0; j<Nj; ++j ) {
                            data(i,j) = raw[i*stride_i + j];
                        }
                    }
                }
                else if constexpr (data.rank() == 3) {
                    ATLAS_ASSERT( extent[0] == data.extent(0));
                    ATLAS_ASSERT( extent[1] == data.extent(1));
                    ATLAS_ASSERT( extent[2] == data.extent(2));
                    size_t Ni = data.extent(0);
                    size_t Nj = data.extent(1);
                    size_t Nk = data.extent(2);
                    size_t stride_i = Nj*Nk;
                    size_t stride_j = Nk;
                    for( size_t i=0; i<Ni; ++i ) {
                        for( size_t j=0; j<Nj; ++j ) {
                            for( size_t k=0; k<Nk; ++k ) {
                                data(i,j,k) = raw[i*stride_i + j*stride_j + k];
                            }
                        }
                    }
                }
                else if constexpr (data.rank() == 4) {
                    ATLAS_ASSERT( extent[0] == data.extent(0));
                    ATLAS_ASSERT( extent[1] == data.extent(1));
                    ATLAS_ASSERT( extent[2] == data.extent(2));
                    ATLAS_ASSERT( extent[4] == data.extent(4));
                    size_t Ni = data.extent(0);
                    size_t Nj = data.extent(1);
                    size_t Nk = data.extent(2);
                    size_t Nl = data.extent(3);
                    size_t stride_i = Nj*Nk*Nl;
                    size_t stride_j = Nk*Nl;
                    size_t stride_k = Nl;
                    for( size_t i=0; i<Ni; ++i ) {
                        for( size_t j=0; j<Nj; ++j ) {
                            for( size_t k=0; k<Nk; ++k ) {
                                for( size_t l=0; l<Nl; ++l ) {
                                    data(i,j,k,l) = raw[i*stride_i + j*stride_j + k*stride_k + l];
                                }
                            }
                        }
                    }
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
            }
            ATLAS_TRACE_SCOPE("deallocate temporary buffer");
            allocator.deallocate(raw, size());
        }
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
            case DataType::KIND_INT32:  return H5T_NATIVE_INT;
            case DataType::KIND_INT64:  return H5T_NATIVE_LONG;
            case DataType::KIND_REAL32: return H5T_NATIVE_FLOAT;
            case DataType::KIND_REAL64: return H5T_NATIVE_DOUBLE;
            default:
                ATLAS_THROW_EXCEPTION("Unsupported datatype " << dt.str());
        }
    }

    std::string name_;
    ::hid_t dataspace_id{-1};
    ::hid_t dataset_id{-1};
};


class HDF5Reader {
public:
  HDF5Reader(const std::string& file, const std::string& dataset) : file_{file}, dataset_{dataset} {}
  util::Metadata read_metadata() {
    util::Metadata metadata;
    int error = 0;
    int mpi_root = 0;
    if (mpi::rank() == mpi_root) {
        auto& h5_dataset = hdf5_dataset();
        if (!h5_dataset.is_open()) {
            error = 1;
            goto endif;
        }
        int nproma = h5_dataset.extents()[h5_dataset.rank()-1];
        int nblk  = h5_dataset.extents()[0];
        int size = h5_dataset.size();
        int nfld;
        if (size == 0) {
            nfld = 0;
        }
        else if (h5_dataset.rank() == 2) {
            nfld = 1;
        }
        else if (h5_dataset.rank() == 3) {
            nfld = int(h5_dataset.extents()[1]);
        }
        else if (h5_dataset.rank() == 4) {
            nfld = int(h5_dataset.extents()[1] * h5_dataset.extents()[2]);
        }
        else if (h5_dataset.rank() == 5) {
            nfld = int(h5_dataset.extents()[1] * h5_dataset.extents()[2] * h5_dataset.extents()[3]);
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
        // int nflds  = h5_dataset.rank() == 3 ? h5_dataset.extents()[1] : 1;
        auto grid = infer_grid_from_nproma_and_nblks(nproma, nblk);
        if (grid.empty()) {
            Log::error() << "ERROR: Could not detect grid" << std::endl;
            error = 1;
            goto endif;
        }
        metadata.set("name", h5_dataset.name());
        metadata.set("datatype", h5_dataset.datatype().str());
        metadata.set("nblk",  nblk);
        metadata.set("nproma", nproma);
        metadata.set("nfld", nfld);
        metadata.set("size", size);
        metadata.set("grid", grid);
    } endif:
    mpi::comm().broadcast(error, mpi_root);
    if (error) {
        ATLAS_THROW_EXCEPTION("Could not extract metadata from HDF5(file="<<file_<<", dataset="<<dataset_<<")");
    }
    metadata.broadcast(mpi_root);
    return metadata;
  }

  void read_field(Field& field) {
    ATLAS_TRACE();
    if (field.size() == 0) {
        return;
    }
    int mpi_root = 0;
    functionspace::BlockStructuredColumns fs{field.functionspace()};
    pluto::scope::push();
    pluto::host::set_default_resource(pluto::host_pool_resource());
    Field field_glb;
    ATLAS_TRACE_SCOPE("allocate field_glb") {
        field_glb = fs.createField(option::name(field.name()) | option::levels(field.levels()) | option::datatype(field.datatype()) | option::global(mpi_root));
    }
    int error = 0;
    if (mpi::rank() == mpi_root) {
        // Get the dataset from the file
        auto& h5_dataset = hdf5_dataset();
        if (!h5_dataset.is_open()) {
            error = 1;
            goto endif;
        }

        Field dataset_field(h5_dataset.name(), h5_dataset.datatype(), h5_dataset.extents());
        h5_dataset.read(dataset_field);

        switch( dataset_field.datatype().kind() ) {
            case DataType::KIND_REAL64 : {
                if (dataset_field.rank() == 2) {
                    auto glb_view = array::make_view<double,1>(field_glb);
                    auto dataset_view =array::make_view<double,2>(dataset_field);
                    idx_t n = 0;
                    for (idx_t jblk=0; jblk < dataset_view.shape(0); ++jblk) {
                        for (idx_t jrof=0; jrof < dataset_view.shape(1); ++jrof, ++n) {
                            if (n < glb_view.shape(0)) {
                                glb_view(n) = dataset_view(jblk, jrof);
                            }
                        }
                    }
                }
                else if (dataset_field.rank() == 3) {
                    auto glb_view = array::make_view<double,2>(field_glb);
                    auto dataset_view =array::make_view<double,3>(dataset_field);
                    idx_t n = 0;
                    for (idx_t jblk=0; jblk < dataset_view.shape(0); ++jblk) {
                        for (idx_t jrof=0; jrof < dataset_view.shape(2); ++jrof, ++n) {
                            if (n < glb_view.shape(0)) {
                                for (idx_t jfld=0; jfld < glb_view.shape(1); ++jfld) {
                                    glb_view(n, jfld) = dataset_view(jblk, jfld, jrof);
                                }
                            }
                        }
                    }
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
                break;
            }
            default:
                ATLAS_THROW_EXCEPTION("Datatype " << dataset_field.datatype().str() << " not implemented");
        }
    } endif:
    mpi::comm().broadcast(error, mpi_root);
    if (error) {
        ATLAS_THROW_EXCEPTION("Could not extract field data from HDF5(file="<<file_<<", dataset="<<dataset_<<")");
    }
    ATLAS_TRACE_SCOPE("scatter") {
        fs.scatter(field_glb, field);
    }
    pluto::scope::pop();
}


  HDF5_Dataset& hdf5_dataset() {
    if (not hdf5_dataset_) {
        hdf5_dataset_ = std::make_unique<HDF5_Dataset>(file_, dataset_);
    }
    return *hdf5_dataset_;
  }
private:
  std::string file_;
  std::string dataset_;
  std::unique_ptr<HDF5_Dataset> hdf5_dataset_;
};

//-----------------------------------------------------------------------------

} // namespace
