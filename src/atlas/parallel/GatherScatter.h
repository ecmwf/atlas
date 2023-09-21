/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "atlas/array/ArrayView.h"
#include "atlas/library/config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"

namespace atlas {
namespace parallel {

//----------------------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
class Field {
public:
    Field() {}

    Field(const DATA_TYPE data_[], const idx_t var_strides_[], const idx_t var_shape_[], const idx_t var_rank_):
        data(const_cast<DATA_TYPE*>(data_)), var_rank(var_rank_) {
        var_strides.assign(var_strides_, var_strides_ + var_rank_);
        var_shape.assign(var_shape_, var_shape_ + var_rank_);
    }

    Field(const DATA_TYPE data_[], const idx_t nb_vars): data(const_cast<DATA_TYPE*>(data_)), var_rank(1) {
        var_strides.assign(1, 1);
        var_shape.assign(1, nb_vars);
    }

    template <typename Value, int Rank>
    Field(const array::ArrayView<Value, Rank>& arr) {
        data = const_cast<DATA_TYPE*>(arr.data());

        var_rank = Rank;
        var_strides.resize(var_rank);
        var_shape.resize(var_rank);

        if (arr.rank() > 1) {
            var_strides[0] = arr.stride(0);
            var_shape[0]   = 1;
            for (int j = 1; j < Rank; ++j) {
                var_strides[j] = arr.stride(j);
                var_shape[j]   = arr.shape(j);
            }
        }
        else {
            var_strides[0] = arr.stride(0);
            var_shape[0]   = 1;
        }
    }

    template <typename Value, int Rank>
    Field(const array::LocalView<Value, Rank>& arr) {
        data = const_cast<DATA_TYPE*>(arr.data());

        var_rank = Rank;
        var_strides.resize(var_rank);
        var_shape.resize(var_rank);

        if (arr.rank() > 1) {
            var_strides[0] = arr.stride(0);
            var_shape[0]   = 1;
            for (int j = 1; j < Rank; ++j) {
                var_strides[j] = arr.stride(j);
                var_shape[j]   = arr.shape(j);
            }
        }
        else {
            var_strides[0] = arr.stride(0);
            var_shape[0]   = 1;
        }
    }

public:
    DATA_TYPE* data;
    std::vector<idx_t> var_strides;
    std::vector<idx_t> var_shape;
    idx_t var_rank;
};

class GatherScatter : public util::Object {
public:
    GatherScatter();
    GatherScatter(const std::string& name);
    virtual ~GatherScatter() {}

public:  // methods
    std::string name() const { return name_; }

    /// @brief Setup
    /// @param [in] mpi_comm     Name of mpi communicator
    /// @param [in] part         List of partitions
    /// @param [in] remote_idx   List of local indices on remote partitions
    /// @param [in] base         values of remote_idx start at "base"
    /// @param [in] glb_idx      List of global indices
    /// @param [in] parsize      size of given lists
    void setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[], const idx_t parsize);

    /// @brief Setup
    /// @param [in] part         List of partitions
    /// @param [in] remote_idx   List of local indices on remote partitions
    /// @param [in] base         values of remote_idx start at "base"
    /// @param [in] glb_idx      List of global indices
    /// @param [in] parsize      size of given lists
    void setup(const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[], const idx_t parsize);

    /// @brief Setup
    /// @param [in] mpi_comm     Name of mpi communicator
    /// @param [in] part         List of partitions
    /// @param [in] remote_idx   List of local indices on remote partitions
    /// @param [in] base         values of remote_idx start at "base"
    /// @param [in] glb_idx      List of global indices
    /// @param [in] parsize      size of given lists
    /// @param [in] mask         Mask indices not to include in the communication
    ///                          pattern (0=include,1=exclude)
    void setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[], const int mask[],
               const idx_t parsize);

    /// @brief Setup
    /// @param [in] part         List of partitions
    /// @param [in] remote_idx   List of local indices on remote partitions
    /// @param [in] base         values of remote_idx start at "base"
    /// @param [in] glb_idx      List of global indices
    /// @param [in] parsize      size of given lists
    /// @param [in] mask         Mask indices not to include in the communication
    ///                          pattern (0=include,1=exclude)
    void setup(const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[], const int mask[],
               const idx_t parsize);

    template <typename DATA_TYPE>
    void gather(const DATA_TYPE ldata[], const idx_t lstrides[], const idx_t lshape[], const idx_t lrank,
                const idx_t lmpl_idxpos[], const idx_t lmpl_rank, DATA_TYPE gdata[], const idx_t gstrides[],
                const idx_t gshape[], const idx_t grank, const idx_t gmpl_idxpos[], const idx_t gmpl_rank,
                const idx_t root) const;

    template <typename DATA_TYPE>
    void gather(const DATA_TYPE ldata[], const idx_t lvar_strides[], const idx_t lvar_shape[], const idx_t lvar_rank,
                DATA_TYPE gdata[], const idx_t gvar_strides[], const idx_t gvar_shape[], const idx_t gvar_rank,
                const idx_t root = 0) const;

    template <typename DATA_TYPE>
    void gather(parallel::Field<DATA_TYPE const> lfields[], parallel::Field<DATA_TYPE> gfields[], const idx_t nb_fields,
                const idx_t root = 0) const;

    template <typename DATA_TYPE, int LRANK, int GRANK>
    void gather(const array::ArrayView<DATA_TYPE, LRANK>& ldata, array::ArrayView<DATA_TYPE, GRANK>& gdata,
                const idx_t root = 0) const;

    template <typename DATA_TYPE>
    void scatter(parallel::Field<DATA_TYPE const> gfields[], parallel::Field<DATA_TYPE> lfields[],
                 const idx_t nb_fields, const idx_t root = 0) const;

    template <typename DATA_TYPE>
    void scatter(const DATA_TYPE gdata[], const idx_t gstrides[], const idx_t gshape[], const idx_t grank,
                 const idx_t gmpl_idxpos[], const idx_t gmpl_rank, DATA_TYPE ldata[], const idx_t lstrides[],
                 const idx_t lshape[], const idx_t lrank, const idx_t lmpl_idxpos[], const idx_t lmpl_rank,
                 const idx_t root) const;

    template <typename DATA_TYPE>
    void scatter(const DATA_TYPE gdata[], const idx_t gvar_strides[], const idx_t gvar_shape[], const idx_t gvar_rank,
                 DATA_TYPE ldata[], const idx_t lvar_strides[], const idx_t lvar_shape[], const idx_t lvar_rank,
                 const idx_t root = 0) const;

    template <typename DATA_TYPE, int GRANK, int LRANK>
    void scatter(const array::ArrayView<DATA_TYPE, GRANK>& gdata, array::ArrayView<DATA_TYPE, LRANK>& ldata,
                 const idx_t root = 0) const;

    gidx_t glb_dof() const { return glbcnt_; }

    idx_t loc_dof() const { return loccnt_; }

    const mpi::Comm& comm() const { return *comm_; }

private:  // methods
    template <typename DATA_TYPE>
    void pack_send_buffer(const parallel::Field<DATA_TYPE const>& field, const std::vector<int>& sendmap,
                          DATA_TYPE send_buffer[]) const;

    template <typename DATA_TYPE>
    void unpack_recv_buffer(const std::vector<int>& recvmap, const DATA_TYPE recv_buffer[],
                            const parallel::Field<DATA_TYPE>& field) const;

    template <typename DATA_TYPE, int RANK>
    void var_info(const array::ArrayView<DATA_TYPE, RANK>& arr, std::vector<idx_t>& varstrides,
                  std::vector<idx_t>& varshape) const;

private:  // data
    std::string name_;
    int loccnt_;
    int glbcnt_;
    std::vector<int> glbcounts_;
    std::vector<int> glbdispls_;
    std::vector<int> locmap_;
    std::vector<int> glbmap_;

    const mpi::Comm* comm_;
    idx_t nproc;
    idx_t myproc;

    bool is_setup_;

    idx_t parsize_;
    friend class Checksum;

    size_t glb_cnt(idx_t root) const { return myproc == root ? glbcnt_ : 0; }
};

//--------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
void GatherScatter::gather(parallel::Field<DATA_TYPE const> lfields[], parallel::Field<DATA_TYPE> gfields[],
                           idx_t nb_fields, const idx_t root) const {
    if (!is_setup_) {
        throw_Exception("GatherScatter was not setup", Here());
    }

    for (idx_t jfield = 0; jfield < nb_fields; ++jfield) {
        const idx_t lvar_size =
            std::accumulate(lfields[jfield].var_shape.data(),
                            lfields[jfield].var_shape.data() + lfields[jfield].var_rank, 1, std::multiplies<idx_t>());
        const idx_t gvar_size =
            std::accumulate(gfields[jfield].var_shape.data(),
                            gfields[jfield].var_shape.data() + gfields[jfield].var_rank, 1, std::multiplies<idx_t>());
        const size_t loc_size = loccnt_ * lvar_size;
        const size_t glb_size = glb_cnt(root) * gvar_size;
        std::vector<DATA_TYPE> loc_buffer(loc_size);
        std::vector<DATA_TYPE> glb_buffer(glb_size);
        std::vector<int> glb_displs(nproc);
        std::vector<int> glb_counts(nproc);

        for (idx_t jproc = 0; jproc < nproc; ++jproc) {
            glb_counts[jproc] = glbcounts_[jproc] * gvar_size;
            glb_displs[jproc] = glbdispls_[jproc] * gvar_size;
        }

        /// Pack

        pack_send_buffer(lfields[jfield], locmap_, loc_buffer.data());

        /// Gather

        ATLAS_TRACE_MPI(GATHER) {
            comm().gatherv(loc_buffer, glb_buffer, glb_counts, glb_displs, root);
        }

        /// Unpack
        if (myproc == root)
            unpack_recv_buffer(glbmap_, glb_buffer.data(), gfields[jfield]);
    }
}

template <typename DATA_TYPE>
void GatherScatter::gather(const DATA_TYPE ldata[], const idx_t lvar_strides[], const idx_t lvar_shape[],
                           const idx_t lvar_rank, DATA_TYPE gdata[], const idx_t gvar_strides[],
                           const idx_t gvar_shape[], const idx_t gvar_rank, const idx_t root) const {
    parallel::Field<DATA_TYPE const> lfield(ldata, lvar_strides, lvar_shape, lvar_rank);
    parallel::Field<DATA_TYPE> gfield(gdata, gvar_strides, gvar_shape, gvar_rank);
    gather(&lfield, &gfield, 1, root);
}

template <typename DATA_TYPE>
void GatherScatter::scatter(parallel::Field<DATA_TYPE const> gfields[], parallel::Field<DATA_TYPE> lfields[],
                            const idx_t nb_fields, const idx_t root) const {
    if (!is_setup_) {
        throw_Exception("GatherScatter was not setup", Here());
    }

    for (idx_t jfield = 0; jfield < nb_fields; ++jfield) {
        const int lvar_size =
            std::accumulate(lfields[jfield].var_shape.data(),
                            lfields[jfield].var_shape.data() + lfields[jfield].var_rank, 1, std::multiplies<int>());
        const int gvar_size =
            std::accumulate(gfields[jfield].var_shape.data(),
                            gfields[jfield].var_shape.data() + gfields[jfield].var_rank, 1, std::multiplies<int>());
        const size_t loc_size = loccnt_ * lvar_size;
        const size_t glb_size = glb_cnt(root) * gvar_size;
        std::vector<DATA_TYPE> loc_buffer(loc_size);
        std::vector<DATA_TYPE> glb_buffer(glb_size);
        std::vector<int> glb_displs(nproc);
        std::vector<int> glb_counts(nproc);

        for (idx_t jproc = 0; jproc < nproc; ++jproc) {
            glb_counts[jproc] = glbcounts_[jproc] * gvar_size;
            glb_displs[jproc] = glbdispls_[jproc] * gvar_size;
        }

        /// Pack
        if (myproc == root)
            pack_send_buffer(gfields[jfield], glbmap_, glb_buffer.data());

        /// Scatter

        ATLAS_TRACE_MPI(SCATTER) {
            comm().scatterv(glb_buffer.begin(), glb_buffer.end(), glb_counts, glb_displs, loc_buffer.begin(),
                            loc_buffer.end(), root);
        }

        /// Unpack
        unpack_recv_buffer(locmap_, loc_buffer.data(), lfields[jfield]);
    }
}

template <typename DATA_TYPE>
void GatherScatter::scatter(const DATA_TYPE gdata[], const idx_t gvar_strides[], const idx_t gvar_shape[],
                            const idx_t gvar_rank, DATA_TYPE ldata[], const idx_t lvar_strides[],
                            const idx_t lvar_shape[], const idx_t lvar_rank, const idx_t root) const {
    parallel::Field<DATA_TYPE const> gfield(gdata, gvar_strides, gvar_shape, gvar_rank);
    parallel::Field<DATA_TYPE> lfield(ldata, lvar_strides, lvar_shape, lvar_rank);
    scatter(&gfield, &lfield, 1, root);
}

template <typename DATA_TYPE>
void GatherScatter::pack_send_buffer(const parallel::Field<DATA_TYPE const>& field, const std::vector<int>& sendmap,
                                     DATA_TYPE send_buffer[]) const {
    const idx_t sendcnt = static_cast<idx_t>(sendmap.size());

    idx_t ibuf              = 0;
    const idx_t send_stride = field.var_strides[0] * field.var_shape[0];

    switch (field.var_rank) {
        case 1:
            for (idx_t p = 0; p < sendcnt; ++p) {
                const idx_t pp = send_stride * sendmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    DATA_TYPE tmp       = field.data[pp + i * field.var_strides[0]];
                    send_buffer[ibuf++] = tmp;
                }
            }
            break;
        case 2:
            for (idx_t p = 0; p < sendcnt; ++p) {
                const idx_t pp = send_stride * sendmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    const idx_t ii = pp + i * field.var_strides[0];
                    for (idx_t j = 0; j < field.var_shape[1]; ++j) {
                        send_buffer[ibuf++] = field.data[ii + j * field.var_strides[1]];
                    }
                }
            }
            break;
        case 3:
            for (idx_t p = 0; p < sendcnt; ++p) {
                const idx_t pp = send_stride * sendmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    const idx_t ii = pp + i * field.var_strides[0];
                    for (idx_t j = 0; j < field.var_shape[1]; ++j) {
                        const idx_t jj = ii + j * field.var_strides[1];
                        for (idx_t k = 0; k < field.var_shape[2]; ++k) {
                            send_buffer[ibuf++] = field.data[jj + k * field.var_strides[2]];
                        }
                    }
                }
            }
            break;
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

template <typename DATA_TYPE>
void GatherScatter::unpack_recv_buffer(const std::vector<int>& recvmap, const DATA_TYPE recv_buffer[],
                                       const parallel::Field<DATA_TYPE>& field) const {
    const idx_t recvcnt = static_cast<idx_t>(recvmap.size());

    size_t ibuf             = 0;
    const idx_t recv_stride = field.var_strides[0] * field.var_shape[0];

    switch (field.var_rank) {
        case 1:
            for (idx_t p = 0; p < recvcnt; ++p) {
                const idx_t pp = recv_stride * recvmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    field.data[pp + i * field.var_strides[0]] = recv_buffer[ibuf++];
                }
            }
            break;
        case 2:
            for (idx_t p = 0; p < recvcnt; ++p) {
                const idx_t pp = recv_stride * recvmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    const idx_t ii = pp + i * field.var_strides[0];
                    for (idx_t j = 0; j < field.var_shape[1]; ++j) {
                        field.data[ii + j * field.var_strides[1]] = recv_buffer[ibuf++];
                    }
                }
            }
            break;
        case 3:
            for (idx_t p = 0; p < recvcnt; ++p) {
                const idx_t pp = recv_stride * recvmap[p];
                for (idx_t i = 0; i < field.var_shape[0]; ++i) {
                    const idx_t ii = pp + i * field.var_strides[0];
                    for (idx_t j = 0; j < field.var_shape[1]; ++j) {
                        const idx_t jj = ii + j * field.var_strides[1];
                        for (idx_t k = 0; k < field.var_shape[2]; ++k) {
                            field.data[jj + k * field.var_strides[2]] = recv_buffer[ibuf++];
                        }
                    }
                }
            }
            break;
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

template <typename DATA_TYPE, int RANK>
void GatherScatter::var_info(const array::ArrayView<DATA_TYPE, RANK>& arr, std::vector<idx_t>& varstrides,
                             std::vector<idx_t>& varshape) const {
    int rank = std::max(1, RANK - 1);
    varstrides.resize(rank);
    varshape.resize(rank);

    if (RANK > 1) {
        idx_t stride = 1;
        for (int j = RANK - 1; j > 0; --j) {
            varstrides[j - 1] = stride;
            varshape[j - 1]   = arr.shape(j);
            stride *= varshape[j - 1];
        }
        //    varstrides.assign(arr.strides()+1,arr.strides()+RANK);
        //    varshape.assign(arr.shape()+1,arr.shape()+RANK);
    }
    else {
        varstrides[0] = 1;
        varshape[0]   = 1;
    }
}

template <typename DATA_TYPE, int LRANK, int GRANK>
void GatherScatter::gather(const array::ArrayView<DATA_TYPE, LRANK>& ldata, array::ArrayView<DATA_TYPE, GRANK>& gdata,
                           const idx_t root) const {
    if (ldata.shape(0) == parsize_ && gdata.shape(0) == idx_t(glb_cnt(root))) {
        std::vector<parallel::Field<DATA_TYPE const>> lfields(1, parallel::Field<DATA_TYPE const>(ldata));
        std::vector<parallel::Field<DATA_TYPE>> gfields(1, parallel::Field<DATA_TYPE>(gdata));
        gather(lfields.data(), gfields.data(), 1, root);
    }
    else {
        ATLAS_NOTIMPLEMENTED;  // Need to implement with parallel ranks > 1
    }
}

template <typename DATA_TYPE, int GRANK, int LRANK>
void GatherScatter::scatter(const array::ArrayView<DATA_TYPE, GRANK>& gdata, array::ArrayView<DATA_TYPE, LRANK>& ldata,
                            const idx_t root) const {
    if (ldata.shape(0) == parsize_ && gdata.shape(0) == idx_t(glb_cnt(root))) {
        std::vector<parallel::Field<DATA_TYPE const>> gfields(1, parallel::Field<DATA_TYPE const>(gdata));
        std::vector<parallel::Field<DATA_TYPE>> lfields(1, parallel::Field<DATA_TYPE>(ldata));
        scatter(gfields.data(), lfields.data(), 1, root);
    }
    else {
        ATLAS_NOTIMPLEMENTED;  // Need to implement with parallel ranks > 1
    }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {
GatherScatter* atlas__GatherScatter__new();
void atlas__GatherScatter__delete(GatherScatter* This);
void atlas__GatherScatter__setup32(GatherScatter* This, int part[], idx_t remote_idx[], int base, int glb_idx[],
                                   int parsize);
void atlas__GatherScatter__setup64(GatherScatter* This, int part[], idx_t remote_idx[], int base, long glb_idx[],
                                   int parsize);
int atlas__GatherScatter__glb_dof(GatherScatter* This);
void atlas__GatherScatter__gather_int(GatherScatter* This, int ldata[], int lvar_strides[], int lvar_shape[],
                                      int lvar_rank, int gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank);
void atlas__GatherScatter__gather_long(GatherScatter* This, long ldata[], int lvar_strides[], int lvar_shape[],
                                       int lvar_rank, long gdata[], int gvar_strides[], int gvar_shape[],
                                       int gvar_rank);
void atlas__GatherScatter__gather_float(GatherScatter* This, float ldata[], int lvar_strides[], int lvar_shape[],
                                        int lvar_rank, float gdata[], int gvar_strides[], int gvar_shape[],
                                        int gvar_rank);
void atlas__GatherScatter__gather_double(GatherScatter* This, double ldata[], int lvar_strides[], int lvar_shape[],
                                         int lvar_rank, double gdata[], int gvar_strides[], int gvar_shape[],
                                         int gvar_rank);
void atlas__GatherScatter__scatter_int(GatherScatter* This, int gdata[], int gvar_strides[], int gvar_shape[],
                                       int gvar_rank, int ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank);
void atlas__GatherScatter__scatter_long(GatherScatter* This, long gdata[], int gvar_strides[], int gvar_shape[],
                                        int gvar_rank, long ldata[], int lvar_strides[], int lvar_shape[],
                                        int lvar_rank);
void atlas__GatherScatter__scatter_float(GatherScatter* This, float gdata[], int gvar_strides[], int gvar_shape[],
                                         int gvar_rank, float ldata[], int lvar_strides[], int lvar_shape[],
                                         int lvar_rank);
void atlas__GatherScatter__scatter_double(GatherScatter* This, double gdata[], int gvar_strides[], int gvar_shape[],
                                          int gvar_rank, double ldata[], int lvar_strides[], int lvar_shape[],
                                          int lvar_rank);
}

//--------------------------------------------------------------------------------------------------

// typedef GatherScatter Gather;

}  // namespace parallel
}  // namespace atlas
