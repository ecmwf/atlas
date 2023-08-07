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

#include "eckit/config/Resource.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/mesh.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace debug {

inline gidx_t global_index(int i = 0) {
    static std::vector<gidx_t> g =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_GLOBAL_INDEX", std::vector<gidx_t>{-1});
    return g[i];
}

inline gidx_t node_global_index(int i = 0) {
    static std::vector<gidx_t> g =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_GLOBAL_INDEX", std::vector<gidx_t>{-1});
    return g[i];
}

inline gidx_t edge_global_index(int i = 0) {
    static std::vector<gidx_t> g =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_EDGE_GLOBAL_INDEX", std::vector<gidx_t>{-1});
    return g[i];
}

inline gidx_t cell_global_index(int i = 0) {
    static std::vector<gidx_t> g =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_CELL_GLOBAL_INDEX", std::vector<gidx_t>{-1});
    return g[i];
}

inline gidx_t node_uid(int i = 0) {
    static std::vector<gidx_t> g =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_UID", std::vector<gidx_t>{-1});
    return g[i];
}

inline bool is_node_global_index(gidx_t x) {
    static std::vector<gidx_t> v =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_GLOBAL_INDEX", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline bool is_global_index(gidx_t x) {
    static std::vector<gidx_t> v =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_GLOBAL_INDEX", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline bool is_edge_global_index(gidx_t x) {
    static std::vector<gidx_t> v =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_EDGE_GLOBAL_INDEX", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline bool is_cell_global_index(gidx_t x) {
    static std::vector<gidx_t> v =
        eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_CELL_GLOBAL_INDEX", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline bool is_node_uid(gidx_t x) {
    static std::vector<gidx_t> v = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_NODE_UID", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline bool is_cell_uid(gidx_t x) {
    static std::vector<gidx_t> v = eckit::Resource<std::vector<gidx_t>>("$ATLAS_DEBUG_CELL_UID", std::vector<gidx_t>());
    for (gidx_t g : v) {
        if (x == g)
            return true;
    }
    return false;
}

inline int mpi_rank(int i = 0) {
    static std::vector<long> g = eckit::Resource<std::vector<long>>("$ATLAS_DEBUG_MPI_RANK", std::vector<long>{-1});
    return g[i];
}

inline int is_mpi_rank() {
    static std::vector<long> v = eckit::Resource<std::vector<long>>("$ATLAS_DEBUG_MPI_RANK", std::vector<long>());
    static int r               = mpi::rank();
    for (long g : v) {
        if (r == g)
            return true;
    }
    return false;
}

inline int is_mpi_rank(int x) {
    static std::vector<long> v = eckit::Resource<std::vector<long>>("$ATLAS_DEBUG_MPI_RANK", std::vector<long>());
    for (long g : v) {
        if (x == g)
            return true;
    }
    return false;
}
inline std::string rank_str() {
    static std::string s = "[" + std::to_string(mpi::rank()) + "] ";
    return s;
}

}  // namespace debug
}  // namespace atlas
