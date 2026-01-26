/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "ScripIO.h"

#include <numeric>

#if ATLAS_HAVE_NETCDF
#include <netcdf>
#endif

#include "atlas/linalg/sparse/SparseMatrixToTriplets.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/runtime/Exception.h"


namespace atlas {

using SparseMatrixStorage = linalg::SparseMatrixStorage;
using scrip_index = int;
using scrip_value = double;
using scrip_size  = size_t;

SparseMatrixStorage ScripIO::read(const std::string& matrix_name) {
#if ATLAS_HAVE_NETCDF == 0
    ATLAS_THROW_EXCEPTION("Cannot read SCRIP files: Atlas not compiled with NetCDF support");
#else
    size_t Nr  = 0;
    size_t Nc  = 0;
    size_t nnz = 0;

    std::vector<scrip_index> rows;
    std::vector<scrip_index> cols;
    std::vector<scrip_value> vals;

    // read SCRIP file
    try {
        netCDF::NcFile f(matrix_name, netCDF::NcFile::read);

        Nr  = f.getDim("dst_grid_size").getSize();
        Nc  = f.getDim("src_grid_size").getSize();
        nnz = f.getDim("num_links").getSize();

        ATLAS_ASSERT(Nr > 0);
        ATLAS_ASSERT(Nc > 0);
        ATLAS_ASSERT(nnz > 0);

        auto var_rows = f.getVar("dst_address");
        ATLAS_ASSERT(var_rows.getDimCount() == 1 && var_rows.getDim(0).getSize() == nnz);
        rows.resize(nnz);  // NOTE: not compressed
        var_rows.getVar(rows.data());

        auto var_cols = f.getVar("src_address");
        ATLAS_ASSERT(var_cols.getDimCount() == 1 && var_cols.getDim(0).getSize() == nnz);
        cols.resize(nnz);
        var_cols.getVar(cols.data());

        auto var_vals = f.getVar("remap_matrix");
        ATLAS_ASSERT(var_vals.getDimCount() == 2 && var_vals.getDim(0).getSize() == nnz && var_vals.getDim(1).getSize() == 1);
        vals.resize(nnz);
        var_vals.getVar(vals.data());
    }
    catch (netCDF::exceptions::NcException& e) {
        ATLAS_THROW_EXCEPTION("SCRIP reading failed : " << e.what());
    }
    constexpr scrip_index scrip_base = 1;
    constexpr bool is_sorted = false;
    return atlas::linalg::make_sparse_matrix_storage_from_rows_columns_values(Nr, Nc, rows, cols, vals, scrip_base, is_sorted);

#endif // ATLAS_HAVE_NETCDF
}


void ScripIO::write(const SparseMatrixStorage& matrix, const std::string& matrix_name) {
#if ATLAS_HAVE_NETCDF == 0
    ATLAS_THROW_EXCEPTION("Cannot write SCRIP file: Atlas not compiled with NetCDF support");
#else
    scrip_size Nr  = matrix.rows();
    scrip_size Nc  = matrix.cols();
    scrip_size nnz = matrix.nnz();

    std::vector<scrip_index> ia(nnz);
    std::vector<scrip_index> ja(nnz);
    std::vector<scrip_value> a(nnz);

    atlas::linalg::sparse_matrix_to_rows_columns_values(matrix, ia, ja, a);
    constexpr scrip_index scrip_index_base = 1;
    for (size_t i = 0; i < matrix.nnz(); ++i) {
        ja[i] += scrip_index_base;
        ia[i] += scrip_index_base;
    }

    try {
        netCDF::NcFile f(matrix_name, netCDF::NcFile::replace);

        f.addDim("dst_grid_size", Nr);
        f.addDim("src_grid_size", Nc);

        std::vector<netCDF::NcDim> dims;
        dims.push_back(f.addDim("num_links", nnz));
        netCDF::NcVar nc_dstaddr = f.addVar("dst_address", netCDF::ncInt, dims);
        netCDF::NcVar nc_srcaddr = f.addVar("src_address", netCDF::ncInt, dims);
        dims.push_back(f.addDim("num_wgts", 1));
        netCDF::NcVar nc_rmatrix = f.addVar("remap_matrix", netCDF::ncDouble, dims);

        nc_dstaddr.putVar(ia.data());
        nc_srcaddr.putVar(ja.data());
        nc_rmatrix.putVar(a.data());
        f.close();
    }
    catch (netCDF::exceptions::NcException& e) {
        ATLAS_THROW_EXCEPTION("SCRIP writing failed : " << e.what());
    }
#endif
}

}
