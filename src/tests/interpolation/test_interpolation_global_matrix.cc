/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/linalg/sparse.h"

#include "atlas/functionspace/CellColumns.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas::test {

//-----------------------------------------------------------------------------

Config config_scheme(std::string scheme_str) {
    Config scheme;
    if (scheme_str == "linear") {
        scheme.set("type", "structured-linear2D");
        // The stencil does not require any halo, but we set it to 1 for pole treatment!
        scheme.set("halo", 1);
    }
    if (scheme_str == "cubic") {
        scheme.set("type", "structured-cubic2D");
        scheme.set("halo", 2);
    }
    if (scheme_str == "quasicubic") {
        scheme.set("type", "structured-quasicubic2D");
        scheme.set("halo", 2);
    }
    if (scheme_str == "conservative") {
        scheme.set("type", "conservative-spherical-polygon");
        scheme.set("src_cell_data", false);
        scheme.set("tgt_cell_data", false);
    }
    if (scheme_str == "finite-element") {
        scheme.set("type", "finite-element");
        scheme.set("halo", 1);
    }

    scheme.set("name", scheme_str);
    return scheme;
}


Config create_fspaces(const std::string& scheme_str, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out) {
    Config scheme = config_scheme(scheme_str);
    auto scheme_type = scheme.getString("type");
    if (scheme_type == "finite-element") {
        auto inmesh = Mesh(input_grid);
        auto outmesh = Mesh(output_grid, grid::MatchingPartitioner(inmesh));
        fs_in = NodeColumns(inmesh, scheme);
        fs_out = NodeColumns(outmesh);
    }
    else if (scheme_type != "conservative-spherical-polygon") {
        fs_in = StructuredColumns(input_grid, scheme);
        fs_out = StructuredColumns(output_grid, grid::MatchingPartitioner(fs_in), scheme);
    }
    return scheme;
}


CASE("test_interpolation_structured using grid API") {

using SparseMatrixStorage = atlas::linalg::SparseMatrixStorage;

    auto  assemble_global_matrix = [&](const std::string scheme_str, const Grid& input_grid, const Grid& output_grid, const int mpi_root) {
        Interpolation interpolator;
        FunctionSpace fs_in;
        FunctionSpace fs_out;
        ATLAS_TRACE_SCOPE("Create interpolation") {
            auto scheme = create_fspaces(scheme_str, input_grid, output_grid, fs_in, fs_out);
            if (scheme_str == "conservative") {
                interpolator = Interpolation{ scheme, input_grid, output_grid };
            }
            else {
                interpolator = Interpolation{ scheme, fs_in, fs_out };
            }
        }

        auto tgt_field = interpolator.target().createField<double>();
        ATLAS_TRACE_SCOPE("parallel interpolation [reference]") {
            auto field_in = interpolator.source().createField<double>();
            auto lonlat_in = array::make_view<double,2>(interpolator.source().lonlat());
            auto view_in = array::make_view<double,1>(field_in);
            for(idx_t j=0; j<field_in.size(); ++j) {
                view_in(j) = util::function::vortex_rollup(lonlat_in(j,0),lonlat_in(j,1), 1.);
            }
            interpolator.execute(field_in, tgt_field);
        }

        SparseMatrixStorage matrix;
        ATLAS_TRACE_SCOPE("assemble global matrix") {
            matrix = interpolation::assemble_global_matrix(interpolator, mpi_root);
        }

        // verification via the gathered field
        //
        std::vector<double> src_data(input_grid.size());
        std::vector<double> tgt_data(output_grid.size());
        ATLAS_TRACE_SCOPE("initialize source") {
            idx_t n{0};
            for (auto p : input_grid.lonlat()) {
                src_data[n++] = util::function::vortex_rollup(p.lon(), p.lat(), 1.);
            }
        }
        // gather the global field from the distributed one
        auto tgt_field_global = interpolator.target().createField<double>(option::global(mpi_root));
        tgt_field.haloExchange();
        interpolator.target().gather(tgt_field, tgt_field_global);

        if (mpi::comm().rank() == mpi_root) {
            ATLAS_TRACE_SCOPE("spmv") {
              auto src = eckit::linalg::Vector(src_data.data(), src_data.size());
              auto tgt = eckit::linalg::Vector(tgt_data.data(), tgt_data.size());
              auto eckit_matrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(matrix);
              eckit::linalg::LinearAlgebraSparse::backend().spmv(eckit_matrix, src, tgt);
            }
            auto tfield_global_v = array::make_view<double,1>(tgt_field_global);
            for (gidx_t i = 0; i < tgt_data.size(); ++i) {
                EXPECT_APPROX_EQ(tgt_data[i], tfield_global_v(i), 1.e-14);
            }
        }

        // avoid deadlocks whilst waiting for proc mpi_root
        mpi::comm().barrier();

        return std::make_tuple(fs_in, fs_out, std::move(matrix));
    };

    auto distribute_global_matrix = [&](const FunctionSpace& fs_in, const FunctionSpace& fs_out, const SparseMatrixStorage& gmatrix, int mpi_root) {
        return interpolation::distribute_global_matrix(fs_in, fs_out, gmatrix, mpi_root);
    };

    auto do_assemble_distribute_matrix = [&](const std::string scheme_str, const Grid& input_grid, const Grid& output_grid, const int mpi_root) {
        FunctionSpace fs_in;
        FunctionSpace fs_out;
        SparseMatrixStorage gmatrix;
        std::tie(fs_in, fs_out, gmatrix) = assemble_global_matrix(scheme_str, input_grid, output_grid, mpi_root);
        SparseMatrixStorage matrix = distribute_global_matrix(fs_in, fs_out, gmatrix, mpi_root);
        
        // createinterpolator

    };

    auto test_matrix_assemble_distribute = [&](const Grid& input_grid, const Grid& output_grid) {
        int mpi_root = 0;
        do_assemble_distribute_matrix("linear", input_grid, output_grid, mpi_root);

    return;

        mpi_root = mpi::size() - 1;
        do_assemble_distribute_matrix("linear", input_grid, output_grid, mpi_root);
        do_assemble_distribute_matrix("cubic", input_grid, output_grid, mpi_root);
        do_assemble_distribute_matrix("quasicubic", input_grid, output_grid, mpi_root);
        do_assemble_distribute_matrix("conservative", input_grid, output_grid, mpi_root);
        do_assemble_distribute_matrix("finite-element", input_grid, output_grid, mpi_root);
    };

    test_matrix_assemble_distribute(Grid("O32"), Grid("O64"));
}

}  // namespace


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
