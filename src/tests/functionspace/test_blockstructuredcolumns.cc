/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas;
using namespace atlas::functionspace;
using namespace atlas::test;
using namespace atlas::util;

namespace {
template<class ValueType>
void run_scatter_gather(const Grid& grid, const BlockStructuredColumns& fs, int nlev, int nvar) {
    auto glb_field = fs.createField<ValueType>( option::global() | option::levels(nlev) | option::variables(nvar) );
    // Only allocated on rank 0
    if( atlas::mpi::comm().rank() != 0 ) {
        EXPECT_EQ( glb_field.shape(0), 0);
    }
    else {
        EXPECT_EQ( glb_field.shape(0), grid.size());
    }


    auto glb_field_v = array::make_view<ValueType, 3>(glb_field);

    for (gidx_t p = 0; p < glb_field.shape(0); p++) {
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                glb_field_v(p, jlev, jvar) = p*nvar*nlev + jlev*nvar + jvar;
            }
        }
    }

    auto field = fs.createField(glb_field);

    EXPECT_EQ(glb_field.horizontal_dimension(), (std::vector<idx_t>{0}));
    EXPECT_EQ(field.horizontal_dimension(), (std::vector<idx_t>{0,3}));

    fs.scatter(glb_field, field);

    auto glb_field_2 = fs.createField(glb_field , option::global(atlas::mpi::comm().size()-1));
    // Only allocated on rank MPI_SIZE-1
    if( atlas::mpi::comm().rank() != atlas::mpi::comm().size() - 1 ) {
        EXPECT_EQ( glb_field_2.shape(0), 0);
    }
    else {
        EXPECT_EQ( glb_field_2.shape(0), grid.size());
    }

    fs.gather(field, glb_field_2);

    auto glb_field_2_v = array::make_view<ValueType, 3>(glb_field_2);
    for (gidx_t p = 0; p < glb_field_2.shape(0); p++) {
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                EXPECT_EQ(glb_field_2_v(p, jlev, jvar), p*nvar*nlev + jlev*nvar + jvar);
            }
        }
    }
}

} // namespace

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


CASE("test_BlockStructuredColumns") {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");
    idx_t halo           = eckit::Resource<int>("--halo",   0);
    idx_t nproma         = eckit::Resource<int>("--nproma", 12);
    idx_t nvar           = 3;
    idx_t nlev           = 4;

    auto grid = StructuredGrid(gridname);

    util::Config config;
    config.set("halo", halo);
    config.set("nproma", nproma);
    config.set("levels", nlev);
    auto fs = functionspace::BlockStructuredColumns(grid, config);

    SECTION("blocking indices parallel without halo") {
        ATLAS_DEBUG_VAR(mpi::comm().size());
        ATLAS_DEBUG_VAR(halo);
        ATLAS_DEBUG_VAR(fs.size());
        ATLAS_DEBUG_VAR(fs.nproma());
        ATLAS_DEBUG_VAR(fs.nblks());
        ATLAS_DEBUG_VAR(fs.block(fs.nblks()-1).size());
        ATLAS_DEBUG_VAR(fs.levels());
        ATLAS_DEBUG_VAR(nvar);

        idx_t ngptot    = fs.size();
        idx_t nblk      = ngptot / nproma;
        idx_t nrof_last = nproma;
        if (ngptot % nproma > 0) {
            // Correction when ngptot is not divisible by nproma
            nrof_last = ngptot - nblk * nproma;
            ++nblk;
        }
        EXPECT_EQ(fs.nproma(), nproma);
        EXPECT_EQ(fs.levels(), nlev);
        EXPECT_EQ(fs.nblks(), nblk);
        EXPECT_EQ(fs.block(nblk-1).size(), nrof_last);

        ATLAS_DEBUG("Test cover full iteration space");
        idx_t jpoint = 0;
        for (idx_t jblk = 0; jblk < nblk; ++jblk) {
            auto blk = fs.block(jblk);
            for (idx_t jrof = 0; jrof < blk.size(); ++jrof, ++jpoint) {
                idx_t h_idx  = fs.index(jblk, jrof);
                idx_t hh_idx  = blk.index(jrof);
                EXPECT_EQ(h_idx, jpoint);
                EXPECT_EQ(h_idx, hh_idx);
            }
        }
        EXPECT_EQ(jpoint, ngptot);

        Field field = fs.createField<gidx_t>(option::name("field") | option::variables(nvar));
        EXPECT_EQ(field.shape(0), nblk );
        EXPECT_EQ(field.shape(1), nvar );
        EXPECT_EQ(field.shape(2), nlev );
        EXPECT_EQ(field.shape(3), nproma );

        // When above is failing, we need to really stop here.
        REQUIRE( not current_test().failed() );

        auto value = array::make_view<gidx_t, 4>(field);
        auto g     = array::make_view<gidx_t, 1>(fs.global_index());

        ATLAS_DEBUG("Set");
        for (idx_t jblk = 0; jblk < nblk; ++jblk) {
            auto blk = fs.block(jblk);
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0; jrof < blk.size(); ++jrof) {
                        idx_t h_idx  = blk.index(jrof);
                        value(jblk,jvar,jlev,jrof) = g(h_idx);
                    }
                }
            }
        }

        ATLAS_DEBUG("Check");
        std::vector<idx_t> counters;
        counters.resize(nvar * nlev);
        std::fill(counters.begin(), counters.end(), 0);
        for (idx_t jblk = 0; jblk < nblk; ++jblk) {
            auto blk = fs.block(jblk);
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0; jrof < blk.size(); ++jrof) {
                        idx_t& icheck = counters[jvar*nlev + jlev];
                        EXPECT_EQ(value(jblk,jvar,jlev,jrof), g(icheck));
                        icheck++;
                    }
                }
            }
        }

        auto xy    = array::make_view<double, 2>(fs.xy());
        auto r     = array::make_view<idx_t, 1>(fs.remote_index());
        auto p     = array::make_view<int, 1>(fs.partition());
    }

    SECTION("test_BlockStructuredColumns scatter/gather") {
        auto fs     = functionspace::BlockStructuredColumns(grid, config);
        run_scatter_gather<idx_t>(grid, fs, nlev, nvar);
        run_scatter_gather<gidx_t>(grid, fs, nlev, nvar);
        run_scatter_gather<float>(grid, fs, nlev, nvar);
        run_scatter_gather<double>(grid, fs, nlev, nvar);
    }
}

//-----------------------------------------------------------------------------


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
