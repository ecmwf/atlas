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

#include <functional>
#include <type_traits>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/StructuredColumns.h"

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class StructuredColumns : public FunctionSpace {
public:
    StructuredColumns();
    StructuredColumns(const FunctionSpace&);
    StructuredColumns(const Grid&, const eckit::Configuration& = util::NoConfig());
    StructuredColumns(const Grid&, const grid::Partitioner&, const eckit::Configuration& = util::NoConfig());
    StructuredColumns(const Grid&, const grid::Distribution&, const eckit::Configuration& = util::NoConfig());
    StructuredColumns(const Grid&, const Vertical&, const eckit::Configuration& = util::NoConfig());
    StructuredColumns(const Grid&, const Vertical&, const grid::Partitioner&,
                      const eckit::Configuration& = util::NoConfig());
    StructuredColumns(const Grid&, const grid::Distribution&, const Vertical&,
                      const eckit::Configuration& = util::NoConfig());

    static std::string type() { return detail::StructuredColumns::static_type(); }

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    idx_t size() const { return functionspace_->size(); }
    idx_t sizeOwned() const { return functionspace_->sizeOwned(); }
    idx_t sizeHalo() const { return functionspace_->sizeHalo(); }

    idx_t levels() const { return functionspace_->levels(); }

    idx_t halo() const { return functionspace_->halo(); }

    const Vertical& vertical() const { return functionspace_->vertical(); }

    const StructuredGrid& grid() const { return functionspace_->grid(); }

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;

    idx_t index(idx_t i, idx_t j) const { return functionspace_->index(i, j); }

    idx_t i_begin(idx_t j) const { return functionspace_->i_begin(j); }
    idx_t i_end(idx_t j) const { return functionspace_->i_end(j); }

    idx_t i_begin_halo(idx_t j) const { return functionspace_->i_begin_halo(j); }
    idx_t i_end_halo(idx_t j) const { return functionspace_->i_end_halo(j); }

    idx_t j_begin() const { return functionspace_->j_begin(); }
    idx_t j_end() const { return functionspace_->j_end(); }

    idx_t j_begin_halo() const { return functionspace_->j_begin_halo(); }
    idx_t j_end_halo() const { return functionspace_->j_end_halo(); }

    idx_t k_begin() const { return functionspace_->k_begin(); }
    idx_t k_end() const { return functionspace_->k_end(); }

    Field xy() const { return functionspace_->xy(); }
    Field partition() const { return functionspace_->partition(); }
    Field global_index() const { return functionspace_->global_index(); }
    Field remote_index() const { return functionspace_->remote_index(); }
    Field index_i() const { return functionspace_->index_i(); }
    Field index_j() const { return functionspace_->index_j(); }
    Field ghost() const { return functionspace_->ghost(); }

    void compute_xy(idx_t i, idx_t j, PointXY& xy) const { return functionspace_->compute_xy(i, j, xy); }
    PointXY compute_xy(idx_t i, idx_t j) const { return functionspace_->compute_xy(i, j); }

    size_t footprint() const { return functionspace_->footprint(); }

    const util::PartitionPolygon& polygon(idx_t halo = 0) const { return functionspace_->polygon(halo); }

    class For {
    public:
        For(const StructuredColumns& fs, const util::Config& config = util::NoConfig()):
            fs_{fs},
            global{config.getBool("global", false)},
            owner{config.getInt("owner", 0)},
            levels{config.getInt("levels", fs_.levels())} {}

    protected:
        const StructuredColumns& fs_;
        bool global;
        idx_t owner;
        idx_t levels;

    public:
#define FunctorArgs(...)                                                                                             \
    typename std::enable_if<std::is_convertible<Functor, std::function<void(__VA_ARGS__)>>::value, Functor>::type* = \
        nullptr


        // Functor: void f(index,i,j,k)
        template <typename Functor, FunctorArgs(idx_t, idx_t, idx_t, idx_t)>
        void operator()(const Functor& f) const {
            ATLAS_ASSERT(levels);
            if (global) {
                auto mpi_rank = mpi::comm(fs_.mpi_comm()).rank();
                if (owner == mpi_rank) {
                    const idx_t ny = fs_.grid().ny();
                    std::vector<idx_t> offset(ny);
                    offset[0] = 0;
                    for (idx_t j = 1; j < ny; ++j) {
                        offset[j] = offset[j - 1] + fs_.grid().nx(j - 1);
                    }
                    atlas_omp_parallel_for(idx_t j = 0; j < ny; ++j) {
                        idx_t index = offset[j];
                        for (auto i = 0; i < fs_.grid().nx(j); ++i, ++index) {
                            for (auto k = 0; k < levels; ++k) {
                                f(index, i, j, k);
                            }
                        }
                    }
                }
            }
            else {
                for (auto j = fs_.j_begin(); j < fs_.j_end(); ++j) {
                    for (auto i = fs_.i_begin(j); i < fs_.i_end(j); ++i) {
                        for (auto k = 0; k < levels; ++k) {
                            f(fs_.index(i, j), i, j, k);
                        }
                    }
                }
            }
        }

        // Functor: void f(index,i,j)
        template <typename Functor, FunctorArgs(idx_t, idx_t, idx_t)>
        void operator()(const Functor& f) const {
            ATLAS_ASSERT(levels == 0);
            if (global) {
                auto mpi_rank = mpi::comm(fs_.mpi_comm()).rank();
                if (owner == mpi_rank) {
                    const idx_t ny = fs_.grid().ny();
                    std::vector<idx_t> offset(ny);
                    offset[0] = 0;
                    for (idx_t j = 1; j < ny; ++j) {
                        offset[j] = offset[j - 1] + fs_.grid().nx(j - 1);
                    }
                    atlas_omp_parallel_for(idx_t j = 0; j < ny; ++j) {
                        idx_t index = offset[j];
                        for (idx_t i = 0; i < fs_.grid().nx(j); ++i, ++index) {
                            f(index, i, j);
                        }
                    }
                }
            }
            else {
                for (auto j = fs_.j_begin(); j < fs_.j_end(); ++j) {
                    for (auto i = fs_.i_begin(j); i < fs_.i_end(j); ++i) {
                        f(fs_.index(i, j), i, j);
                    }
                }
            }
        }

        // Functor: void f(index,k)
        template <typename Functor, FunctorArgs(idx_t, idx_t)>
        void operator()(const Functor& f) const {
            ATLAS_ASSERT(levels);
            if (global) {
                auto mpi_rank = mpi::comm(fs_.mpi_comm()).rank();
                if (owner == mpi_rank) {
                    const idx_t size = fs_.grid().size();
                    atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) {
                        for (idx_t k = 0; k < levels; ++k) {
                            f(n, k);
                        }
                    }
                }
            }
            else {
                const idx_t size = fs_.sizeOwned();
                atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) {
                    for (idx_t k = 0; k < levels; ++k) {
                        f(n, k);
                    }
                }
            }
        }


        // Functor: void f(index)
        template <typename Functor, FunctorArgs(idx_t)>
        void operator()(const Functor& f) const {
            ATLAS_ASSERT(levels == 0);
            if (global) {
                auto mpi_rank = mpi::comm(fs_.mpi_comm()).rank();
                if (owner == mpi_rank) {
                    const idx_t size = fs_.grid().size();
                    atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) { f(n); }
                }
            }
            else {
                const idx_t size = fs_.sizeOwned();
                atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) { f(n); }
            }
        }

#undef FunctorArgs
    };

    template <typename Functor>
    void parallel_for(const Functor& f) const {
        For(*this, util::NoConfig())(f);
    }
    template <typename Functor>
    void parallel_for(const util::Config& config, const Functor& f) const {
        For(*this, config)(f);
    }

private:
    const detail::StructuredColumns* functionspace_;
    void setup(const Grid& grid, const Vertical& vertical, const grid::Distribution& distribution,
               const eckit::Configuration& config);
};

// -------------------------------------------------------------------


}  // namespace functionspace
}  // namespace atlas
