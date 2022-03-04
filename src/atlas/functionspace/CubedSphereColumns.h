/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <functional>
#include <type_traits>

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/detail/CubedSphereStructure.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
class Mesh;
}

namespace atlas {
namespace functionspace {

/// Extend NodeColumns and CellColumns so that they can exploit CubedSphere structure.
template <typename BaseFunctionSpace>
class CubedSphereColumns : public BaseFunctionSpace {
public:
    /// Constructors.
    CubedSphereColumns();
    CubedSphereColumns(const FunctionSpace& functionSpace);
    CubedSphereColumns(const Mesh& mesh, const eckit::Configuration& configuration);
    CubedSphereColumns(const Mesh& mesh);

    /// Invalid index.
    idx_t invalid_index() const;

    /// Get number of owned elements.
    idx_t sizeOwned() const;

    /// i lower bound for tile t (including halo)
    idx_t i_begin(idx_t t) const;
    /// i lower bound for tile t (including halo)
    idx_t i_end(idx_t t) const;

    /// j lower bound for tile t (including halo)
    idx_t j_begin(idx_t t) const;
    /// j lower bound for tile t (including halo)
    idx_t j_end(idx_t t) const;

    /// Return array_view index for (t, i, j).
    idx_t index(idx_t t, idx_t i, idx_t j) const;

    /// Return true if (t, i, j) is a valid index.
    bool is_valid_index(idx_t t, idx_t i, idx_t j) const;

    /// Return tij field.
    Field tij() const;

private:
    class For {
    public:
        For(const CubedSphereColumns<BaseFunctionSpace>& functionSpace, const util::Config& config = util::NoConfig()):
            functionSpace_{functionSpace},
            indexMax_{config.getBool("include_halo", false) ? functionSpace.size() : functionSpace.sizeOwned()},
            levels_{config.getInt("levels", functionSpace_.levels())},
            tijView_(array::make_view<idx_t, 2>(functionSpace_.tij())) {}

        // Define template to disable invalid functors.
        template <typename FuncType, typename... ArgTypes>
        using EnableFunctor =
            typename std::enable_if<std::is_convertible<FuncType, std::function<void(ArgTypes...)>>::value>::type*;

        // Functor: void f(index, t, i, j, k)
        template <typename Functor, EnableFunctor<Functor, idx_t, idx_t, idx_t, idx_t, idx_t> = nullptr>
        void operator()(const Functor& f) const {
            using namespace meshgenerator::detail::cubedsphere;

            // Loop over elements.
            atlas_omp_parallel_for(idx_t index = 0; index < indexMax_; ++index) {
                const idx_t t = tijView_(index, Coordinates::T);
                const idx_t i = tijView_(index, Coordinates::I);
                const idx_t j = tijView_(index, Coordinates::J);
                for (idx_t k = 0; k < levels_; ++k) {
                    f(index, t, i, j, k);
                }
            }
        }

        // Functor: void f(index, t, i, j)
        template <typename Functor, EnableFunctor<Functor, idx_t, idx_t, idx_t, idx_t> = nullptr>
        void operator()(const Functor& f) const {
            using namespace meshgenerator::detail::cubedsphere;

            // Loop over elements.
            atlas_omp_parallel_for(idx_t index = 0; index < indexMax_; ++index) {
                const idx_t t = tijView_(index, Coordinates::T);
                const idx_t i = tijView_(index, Coordinates::I);
                const idx_t j = tijView_(index, Coordinates::J);
                f(index, t, i, j);
            }
        }

        // Functor: void f(index, k)
        template <typename Functor, EnableFunctor<Functor, idx_t, idx_t> = nullptr>
        void operator()(const Functor& f) const {
            using namespace meshgenerator::detail::cubedsphere;

            // Loop over elements.
            atlas_omp_parallel_for(idx_t index = 0; index < indexMax_; ++index) {
                for (idx_t k = 0; k < levels_; ++k) {
                    f(index, k);
                }
            }
        }

        // Functor: void f(index )
        template <typename Functor, EnableFunctor<Functor, idx_t> = nullptr>
        void operator()(const Functor& f) const {
            using namespace meshgenerator::detail::cubedsphere;

            // Loop over elements.
            atlas_omp_parallel_for(idx_t index = 0; index < indexMax_; ++index) { f(index); }
        }

    private:
        const CubedSphereColumns<BaseFunctionSpace>& functionSpace_;
        idx_t indexMax_;
        idx_t levels_;
        array::ArrayView<idx_t, 2> tijView_;
    };

public:
    /// Visit each element index and apply functor with signature:
    ///   f(index, t, i, j, k), f(index t, i, j), f(index, k) or f(index).
    template <typename Functor>
    void parallel_for(const Functor& f) const {
        For(*this, util::NoConfig())(f);
    }

    /// Visit each element index and apply functor with signature:
    ///   f(index, t, i, j, k), f(index t, i, j), f(index, k) or f(index).
    /// Can also specify "idx_t levels" and "bool include_halo".
    template <typename Functor>
    void parallel_for(const util::Config& config, const Functor& f) const {
        For(*this, config)(f);
    }

private:
    // Object hanldle for CubedSphereStructure
    util::ObjectHandle<detail::CubedSphereStructure> cubedSphereColumnsHandle_;
};

using CubedSphereNodeColumns = CubedSphereColumns<NodeColumns>;
using CubedSphereCellColumns = CubedSphereColumns<CellColumns>;

}  // namespace functionspace
}  // namespace atlas
