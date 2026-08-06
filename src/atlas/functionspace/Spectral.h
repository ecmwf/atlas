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

#include "atlas/array/LocalView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/library/config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
class Field;
class FieldSet;
}  // namespace atlas

namespace atlas {
namespace trans {
class Trans;
}
}  // namespace atlas

namespace atlas {
namespace functionspace {
namespace detail {

// -------------------------------------------------------------------

class Spectral : public functionspace::FunctionSpaceImpl {
    /*
  Spectral data is organised as:
     m = zonal wavenumber
     n = total wavenumber

  const auto zonal_wavenumbers = spectral.zonal_wavenumbers();
  const int truncation = spectral.truncation();
  idx_t index = 0;
  for( idx_t jm=0; jm<zonal_wavenumbers.size(); ++jm ) {
      const int m = zonal_wavenumbers(jm);
      for( int n=m; n<=truncation; ++n ) {
          data( index,   level ) = func_real_part(n,m,level);
          data( index+1, level ) = func_imag_part(n,m,level);
          index += 2;
      }
  }

  Alternatively, offsets_by_zonal_wavenumber() can be used as the base offset
  for each zonal wavenumber.

  const auto offsets_by_zonal_wavenumber = spectral.offsets_by_zonal_wavenumber();
  for( idx_t jm=0; jm<zonal_wavenumbers.size(); ++jm ) {
      const int m = zonal_wavenumbers(jm);
      const idx_t offset_for_zonal_wavenumber = offsets_by_zonal_wavenumber[m];
      for( int n=m; n<=truncation; ++n ) {
          const idx_t index = offset_for_zonal_wavenumber + 2 * (n - m);
          data( index,   level ) = func_real_part(n,m,level);
          data( index+1, level ) = func_imag_part(n,m,level);
      }
  }

  Or inverting the loop order, less efficient.

  for( int n=0; n<=truncation; ++n ) {
      for( idx_t jm=0; jm<zonal_wavenumbers.size(); ++jm ) {
          const int m = zonal_wavenumbers(jm);
          if( m > n ) {
              continue;
          }
          const idx_t offset_for_zonal_wavenumber = offsets_by_zonal_wavenumber[m];
          const idx_t index = offset_for_zonal_wavenumber + 2 * (n - m);
          data( index,   level ) = func_real_part(n,m,level);
          data( index+1, level ) = func_imag_part(n,m,level);
      }
  }

  The same storage order can also be accessed with parallel_for().

  spectral.parallel_for([&]( idx_t real, idx_t imag, int n, int m ) {
      data( real, level ) = func_real_part(n,m,level);
      data( imag, level ) = func_imag_part(n,m,level);
  });
*/

public:
    Spectral(const eckit::Configuration&);

    Spectral(const int truncation, const eckit::Configuration& = util::NoConfig());

    ~Spectral() override;

    std::string type() const override { return "Spectral"; }

    std::string distribution() const override;

    idx_t part() const override;

    idx_t nb_parts() const override;

    using FunctionSpaceImpl::createField;
    Field createField(const eckit::Configuration&) const override;
    Field createField(const Field&, const eckit::Configuration&) const override;

    using FunctionSpaceImpl::gather;
    void gather(const FieldSet&, FieldSet&) const override;
    void gather(const Field&, Field&) const override;

    using FunctionSpaceImpl::scatter;
    void scatter(const FieldSet&, FieldSet&) const override;
    void scatter(const Field&, Field&) const override;

    std::string checksum(const FieldSet&) const override;
    std::string checksum(const Field&) const override;

    void norm(const Field&, double& norm, int rank = 0) const;
    void norm(const Field&, double norm_per_level[], int rank = 0) const;
    void norm(const Field&, std::vector<double>& norm_per_level, int rank = 0) const;

    array::LocalView<const int, 1> zonal_wavenumbers() const;  // zero-based
    array::LocalView<const int, 1> offsets_by_zonal_wavenumber() const;  // zero-based nasm0 for use in C++

    idx_t levels() const { return nb_levels_; }

    class For {
    public:
        For(const Spectral& fs, const util::Config& config = util::NoConfig()):
            truncation{fs.truncation()},
            zonal_wavenumbers{fs.zonal_wavenumbers()},
            offsets_by_zonal_wavenumber{fs.offsets_by_zonal_wavenumber()},
            global{config.getBool("global", false)},
            owner{config.getInt("owner", 0)} {}

    protected:
        using View = const array::LocalView<const int, 1>;
        int truncation;
        View zonal_wavenumbers;
        View offsets_by_zonal_wavenumber;
        bool global;
        idx_t owner;

    public:
#define FunctorArgs(...)                                                                                             \
    typename std::enable_if<std::is_convertible<Functor, std::function<void(__VA_ARGS__)>>::value, Functor>::type* = \
        nullptr

        // Functor: void f(real,imag,n,m)
        template <typename Functor, FunctorArgs(idx_t, idx_t, int, int)>
        void operator()(const Functor& f) const {
            if (global) {
                if (owner == mpi::rank()) {
                    atlas_omp_parallel_for(int m = 0; m <= truncation; ++m) {
                        idx_t index = global_offset_by_zonal_wavenumber(m);
                        for (int n = m; n <= truncation; ++n, index += 2) {
                            f(index, index + 1, n, m);
                        }
                    }

                }
            }
            else {
                const int nb_zonal_wavenumbers{static_cast<int>(zonal_wavenumbers.size())};
                atlas_omp_parallel_for(int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
                    const int m = zonal_wavenumbers(jm);
                    idx_t index = offsets_by_zonal_wavenumber[m];
                    for (int n = m; n <= truncation; ++n, index += 2) {
                        f(index, index + 1, n, m);
                    }
                }
            }
        }

        // Functor: void f(real,imag,n)
        template <typename Functor, FunctorArgs(idx_t, idx_t, int)>
        void operator()(const Functor& f) const {
            if (global) {
                if (owner == mpi::rank()) {
                    atlas_omp_parallel_for(int m = 0; m <= truncation; ++m) {
                        idx_t index = global_offset_by_zonal_wavenumber(m);
                        for (int n = m; n <= truncation; ++n, index += 2) {
                            f(index, index + 1, n);
                        }
                    }
                }
            }
            else {
                const int nb_zonal_wavenumbers{static_cast<int>(zonal_wavenumbers.size())};
                atlas_omp_parallel_for(int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
                    const int m = zonal_wavenumbers(jm);
                    idx_t index = offsets_by_zonal_wavenumber[m];
                    for (int n = m; n <= truncation; ++n, index += 2) {
                        f(index, index + 1, n);
                    }
                }
            }
        }
        idx_t global_offset_by_zonal_wavenumber(const int m) const {
            // For the global packed triangular layout:
            // idx_t index = 0;
            // for(int m = 0; m <= truncation; ++m) {
            //     for (int n = m; n <= truncation; ++n) {
            //         f(index, index + 1, n);
            //         index += 2;
            //     }
            // }

            // The value `index` can be computed directly rather than accumulated:
            //    index = offset(m) + 2 * (n - m)
            // where
            //    offset(m) = sum_{k=0}^{m-1} 2 * (truncation - k + 1)
            // which simplifies to:
            //    offset(m) = m * (2 * truncation + 3 - m)
            return static_cast<idx_t>(m) * (static_cast<idx_t>(truncation) * 2 + 3 - m);
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

public:  // methods
    idx_t nb_spectral_coefficients() const;
    idx_t nb_spectral_coefficients_global() const;
    int truncation() const { return truncation_; }

    idx_t size() const override { return nb_spectral_coefficients(); }

private:  // methods
    array::DataType config_datatype(const eckit::Configuration&) const;
    std::string config_name(const eckit::Configuration&) const;
    idx_t config_size(const eckit::Configuration&) const;
    idx_t config_levels(const eckit::Configuration&) const;
    void set_field_metadata(const eckit::Configuration&, Field&) const;
    size_t footprint() const override;


private:  // Fortran access
    friend struct SpectralFortranAccess;
    int nump() const;                               // Number of zonal wave numbers m on THIS rank, equivalent to nmyms().size()
    array::LocalView<const int, 1> nvalue() const;  // Return wave number n for a given index
    array::LocalView<const int, 1> nmyms() const;   // Array of actual m values (zonal wave numbers) on this rank (size nump)
    array::LocalView<const int, 1> nasm0_base0() const;   // Base offset in memory for this zonal wave number
    array::LocalView<const int, 1> nasm0_base1() const;   // Base offset in memory for this zonal wave number

private:  // data
    idx_t nb_levels_;
    int truncation_;

    class Parallelisation;
    friend class Parallelisation_ectrans;
    friend class Parallelisation_local;
    std::unique_ptr<Parallelisation> parallelisation_;
};

}  // namespace detail

// -------------------------------------------------------------------

class Spectral : public FunctionSpace {
public:
    Spectral();
    Spectral(const FunctionSpace&);
    Spectral(const eckit::Configuration&);
    Spectral(const int truncation, const eckit::Configuration& = util::NoConfig());

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;

    void norm(const Field&, double& norm, int rank = 0) const;
    void norm(const Field&, double norm_per_level[], int rank = 0) const;
    void norm(const Field&, std::vector<double>& norm_per_level, int rank = 0) const;

    array::LocalView<const int, 1> zonal_wavenumbers() const;  // zero-based
    array::LocalView<const int, 1> offsets_by_zonal_wavenumber() const;  // zero-based nasm0 for use in C++

    idx_t nb_spectral_coefficients() const;
    idx_t nb_spectral_coefficients_global() const;
    int truncation() const;
    idx_t levels() const { return functionspace_->levels(); }

    template <typename Functor>
    void parallel_for(const Functor& f) const {
        functionspace_->parallel_for(f);
    }
    template <typename Functor>
    void parallel_for(const util::Config& config, const Functor& f) const {
        functionspace_->parallel_for(config, f);
    }

private:
    const detail::Spectral* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas
