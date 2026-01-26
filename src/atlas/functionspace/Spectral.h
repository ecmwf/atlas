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
  
  const auto zonal_wavenumbers = Spectral::zonal_wavenumbers();
  const int truncation = Spectral::truncation();
  idx_t jc=0;
  for( int jm=0; jm<zonal_wavenumbers.size(); ++jm ) {
    int m = zonal_wavenumbers(jm);
    for( int n=m; m<=truncation; ++n ) {
      data( jc++, jfld ) = func_real_part(m,n);
      data( jc++, jfld ) = func_imag_part(m,n);
    }
  }
  
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

    std::string checksum(const FieldSet&) const;
    std::string checksum(const Field&) const;

    void norm(const Field&, double& norm, int rank = 0) const;
    void norm(const Field&, double norm_per_level[], int rank = 0) const;
    void norm(const Field&, std::vector<double>& norm_per_level, int rank = 0) const;

    array::LocalView<const int, 1> zonal_wavenumbers() const;  // zero-based, OK

    idx_t levels() const { return nb_levels_; }

    class For {
    public:
        For(const Spectral& fs, const util::Config& config = util::NoConfig()):
            truncation{fs.truncation()},
            zonal_wavenumbers{fs.zonal_wavenumbers()},
            global{config.getBool("global", false)},
            owner{config.getInt("owner", 0)} {}

    protected:
        using View = const array::LocalView<const int, 1>;
        int truncation;
        View zonal_wavenumbers;
        bool global;
        idx_t owner;

    public:
#define FunctorArgs(...)                                                                                             \
    typename std::enable_if<std::is_convertible<Functor, std::function<void(__VA_ARGS__)>>::value, Functor>::type* = \
        nullptr

        // Functor: void f(real,imag,n,m)
        template <typename Functor, FunctorArgs(idx_t, idx_t, int, int)>
        void operator()(const Functor& f) const {
            idx_t index = 0;
            if (global) {
                if (owner == mpi::rank()) {
                    for(int m = 0; m <= truncation; ++m) {
                        for (int n = m; n <= truncation; ++n) {
                            f(index, index + 1, n, m);
                            index += 2;
                        }
                    }
                }
            }
            else {
                const int nb_zonal_wavenumbers{static_cast<int>(zonal_wavenumbers.size())};
                for(int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
                    const int m = zonal_wavenumbers(jm);
                    for (int n = m; n <= truncation; ++n) {
                        f(index, index + 1, n, m);
                        index += 2;
                    }
                }
            }
        }

        // Functor: void f(real,imag,n)
        template <typename Functor, FunctorArgs(idx_t, idx_t, int)>
        void operator()(const Functor& f) const {
            idx_t index = 0;
            if (global) {
                if (owner == mpi::rank()) {
                    for(int m = 0; m <= truncation; ++m) {
                        for (int n = m; n <= truncation; ++n) {
                            f(index, index + 1, n);
                            index += 2;
                        }
                    }
                }
            }
            else {
                const int nb_zonal_wavenumbers{static_cast<int>(zonal_wavenumbers.size())};
                for(int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
                    const int m = zonal_wavenumbers(jm);
                    for (int n = m; n <= truncation; ++n) {
                        f(index, index + 1, n);
                        index += 2;
                    }
                }
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
    int nump() const;                               // equivalent to nmyms().size()
    array::LocalView<const int, 1> nvalue() const;  // Return wave number n for a given index
    array::LocalView<const int, 1> nmyms() const;   // Return list of local zonal wavenumbers "m"
    array::LocalView<const int, 1> nasm0() const;

private:  // data
    idx_t nb_levels_;
    int truncation_;

    class Parallelisation;
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

    array::LocalView<const int, 1> zonal_wavenumbers() const;  // zero-based, OK

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
