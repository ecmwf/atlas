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

#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/library/config.h"
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

class Spectral : public FunctionSpaceImpl {
public:
    Spectral( const eckit::Configuration& );

    Spectral( const int truncation, const eckit::Configuration& = util::NoConfig() );

    Spectral( const trans::Trans&, const eckit::Configuration& = util::NoConfig() );

    virtual ~Spectral() override;

    virtual std::string type() const override { return "Spectral"; }

    virtual std::string distribution() const override;

    /// @brief Create a spectral field
    using FunctionSpaceImpl::createField;
    virtual Field createField( const eckit::Configuration& ) const override;
    virtual Field createField( const Field&, const eckit::Configuration& ) const override;

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;

    void norm( const Field&, double& norm, int rank = 0 ) const;
    void norm( const Field&, double norm_per_level[], int rank = 0 ) const;
    void norm( const Field&, std::vector<double>& norm_per_level, int rank = 0 ) const;

public:  // methods
    idx_t nb_spectral_coefficients() const;
    idx_t nb_spectral_coefficients_global() const;
    int truncation() const { return truncation_; }

    virtual idx_t size() const override { return nb_spectral_coefficients(); }

private:  // methods
    array::DataType config_datatype( const eckit::Configuration& ) const;
    std::string config_name( const eckit::Configuration& ) const;
    idx_t config_size( const eckit::Configuration& ) const;
    idx_t config_levels( const eckit::Configuration& ) const;
    void set_field_metadata( const eckit::Configuration&, Field& ) const;
    size_t footprint() const override;

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
    Spectral( const FunctionSpace& );
    Spectral( const eckit::Configuration& );
    Spectral( const int truncation, const eckit::Configuration& = util::NoConfig() );
    Spectral( const trans::Trans&, const eckit::Configuration& = util::NoConfig() );

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;

    void norm( const Field&, double& norm, int rank = 0 ) const;
    void norm( const Field&, double norm_per_level[], int rank = 0 ) const;
    void norm( const Field&, std::vector<double>& norm_per_level, int rank = 0 ) const;

    idx_t nb_spectral_coefficients() const;
    idx_t nb_spectral_coefficients_global() const;
    int truncation() const;

private:
    const detail::Spectral* functionspace_;
};

}  // namespace functionspace
}  // namespace atlas

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
namespace atlas {
namespace field {
class FieldSetImpl;
class FieldImpl;
}  // namespace field
namespace trans {
class TransImpl;
}
namespace functionspace {

extern "C" {
const detail::Spectral* atlas__SpectralFunctionSpace__new__config( const eckit::Configuration* config );
const detail::Spectral* atlas__SpectralFunctionSpace__new__trans( trans::TransImpl* trans,
                                                                  const eckit::Configuration* config );
void atlas__SpectralFunctionSpace__delete( detail::Spectral* This );
field::FieldImpl* atlas__fs__Spectral__create_field( const detail::Spectral* This,
                                                     const eckit::Configuration* options );
void atlas__SpectralFunctionSpace__gather( const detail::Spectral* This, const field::FieldImpl* local,
                                           field::FieldImpl* global );
void atlas__SpectralFunctionSpace__gather_fieldset( const detail::Spectral* This, const field::FieldSetImpl* local,
                                                    field::FieldSetImpl* global );
void atlas__SpectralFunctionSpace__scatter( const detail::Spectral* This, const field::FieldImpl* global,
                                            field::FieldImpl* local );
void atlas__SpectralFunctionSpace__scatter_fieldset( const detail::Spectral* This, const field::FieldSetImpl* global,
                                                     field::FieldSetImpl* local );
void atlas__SpectralFunctionSpace__norm( const detail::Spectral* This, const field::FieldImpl* field, double norm[],
                                         int rank );
}

}  // namespace functionspace
}  // namespace atlas
