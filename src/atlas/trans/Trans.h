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

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/trans/Cache.h"
#include "atlas/util/Config.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
class Domain;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class TransImpl : public eckit::Owned {
public:
    virtual ~TransImpl() = 0;

    virtual int truncation() const = 0;

    virtual size_t spectralCoefficients() const = 0;

    virtual const Grid& grid() const = 0;

    virtual void dirtrans( const Field& gpfield, Field& spfield,
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void dirtrans( const FieldSet& gpfields, FieldSet& spfields,
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                       const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void invtrans( const Field& spfield, Field& gpfield,
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void invtrans( const FieldSet& spfields, FieldSet& gpfields,
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void invtrans_grad( const Field& spfield, Field& gradfield,
                                const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void invtrans_grad( const FieldSet& spfields, FieldSet& gradfields,
                                const eckit::Configuration& = util::NoConfig() ) const = 0;

    virtual void invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                       const eckit::Configuration& = util::NoConfig() ) const = 0;

    // -- IFS type fields --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
   * @brief invtrans
   * @param nb_scalar_fields
   * @param scalar_spectra
   * @param nb_vordiv_fields
   * @param vorticity_spectra
   * @param divergence_spectra
   * @param gp_fields
   */
    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                           const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    /*!
   * @brief invtrans
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    virtual void invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    /*!
   * @brief Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    virtual void invtrans( const int nb_vordiv_fields, const double vorticity_spectra[],
                           const double divergence_spectra[], double gp_fields[],
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    /*!
   * @brief Direct transform of scalar fields
   */
    virtual void dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                           const eckit::Configuration& = util::NoConfig() ) const = 0;

    /*!
   * @brief Direct transform of wind(U/V) to vorticity/divergence
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    virtual void dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                           double divergence_spectra[], const eckit::Configuration& = util::NoConfig() ) const = 0;
};

// ------------------------------------------------------------------

class TransFactory {
protected:
    using Trans_t = const TransImpl;

public:
    /*!
   * \brief build Trans
   * \return TransImpl
   */
    static Trans_t* build( const FunctionSpace& gp, const FunctionSpace& sp,
                           const eckit::Configuration& = util::Config() );
    static Trans_t* build( const Grid&, int truncation, const eckit::Configuration& = util::Config() );

    static Trans_t* build( const Grid&, const Domain&, int truncation, const eckit::Configuration& = util::Config() );

    static Trans_t* build( const Cache&, const FunctionSpace& gp, const FunctionSpace& sp,
                           const eckit::Configuration& = util::Config() );

    static Trans_t* build( const Cache&, const Grid&, int truncation, const eckit::Configuration& = util::Config() );

    static Trans_t* build( const Cache&, const Grid&, const Domain&, int truncation,
                           const eckit::Configuration& = util::Config() );

    /*!
   * \brief list all registered trans implementations
   */
    static void list( std::ostream& );

    static bool has( const std::string& name );

    static void backend( const std::string& );

    static std::string backend();

    static void config( const eckit::Configuration& );

    static const eckit::Configuration& config();

private:
    std::string name_;
    static util::Config default_options_;
    virtual Trans_t* make( const Cache&, const FunctionSpace& gp, const FunctionSpace& sp,
                           const eckit::Configuration& ) {
        return nullptr;
    }
    virtual Trans_t* make( const Cache&, const Grid& gp, const Domain&, int truncation, const eckit::Configuration& ) {
        return nullptr;
    }

protected:
    TransFactory( const std::string& );
    virtual ~TransFactory();
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class TransBuilderFunctionSpace : public TransFactory {
    virtual Trans_t* make( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
                           const eckit::Configuration& config ) {
        return new T( cache, gp, sp, config );
    }
    virtual Trans_t* make( const Cache&, const Grid&, const Domain&, int, const eckit::Configuration& ) {
        throw eckit::SeriousBug( "This function should not be called", Here() );
    }

public:
    TransBuilderFunctionSpace( const std::string& name ) : TransFactory( name ) {}
};

template <class T>
class TransBuilderGrid : public TransFactory {
    virtual Trans_t* make( const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
                           const eckit::Configuration& config ) {
        return new T( cache, grid, domain, truncation, config );
    }
    virtual Trans_t* make( const Cache&, const FunctionSpace&, const FunctionSpace&, const eckit::Configuration& ) {
        throw eckit::SeriousBug( "This function should not be called", Here() );
    }

public:
    TransBuilderGrid( const std::string& name ) : TransFactory( name ) {}
};

//----------------------------------------------------------------------------------------------------------------------

class Trans {
public:
    using Implementation = const TransImpl;

private:
    eckit::SharedPtr<Implementation> impl_;

public:
    static bool hasBackend( const std::string& );
    static void backend( const std::string& );
    static std::string backend();
    static void config( const eckit::Configuration& );
    static const eckit::Configuration& config();

    Trans();
    Trans( Implementation* );
    Trans( const Trans& );

    Trans( const FunctionSpace& gp, const FunctionSpace& sp, const eckit::Configuration& = util::NoConfig() );
    Trans( const Grid&, int truncation, const eckit::Configuration& = util::NoConfig() );
    Trans( const Grid&, const Domain&, int truncation, const eckit::Configuration& = util::NoConfig() );

    Trans( const Cache&, const FunctionSpace& gp, const FunctionSpace& sp,
           const eckit::Configuration& = util::NoConfig() );
    Trans( const Cache&, const Grid&, int truncation, const eckit::Configuration& = util::NoConfig() );
    Trans( const Cache&, const Grid&, const Domain&, int truncation, const eckit::Configuration& = util::NoConfig() );

    void hash( eckit::Hash& ) const;
    Implementation* get() const { return impl_.get(); }
    operator bool() const { return impl_.owners(); }

    int truncation() const;
    size_t spectralCoefficients() const;
    const Grid& grid() const;

    void dirtrans( const Field& gpfield, Field& spfield, const eckit::Configuration& = util::NoConfig() ) const;

    void dirtrans( const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& = util::NoConfig() ) const;

    void dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                               const eckit::Configuration& = util::NoConfig() ) const;

    void invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& = util::NoConfig() ) const;

    void invtrans( const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& = util::NoConfig() ) const;

    void invtrans_grad( const Field& spfield, Field& gradfield, const eckit::Configuration& = util::NoConfig() ) const;

    void invtrans_grad( const FieldSet& spfields, FieldSet& gradfields,
                        const eckit::Configuration& = util::NoConfig() ) const;

    void invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                               const eckit::Configuration& = util::NoConfig() ) const;

    // -- IFS type fields --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
   * @brief invtrans
   * @param nb_scalar_fields
   * @param scalar_spectra
   * @param nb_vordiv_fields
   * @param vorticity_spectra
   * @param divergence_spectra
   * @param gp_fields
   */
    void invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                   const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                   const eckit::Configuration& = util::NoConfig() ) const;

    /*!
   * @brief invtrans
   * @param nb_fields
   * @param scalar_spectra
   * @param scalar_fields
   */
    void invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                   const eckit::Configuration& = util::NoConfig() ) const;

    /*!
   * @brief Inverse transform of vorticity/divergence to wind(U/V)
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    void invtrans( const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                   double gp_fields[], const eckit::Configuration& = util::NoConfig() ) const;

    /*!
   * @brief Direct transform of scalar fields
   */
    void dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                   const eckit::Configuration& = util::NoConfig() ) const;

    /*!
   * @brief Direct transform of wind(U/V) to vorticity/divergence
   * @param nb_fields [in] Number of fields ( both components of wind count as 1
   * )
   */
    void dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                   double divergence_spectra[], const eckit::Configuration& = util::NoConfig() ) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
