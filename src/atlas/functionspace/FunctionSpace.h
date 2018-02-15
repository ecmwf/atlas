/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include "atlas/field/Field.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {
namespace functionspace {

#define FunctionspaceT_nonconst typename FunctionSpaceImpl::remove_const<FunctionSpaceT>::type
#define FunctionspaceT_const typename FunctionSpaceImpl::add_const<FunctionSpaceT>::type

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpaceImpl : public eckit::Owned {
private:
    template <typename T>
    struct remove_const {
        typedef T type;
    };
    template <typename T>
    struct remove_const<T const> {
        typedef T type;
    };

    template <typename T>
    struct add_const {
        typedef const typename remove_const<T>::type type;
    };
    template <typename T>
    struct add_const<T const> {
        typedef const T type;
    };

public:
    FunctionSpaceImpl() {}
    virtual ~FunctionSpaceImpl()     = 0;
    virtual std::string type() const = 0;
    virtual operator bool() const { return true; }
    virtual size_t footprint() const = 0;

    virtual atlas::Field createField( const eckit::Configuration& ) const = 0;

    virtual atlas::Field createField( const atlas::Field&, const eckit::Configuration& ) const = 0;

    atlas::Field createField( const atlas::Field& ) const;

    template <typename DATATYPE>
    atlas::Field createField( const eckit::Configuration& ) const;

    template <typename DATATYPE>
    atlas::Field createField() const;

    const util::Metadata& metadata() const { return metadata_; }
    util::Metadata& metadata() { return metadata_; }

    template <typename FunctionSpaceT>
    FunctionspaceT_nonconst* cast();

    template <typename FunctionSpaceT>
    FunctionspaceT_const* cast() const;

    virtual std::string distribution() const = 0;

private:
    util::Metadata metadata_;
};

inline FunctionSpaceImpl::~FunctionSpaceImpl() {}

template <typename DATATYPE>
inline Field FunctionSpaceImpl::createField( const eckit::Configuration& options ) const {
    return createField( option::datatypeT<DATATYPE>() | options );
}

template <typename DATATYPE>
inline Field FunctionSpaceImpl::createField() const {
    return createField( option::datatypeT<DATATYPE>() );
}

template <typename FunctionSpaceT>
inline FunctionspaceT_nonconst* FunctionSpaceImpl::cast() {
    return dynamic_cast<FunctionspaceT_nonconst*>( this );
}

template <typename FunctionSpaceT>
inline FunctionspaceT_const* FunctionSpaceImpl::cast() const {
    return dynamic_cast<FunctionspaceT_const*>( this );
}

#undef FunctionspaceT_const
#undef FunctionspaceT_nonconst

//------------------------------------------------------------------------------------------------------

/// @brief Dummy Functionspace class that evaluates to false
class NoFunctionSpace : public FunctionSpaceImpl {
public:
    NoFunctionSpace() {}
    virtual ~NoFunctionSpace() {}
    virtual std::string type() const { return "NoFunctionSpace"; }
    virtual operator bool() const { return false; }
    virtual size_t footprint() const { return sizeof( *this ); }
    virtual std::string distribution() const { return std::string(); }

    virtual Field createField( const eckit::Configuration& ) const;
    virtual Field createField( const Field&, const eckit::Configuration& ) const;
};

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__FunctionSpace__delete( FunctionSpaceImpl* This );
void atlas__FunctionSpace__name( const FunctionSpaceImpl* This, char*& name, int& size );
field::FieldImpl* atlas__FunctionSpace__create_field( const FunctionSpaceImpl* This,
                                                      const eckit::Configuration* options );
field::FieldImpl* atlas__FunctionSpace__create_field_template( const FunctionSpaceImpl* This,
                                                               const field::FieldImpl* field_template,
                                                               const eckit::Configuration* options );
}

//------------------------------------------------------------------------------------------------------

}  // namespace functionspace

//------------------------------------------------------------------------------------------------------

class FunctionSpace {
public:
    using Implementation = functionspace::FunctionSpaceImpl;

private:
    eckit::SharedPtr<const Implementation> functionspace_;

public:
    FunctionSpace();
    FunctionSpace( const Implementation* );
    FunctionSpace( const FunctionSpace& );

    std::string type() const;
    operator bool() const;
    size_t footprint() const;
    std::string distribution() const;

    const Implementation* get() const { return functionspace_.get(); }

    atlas::Field createField( const eckit::Configuration& ) const;

    atlas::Field createField( const atlas::Field&, const eckit::Configuration& = util::NoConfig() ) const;

    template <typename DATATYPE>
    Field createField( const eckit::Configuration& = util::NoConfig() ) const;
};

template <typename DATATYPE>
Field FunctionSpace::createField( const eckit::Configuration& options ) const {
    return functionspace_->createField<DATATYPE>( options );
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
