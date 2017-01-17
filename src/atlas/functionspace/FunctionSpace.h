/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_FunctionSpace_h
#define atlas_FunctionSpace_h

#include <string>
#include "eckit/memory/Owned.h"

namespace atlas {
namespace functionspace {

#define FunctionspaceT_nonconst typename FunctionSpace::remove_const<FunctionSpaceT>::type
#define FunctionspaceT_const    typename FunctionSpace::add_const<FunctionSpaceT>::type

/// @brief FunctionSpace class helps to interprete Fields.
/// @note  Abstract base class
class FunctionSpace : public eckit::Owned
{
private:
    template<typename T> struct remove_const          { typedef T type; };
    template<typename T> struct remove_const<T const> { typedef T type; };

    template<typename T> struct add_const          { typedef const typename remove_const<T>::type type; };
    template<typename T> struct add_const<T const> { typedef const T type; };

public:
    FunctionSpace() {}
    virtual ~FunctionSpace() = 0;
    virtual std::string name() const = 0;
    virtual operator bool() const { return true; }
    virtual size_t footprint() const = 0;

    template <typename FunctionSpaceT>
    FunctionspaceT_nonconst *cast();

    template <typename FunctionSpaceT>
    FunctionspaceT_const *cast() const;

};

inline FunctionSpace::~FunctionSpace() {}

template <typename FunctionSpaceT>
inline FunctionspaceT_nonconst *FunctionSpace::cast()
{
  return dynamic_cast<FunctionspaceT_nonconst*> (this);
}

template <typename FunctionSpaceT>
inline FunctionspaceT_const *FunctionSpace::cast() const
{
  return dynamic_cast<FunctionspaceT_const*> (this);
}

#undef FunctionspaceT_const
#undef FunctionspaceT_nonconst
//------------------------------------------------------------------------------------------------------

/// @brief Dummy Functionspace class that evaluates to false
class NoFunctionSpace : public FunctionSpace
{
public:
    NoFunctionSpace() {}
    virtual ~NoFunctionSpace() {}
    virtual std::string name() const { return "NoFunctionSpace"; }
    virtual operator bool() const { return false; }
    virtual size_t footprint() const { return sizeof(*this); }
};

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C"
{
    void atlas__FunctionSpace__delete (FunctionSpace* This);
    const char* atlas__FunctionSpace__name (FunctionSpace* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

#endif // atlas_FunctionSpace_h
