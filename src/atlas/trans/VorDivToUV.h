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

#include <memory>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class VorDivToUVImpl : public util::Object {
public:
    virtual ~VorDivToUVImpl() = 0;

    virtual int truncation() const = 0;

    // -- IFS type fields --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    /*!
   * @brief Compute spectral wind (U/V) from spectral vorticity/divergence
   *
   * U = u*cos(lat)
   * V = v*cos(lat)
   *
   * @param nb_fields [in] Number of fields
   * @param vorticity [in] Spectral vorticity
   * @param divergence [in] Spectral divergence
   * @param U [out] Spectral wind U = u*cos(lat)
   * @param V [out] Spectral wind V = v*cos(lat)
   */
    virtual void execute(const int nb_coeff, const int nb_fields, const double vorticity[], const double divergence[],
                         double U[], double V[], const eckit::Configuration& = util::NoConfig()) const = 0;
};

// ------------------------------------------------------------------

class VorDivToUVFactory {
public:
    /*!
   * \brief build VorDivToUV
   * \return VorDivToUVImpl
   */
    static VorDivToUVImpl* build(const FunctionSpace& sp, const eckit::Configuration& = util::NoConfig());
    static VorDivToUVImpl* build(int truncation, const eckit::Configuration& = util::NoConfig());

    /*!
   * \brief list all registered trans implementations
   */
    static void list(std::ostream&);

    static bool has(const std::string& name);

private:
    std::string name_;
    virtual VorDivToUVImpl* make(const FunctionSpace& sp, const eckit::Configuration&) = 0;
    virtual VorDivToUVImpl* make(int truncation, const eckit::Configuration&)          = 0;

protected:
    VorDivToUVFactory(const std::string&);
    virtual ~VorDivToUVFactory();
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class VorDivToUVBuilder : public VorDivToUVFactory {
    virtual VorDivToUVImpl* make(const FunctionSpace& sp, const eckit::Configuration& config) {
        return new T(sp, config);
    }
    virtual VorDivToUVImpl* make(int truncation, const eckit::Configuration& config) {
        return new T(truncation, config);
    }

public:
    VorDivToUVBuilder(const std::string& name): VorDivToUVFactory(name) {}
};

//----------------------------------------------------------------------------------------------------------------------

class VorDivToUV : DOXYGEN_HIDE(public util::ObjectHandle<VorDivToUVImpl>) {
public:
    using Handle::Handle;
    VorDivToUV() = default;

    VorDivToUV(const FunctionSpace& sp, const eckit::Configuration& = util::NoConfig());
    VorDivToUV(int truncation, const eckit::Configuration& = util::NoConfig());

    void hash(eckit::Hash&) const;

    int truncation() const;

    // -- IFS type fields --
    // These fields have special interpretation required. You need to know what
    // you're doing.
    // See IFS trans library.

    virtual void execute(const int nb_coeff, const int nb_fields, const double vorticity[], const double divergence[],
                         double U[], double V[], const eckit::Configuration& = util::NoConfig()) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
