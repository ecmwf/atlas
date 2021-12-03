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

#include <string>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"

//----------------------------------------------------------------------------------------------------------------------

// forward declarations
namespace atlas {
class Grid;
class FunctionSpace;
class Domain;
namespace trans {
class TransImpl;
class Cache;
}  // namespace trans
}  // namespace atlas

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace trans {

//----------------------------------------------------------------------------------------------------------------------

class TransFactory : public util::Factory<TransFactory> {
public:
    static std::string className() { return "TransFactory"; }

public:
    /*!
   * \brief build Trans
   * \return TransImpl
   */
    static const TransImpl* build(const FunctionSpace& gp, const FunctionSpace& sp,
                                  const eckit::Configuration& = util::Config());
    static const TransImpl* build(const Grid&, int truncation, const eckit::Configuration& = util::Config());

    static const TransImpl* build(const Grid&, const Domain&, int truncation,
                                  const eckit::Configuration& = util::Config());

    static const TransImpl* build(const Cache&, const FunctionSpace& gp, const FunctionSpace& sp,
                                  const eckit::Configuration& = util::Config());

    static const TransImpl* build(const Cache&, const Grid&, int truncation,
                                  const eckit::Configuration& = util::Config());

    static const TransImpl* build(const Cache&, const Grid&, const Domain&, int truncation,
                                  const eckit::Configuration& = util::Config());

    static void list(std::ostream& out);

    static bool has(const std::string& backend);

    static void backend(const std::string&);

    static std::string backend();

    static const eckit::Configuration& config();

    static void config(const eckit::Configuration&);

public:
    TransFactory(const std::string& name, const std::string& backend);
    virtual ~TransFactory();

    TransFactory() {}

private:
    std::string name_;
    std::string backend_;

    virtual const TransImpl* make(const Cache&, const FunctionSpace& /*gp*/, const FunctionSpace& /*sp*/,
                                  const eckit::Configuration&) {
        return nullptr;
    }
    virtual const TransImpl* make(const Cache&, const Grid& /*gp*/, const Domain&, int /*truncation*/,
                                  const eckit::Configuration&) {
        return nullptr;
    }
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class TransBuilderFunctionSpace : public TransFactory {
    virtual const TransImpl* make(const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
                                  const eckit::Configuration& config) override {
        return new T(cache, gp, sp, config);
    }
    virtual const TransImpl* make(const Cache&, const Grid&, const Domain&, int, const eckit::Configuration&) override {
        throw_Exception("This function should not be called", Here());
    }

public:
    TransBuilderFunctionSpace(const std::string& name, const std::string& backend): TransFactory(name, backend) {}

    TransBuilderFunctionSpace() {}
};

template <class T>
class TransBuilderGrid : public TransFactory {
    virtual const TransImpl* make(const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
                                  const eckit::Configuration& config) override {
        return new T(cache, grid, domain, truncation, config);
    }
    virtual const TransImpl* make(const Cache&, const FunctionSpace&, const FunctionSpace&,
                                  const eckit::Configuration&) override {
        throw_Exception("This function should not be called", Here());
    }

public:
    TransBuilderGrid(const std::string& name, const std::string& backend): TransFactory(name, backend) {}

    TransBuilderGrid() {}
};

//----------------------------------------------------------------------------------------------------------------------


}  // namespace trans
}  // namespace atlas
