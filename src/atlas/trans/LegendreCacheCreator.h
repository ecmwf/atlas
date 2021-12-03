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

#include "atlas/trans/Trans.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Grid;
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class LegendreCacheCreatorImpl : public util::Object {
public:
    virtual ~LegendreCacheCreatorImpl() = 0;

    virtual bool supported() const = 0;

    virtual std::string uid() const = 0;

    virtual void create(const std::string& path) const = 0;

    virtual Cache create() const = 0;

    virtual size_t estimate() const = 0;
};

// ------------------------------------------------------------------

class LegendreCacheCreator : DOXYGEN_HIDE(public util::ObjectHandle<LegendreCacheCreatorImpl>) {
public:
    using Handle::Handle;
    LegendreCacheCreator() = default;
    LegendreCacheCreator(const Grid&, int truncation, const eckit::Configuration& = util::NoConfig());

    bool supported() const;

    std::string uid() const;

    void create(const std::string& path) const;

    Cache create() const;

    size_t estimate() const;
};

//----------------------------------------------------------------------------------------------------------------------

class LegendreCacheCreatorFactory {
public:
    /*!
   * \brief build Trans
   * \return TransImpl
   */
    static LegendreCacheCreatorImpl* build(const Grid&, int truncation, const eckit::Configuration& = util::Config());

    /*!
   * \brief list all registered trans implementations
   */
    static void list(std::ostream&);

    static bool has(const std::string& name);

private:
    std::string name_;
    virtual LegendreCacheCreatorImpl* make(const Grid& /*gp*/, int truncation, const eckit::Configuration&) {
        return nullptr;
    }

protected:
    LegendreCacheCreatorFactory(const std::string&);
    virtual ~LegendreCacheCreatorFactory();
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class LegendreCacheCreatorBuilder : public LegendreCacheCreatorFactory {
    virtual LegendreCacheCreatorImpl* make(const Grid& grid, int truncation, const eckit::Configuration& config) {
        return new T(grid, truncation, config);
    }

public:
    LegendreCacheCreatorBuilder(const std::string& name): LegendreCacheCreatorFactory(name) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
