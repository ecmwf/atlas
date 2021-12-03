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

#include "atlas/util/Config.h"
#include "atlas/util/Registry.h"

namespace atlas {
namespace grid {

using Spec = util::Config;

class SpecRegistry {
public:
    static void add(const std::string& id, Spec&& spec) { registry().add(id, std::move(spec)); }

    static bool has(const std::string& id) { return registry().has(id); }

    static Spec get(const std::string& id) { return registry().get(id); }

    static void list(std::ostream& out) { return registry().list(out); }

    static std::vector<std::string> keys() { return registry().keys(); }

private:
    static util::Registry<Spec>& registry();
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
