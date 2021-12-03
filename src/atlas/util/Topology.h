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

#include "atlas/util/Bitflags.h"

namespace atlas {
namespace util {

class Topology : public util::Bitflags {
public:
    enum
    {
        NONE     = 0,
        GHOST    = (1 << 1),
        PERIODIC = (1 << 2),
        BC       = (1 << 3),
        WEST     = (1 << 4),
        EAST     = (1 << 5),
        NORTH    = (1 << 6),
        SOUTH    = (1 << 7),
        PATCH    = (1 << 8),
        POLE     = (1 << 9),
        LAND     = (1 << 10),
        WATER    = (1 << 11),
        INVALID  = (1 << 12),
    };
};

}  // namespace util
}  // namespace atlas
