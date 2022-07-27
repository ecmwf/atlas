/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/defines.h"

namespace atlas {
namespace io {

enum class Endian
{
    little = 0,
    big    = 1,
#if ATLAS_BIG_ENDIAN
    native  = big,
    swapped = little
#elif ATLAS_LITTLE_ENDIAN
    native  = little,
    swapped = big
#else
#error Neither ATLAS_BIG_ENDIAN nor ATLAS_LITTLE_ENDIAN equals true
#endif
};


}  // namespace io
}  // namespace atlas
