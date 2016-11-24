/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#pragma once
#include "atlas/internals/atlas_defines.h"
#include "atlas/array/Array.h"

namespace atlas {
namespace array {
#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

  template <typename Value, unsigned int NDims, bool ReadOnly = false>
  inline data_view_tt<Value, NDims>
  make_gt_host_view(const Array& array);
#endif

}
}
