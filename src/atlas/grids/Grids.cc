/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <eckit/log/Log.h>
#include "atlas/grids/Grids.h"

extern "C"
{
  void atlas__grids__load()
  {
    eckit::Log::info() << "Loading library [atlas::grids]" << std::endl;
  }
}
