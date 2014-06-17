/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <numeric>

#include "atlas/meshgen/RGG.hpp"

#define DEBUG_OUTPUT 0

namespace atlas {
namespace meshgen {

int RGG::ngptot() const
{
  return std::accumulate(lon_.data(),lon_.data()+lon_.size(),0);
}  

} // namespace meshgen
} // namespace atlas

