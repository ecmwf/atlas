/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/grid/global/gaussian/ClassicGaussian.h"
#include "atlas/grid/global/gaussian/latitudes/Latitudes.h"

namespace atlas {
namespace grid {
namespace global {
namespace gaussian {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,ClassicGaussian,ClassicGaussian::grid_type_str());

std::string ClassicGaussian::className()
{
  return "atlas.grid.global.gaussian.ClassicGaussian";
}

void ClassicGaussian::set_typeinfo()
{
  std::stringstream s;
  s << "N" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

ClassicGaussian::ClassicGaussian()
{
}

ClassicGaussian::ClassicGaussian( const size_t N, const long nlons[], const Domain& d)
  : ReducedGaussianGrid(N,nlons,d)
{
  set_typeinfo();
}

ClassicGaussian::ClassicGaussian(const eckit::Parametrisation& params)
  : ReducedGaussianGrid(params)
{
  set_typeinfo();
}

//-----------------------------------------------------------------------------

} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas
