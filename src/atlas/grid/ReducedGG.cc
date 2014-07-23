/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/value/Value.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/ReducedGG.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,ReducedGG> ReducedGG_builder("reduced_gg");

//void ReducedGG::constructFrom(const GridSpec& grid_spec)
//{
//   if (grid_spec.has("gaussianNumber")) gaussianNumber_ = grid_spec.get("gaussianNumber");
//   if (grid_spec.has("hash"))           hash_ = (std::string)grid_spec.get("hash");
//   grid_spec.get_bounding_box(bbox_);
//   grid_spec.get_rgspec(rgSpec_);
//   grid_spec.get_latitudes(latitudes_);
//   grid_spec.get_points(points_);
//}

ReducedGG::ReducedGG( const eckit::Params& p )
{
	NOTIMP;
}

ReducedGG::~ReducedGG()
{
}

void ReducedGG::coordinates( std::vector<double>& pts ) const
{
	NOTIMP;
}

void ReducedGG::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

string ReducedGG::gridType() const
{
	return ReducedGG::gridTypeStr();
}

GridSpec* ReducedGG::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "QG" << gaussianNumber_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("gaussianNumber",eckit::Value(gaussianNumber_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_rgspec(rgSpec_);
   grid_spec->set_latitudes(latitudes_);
   grid_spec->set_points(points_);

   return grid_spec;
}

bool ReducedGG::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const ReducedGG&>(grid).gaussianNumber_ != gaussianNumber_) return false;
   if ( static_cast<const ReducedGG&>(grid).hash_ != hash_) return false;
   if ( static_cast<const ReducedGG&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const ReducedGG&>(grid).rgSpec_ != rgSpec_) return false;
   if ( static_cast<const ReducedGG&>(grid).latitudes_ != latitudes_) return false;
   if ( static_cast<const ReducedGG&>(grid).points_ != points_) return false;

   return true;
}


//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
