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
#include "atlas/grid/RegularGG.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,RegularGG> RegularGG_builder;

//void RegularGG::constructFrom(const GridSpec& grid_spec)
//{
//   if (grid_spec.has("Nj"))           nj_ = grid_spec.get("Nj");
//   if (grid_spec.has("gaussianNumber")) gaussianNumber_ = grid_spec.get("gaussianNumber");
//   if (grid_spec.has("hash"))           hash_ = (std::string)grid_spec.get("hash");
//   grid_spec.get_bounding_box(bbox_);
//   grid_spec.get_latitudes(latitudes_);
//   grid_spec.get_points(points_);
//}

RegularGG::RegularGG( const eckit::Params& p )
{
	NOTIMP;
}

RegularGG::~RegularGG()
{
}

Grid::Point RegularGG::latLon(size_t the_i, size_t the_j) const
{
    /// @todo this function is VERY inneficient -- please rewrite it!

    long nptsWE = 4 * gaussianNumber_ ;
    long weIncrement = 360.0 / nptsWE;
    for(size_t i = 0 ; i < latitudes_.size(); i++) {

        double plon = bbox_.bottom_left().lon(); // west_;
        for( size_t j = 0; j < nptsWE; ++j) {
            if ( i== the_i && j== the_j) {
                return  Point( latitudes_[i], plon );
            }
            plon += weIncrement;
        }
    }

    return Grid::Point();
}

void RegularGG::coordinates(std::vector<double>& r ) const
{
	NOTIMP;
}

void RegularGG::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

string RegularGG::gridType() const
{
	return RegularGG::gridTypeStr();
}

GridSpec* RegularGG::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "GG" << gaussianNumber_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("Nj",eckit::Value(nj_));
   grid_spec->set("gaussianNumber",eckit::Value(gaussianNumber_));

   grid_spec->set("hash",eckit::Value(hash_));
   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_latitudes(latitudes_);
   grid_spec->set_points(points_);

   return grid_spec;
}

bool RegularGG::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const RegularGG&>(grid).gaussianNumber_ != gaussianNumber_) return false;
   if ( static_cast<const RegularGG&>(grid).nj_ != nj_) return false;
   if ( static_cast<const RegularGG&>(grid).hash_ != hash_) return false;
   if ( static_cast<const RegularGG&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const RegularGG&>(grid).latitudes_ != latitudes_) return false;
   if ( static_cast<const RegularGG&>(grid).points_ != points_) return false;

   return true;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
