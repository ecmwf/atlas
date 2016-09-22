#ifndef atlas_LonLatProjection_H
#define atlas_LonLatProjection_H

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class LonLatProjection: public Projection {


	public:
		
		// constructor
		LonLatProjection(const eckit::Parametrisation& p);
		
		// class name
		static std::string className() { return "atlas.LonLatProjection"; }
		static std::string projection_type_str() {return "lonlat";}

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);
		
		// purely regional? - no!
		bool isRegional() { return false; }
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
