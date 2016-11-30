#ifndef atlas_RotatedLonLatProjection_H
#define atlas_RotatedLonLatProjection_H

#include "atlas/grid/projection/LonLatProjection.h"

namespace atlas {
namespace grid {
namespace projection {

class RotatedLonLatProjection: public LonLatProjection {
	public:
		// constructor
		RotatedLonLatProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.RotatedLonLatProjection"; }
		static std::string projection_type_str() {return "rotatedLonLat";}

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		// purely regional? - no!
		bool isRegional() { return false; }	// lonlat can be global
	
	private:
		eckit::geometry::LLPoint2 pole_;		// pole
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
