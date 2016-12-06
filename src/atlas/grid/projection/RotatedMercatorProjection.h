#ifndef atlas_RotatedMercatorProjection_H
#define atlas_RotatedMercatorProjection_H

#include "atlas/grid/projection/MercatorProjection.h"

namespace atlas {
namespace grid {
namespace projection {

class RotatedMercatorProjection: public MercatorProjection {
	public:
	
		// constructor
		RotatedMercatorProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.RotatedMercatorProjection"; }
		static std::string projection_type_str() {return "mercator";}

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		bool isRegional() { return true; }	// lambert projection cannot be used for global grids

	private:

		eckit::geometry::LLPoint2 pole_;		// pole
		void setup(const eckit::Parametrisation & p);
		
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
