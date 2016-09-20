#ifndef atlas_LambertProjection_H
#define atlas_LambertProjection_H

#include "atlas/grid/projections/Projection.h"

namespace atlas {
namespace grid {

class LambertProjection: public Projection {
	public:
	
		// constructor
		LambertProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.LambertProjection"; }

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		bool isRegional() { return true; }	// lambert projection cannot be used for global grids

	private:

		double lat1_, lat2_;		// secant latitudes
		bool isTangent_;
		double lon0_;						// central longitude
		double radius_;					// sphere radius
		double n_, F_, rho0_;		// projection constants
		
		void setup();
};

}  // namespace grid
}  // namespace atlas


#endif
