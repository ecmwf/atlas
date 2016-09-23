#ifndef atlas_RotatedTiltedSchmidtProjection_H
#define atlas_RotatedTiltedSchmidtProjection_H

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class RotatedTiltedSchmidtProjection: public Projection {
	public:
		// constructor
		RotatedTiltedSchmidtProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.RotatedTiltedSchmidtProjection"; }
		static std::string projection_type_str() {return "rotatedTiltedSchmidt";}

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		// purely regional? - no!
		bool isRegional() { return false; }	// schmidt is global grid
	
	private:
		double c_;				// stretching factor
		eckit::geometry::LLPoint2 pole_;		// pole coordinates
		
		void rotate_tilt_( eckit::geometry::LLPoint2 & );	// auxiliary function for rotating/tilting
		void un_rotate_tilt_( eckit::geometry::LLPoint2 & );	// auxiliary function for un rotating/tilting
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
