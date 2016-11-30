#ifndef atlas_SchmidtProjection_H
#define atlas_SchmidtProjection_H

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class SchmidtProjection: public Projection {
	public:
		// constructor
		SchmidtProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.SchmidtProjection"; }
		static std::string projection_type_str() {return "schmidt";}

		// projection and inverse projection
		virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		// purely regional? - no!
		bool isRegional() { return false; }	// schmidt is global grid
	
	private:
		double c_;		// stretching factor
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
