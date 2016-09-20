#ifndef atlas_SchmidtProjection_H
#define atlas_SchmidtProjection_H

#include "atlas/grid/projections/Projection.h"

namespace atlas {
namespace grid {

class SchmidtProjection: public Projection {
	public:
		// constructor
		SchmidtProjection(const eckit::Parametrisation& p);

		// class name
		static std::string className() { return "atlas.SchmidtProjection"; }

		// projection and inverse projection
		eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2);
		eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2);

		// purely regional? - no!
		bool isRegional() { return false; }	// schmidt is global grid
	
	private:
		double c_;		// stretching factor
};

}  // namespace grid
}  // namespace atlas


#endif
