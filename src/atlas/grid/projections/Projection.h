#ifndef atlas_Projection_H
#define atlas_Projection_H

#include "eckit/geometry/KPoint.h"
#include "eckit/geometry/Point2.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {

class Projection {

	public:
		typedef const eckit::Parametrisation& ARG1;
		typedef eckit::BuilderT1<Projection> builder_t;

	public:
	
		static Projection* create() {
			// default: no projection, i.e. stay in (lon,lat)-space
			util::Config projParams;
			projParams.set("projectionType","atlas.LonLatProjection");
			return Projection::create(projParams);
		};
		
		static Projection* create(const eckit::Parametrisation& p) {
			std::string projectionType;
			if (p.get("projectionType",projectionType)) {
				return eckit::Factory<Projection>::instance().get(projectionType).create(p);
			}

			// should return error here
	    throw eckit::BadParameter("projectionType missing in Params",Here());
			return NULL;
		}
		
		static std::string className() {return "atlas.Projection";}
		
		// purely virtual functions: must be implemented by inheriting classes
		virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2)=0;
		virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2)=0;
		virtual bool isRegional()=0;
};

}  // namespace grid
}  // namespace atlas


#endif
