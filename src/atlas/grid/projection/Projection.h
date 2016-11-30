#ifndef atlas_Projection_H
#define atlas_Projection_H

#include "eckit/geometry/KPoint.h"
#include "eckit/geometry/Point2.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace projection {

class Projection {

	public:
		typedef const eckit::Parametrisation& ARG1;
		typedef eckit::BuilderT1<Projection> builder_t;

	public:
	
		static Projection* create();
		static Projection* create(const eckit::Parametrisation& p);
		
		/* 
		// implementations moved to Projection.cc
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
		*/
		
		static std::string className() {return "atlas.Projection";}
		static std::string projection_type_str() {return "projection";}
		
		// purely virtual functions: must be implemented by inheriting classes
		virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2)=0;
		virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2)=0;
		virtual bool isRegional()=0;
	
	protected:
		void rotate_(eckit::geometry::LLPoint2 &P,const eckit::geometry::LLPoint2 &pole);			// coordinates of the point on a rotated sphere with specified pole
		void unrotate_(eckit::geometry::LLPoint2 &P,const eckit::geometry::LLPoint2 &pole);		// inverse operation of rotate

};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
