#ifndef atlas_Spacing_H
#define atlas_Spacing_H

#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {

class Spacing {

	public:
		typedef const eckit::Parametrisation& ARG1;
		typedef eckit::BuilderT1<Spacing> builder_t;

	public:
	
		static Spacing* create() {
			// default: uniform spacing
			util::Config params;
			params.set("spacingType","atlas.UniformSpacing");
			return Spacing::create(params);
		};
		
		static Spacing* create(const eckit::Parametrisation& params) {
			std::string spacingType;
			if (params.get("spacingType",spacingType)) {
				return eckit::Factory<Spacing>::instance().get(spacingType).create(params);
			}

			// should return error here
	    throw eckit::BadParameter("spacingType missing in params",Here());
			return NULL;
		}
		
		static std::string className() {return "atlas.Spacing";}
		
		// purely virtual functions: must be implemented by inheriting classes
		virtual std::vector<double>         generate(double xmin, double xmax, size_t N)=0;

};

}  // namespace grid
}  // namespace atlas


#endif
