#ifndef atlas_Spacing_H
#define atlas_Spacing_H

#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace spacing {

class Spacing {

	public:
		typedef const eckit::Parametrisation& ARG1;
		typedef eckit::BuilderT1<Spacing> builder_t;

	public:
	
		static Spacing* create() {
			// default: uniform spacing
			util::Config params;
			params.set("spacingType","uniform");
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
		static std::string spacing_type_str() {return "spacing";}
		
		// purely virtual functions: must be implemented by inheriting classes
		virtual void generate(double xmin, double xmax, size_t N, std::vector<double>& x) const =0;

};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas


#endif
