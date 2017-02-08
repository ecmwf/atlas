#ifndef atlas_grid_regular_RegionalLonLat_h
#define atlas_grid_regular_RegionalLonLat_h

#include "atlas/grid/regular/RegularRegional.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class RegionalLonLat: public RegularRegional {

public:

    static std::string grid_type_str();

    static std::string className();
    virtual std::string shortName() const;
    virtual std::string gridType() const { return "regional_lonlat"; }

    RegionalLonLat();
    RegionalLonLat(const util::Config& params);

protected:

    void setup(const util::Config& params);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);

    bool shiftLon_;
    bool shiftLat_;

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
