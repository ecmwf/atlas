#ifndef atlas_grid_regular_RegularLonLat_h
#define atlas_grid_regular_RegularLonLat_h

#include "atlas/grid/regular/GlobalLonLat.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class RegularLonLat: public GlobalLonLat {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "regular_lonlat"; }

    RegularLonLat(const util::Config& params);
    RegularLonLat(long nlon, long nlat);
    RegularLonLat(long N);

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
