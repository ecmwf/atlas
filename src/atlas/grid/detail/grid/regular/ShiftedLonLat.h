#pragma once

#include "atlas/grid/detail/grid/regular/GlobalLonLat.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class ShiftedLonLat: public GlobalLonLat {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "shifted_lonlat"; }

    ShiftedLonLat( const Config& );
    ShiftedLonLat( long nlon, long nlat );

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
