#pragma once

#include "atlas/grid/detail/grid/regular/GlobalLonLat.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class ShiftedLonLat: public GlobalLonLat {

public:

    static std::string static_type();

    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

    ShiftedLonLat( const Config& );
    ShiftedLonLat( long nlon, long nlat );

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
