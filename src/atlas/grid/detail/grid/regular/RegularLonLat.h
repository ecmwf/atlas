#pragma once

#include "atlas/grid/detail/grid/regular/GlobalLonLat.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class RegularLonLat: public GlobalLonLat {

public:

    static std::string static_type();

    virtual std::string name() const;

    RegularLonLat( const Config& );
    RegularLonLat( long nlon, long nlat );
    RegularLonLat( long N );

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
