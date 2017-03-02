#pragma once

#include "atlas/grid/detail/grid/regular/GlobalLonLat.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class ShiftedLat: public GlobalLonLat {

public:

    static std::string static_type();

    virtual std::string name() const;

    ShiftedLat( const Config& );
    ShiftedLat( long nlon, long nlat );


protected:

    void setup();

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
