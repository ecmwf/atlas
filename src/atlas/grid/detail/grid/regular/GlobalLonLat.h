#pragma once

#include "atlas/grid/detail/grid/Regular.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {


struct Shift {

    enum Bits {
        NONE = 0,
        LAT  = (1<<1),
        LON  = (1<<2)
    };

    Shift(int bits=NONE) : bits_(bits) {
    }

    Shift(bool shift_lon, bool shift_lat) : bits_((shift_lon? LON:NONE) | (shift_lat? LAT:NONE)) {
    }

    bool operator()(int bits) const {
        return (bits_ & bits) == bits;
    }

    const int bits_;

};

// Uninstantiatable class
class GlobalLonLat: public Regular {

public:

    bool isShiftedLon() const { return shifted_x_; }
    bool isShiftedLat() const { return shifted_y_; }

protected:

    GlobalLonLat();

    void setup(long nx, long ny, Shift);

    virtual eckit::Properties spec() const;

    bool shifted_x_;
    bool shifted_y_;
};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
