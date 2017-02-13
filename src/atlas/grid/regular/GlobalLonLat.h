#pragma once

#include "atlas/grid/regular/Regular.h"

namespace atlas {
namespace grid {
namespace regular {

// Uninstantiatable class
class GlobalLonLat: public Regular {

protected:

    GlobalLonLat();

    bool isShiftedLon() const { return shiftLon_; }
    bool isShiftedLat() const { return shiftLat_; }

protected:

    void setup(long nlon, long nlat);

    bool shiftLon_;
    bool shiftLat_;

    virtual eckit::Properties spec() const;

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas
