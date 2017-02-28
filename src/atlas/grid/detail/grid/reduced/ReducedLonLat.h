#pragma once

#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class ReducedLonLat: public Structured {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "reduced_lonlat"; }

    ReducedLonLat( const Config& );

    Spec spec() const;

protected:

    void setup(size_t nlat, long pl[]);

};


}  // namespace reduced
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
