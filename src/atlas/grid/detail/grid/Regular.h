#pragma once

#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Regular: public Structured {

public:

    static std::string static_type();

    static std::string className();

    virtual Spec spec() const;

    size_t nlon()           { return nlonmin(); }               // same for all latitudes
    double lon(size_t jlon) { return Structured::lon(0,jlon); } // same for all latitudes

protected:

    Regular();

};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
