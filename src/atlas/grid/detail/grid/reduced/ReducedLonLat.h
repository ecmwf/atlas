#pragma once

#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class ReducedLonLat: public Structured {

public:

    static std::string static_type();

    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

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
