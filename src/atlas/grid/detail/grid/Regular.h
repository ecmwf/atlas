#pragma once

#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Regular: public Structured {

public:

    static std::string static_type();

    virtual Spec spec() const;

    size_t nx()        const { return nxmin(); }            // same for every y
    double x(size_t i) const { return Structured::x(i,0); } // same for every y

protected:

    Regular();

};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
