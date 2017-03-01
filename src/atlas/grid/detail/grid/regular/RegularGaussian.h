#pragma once

#include "atlas/grid/detail/grid/Regular.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class RegularGaussian: public Regular {

public:

    static std::string static_type();

    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

    RegularGaussian( const Config& );
    RegularGaussian( long N );

protected:

    void setup(long N);

};


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
