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
