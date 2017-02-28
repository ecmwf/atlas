#pragma once

#include "atlas/grid/detail/grid/Regular.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

class RegularGaussian: public Regular {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "regular_gaussian"; }

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
