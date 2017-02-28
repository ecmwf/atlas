#pragma once

#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class ClassicGaussian: public ReducedGaussian {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "classic_gaussian"; }

    ClassicGaussian(): ReducedGaussian() {}
    ClassicGaussian( const Config& );
    ClassicGaussian( size_t N );

protected:

    void setup(size_t N);

 };


}  // namespace reduced
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
