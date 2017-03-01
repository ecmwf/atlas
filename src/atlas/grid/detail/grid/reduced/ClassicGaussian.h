#pragma once

#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class ClassicGaussian: public ReducedGaussian {

public:

    static std::string static_type();

    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

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
