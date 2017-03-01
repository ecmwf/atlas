#pragma once

#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class OctahedralGaussian: public ReducedGaussian {

public:

    static std::string static_type();

    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

    static std::vector<long> computePL(const size_t N, const size_t start=20);

    OctahedralGaussian( const Config& );

    virtual Spec spec() const;

protected:

    void setup( const size_t N, const size_t start=20 );

};


}  // namespace reduced
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
