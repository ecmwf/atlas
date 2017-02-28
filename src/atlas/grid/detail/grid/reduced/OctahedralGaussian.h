#pragma once

#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

class OctahedralGaussian: public ReducedGaussian {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "octahedral_gaussian"; }

    static std::vector<long> computePL(const size_t N, const size_t start=20);

    OctahedralGaussian( const Config& );

    virtual eckit::Properties spec() const;

protected:

    void setup( const size_t N, const size_t start=20 );

};


}  // namespace reduced
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
