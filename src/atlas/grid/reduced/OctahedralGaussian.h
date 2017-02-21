#ifndef atlas_grid_reduced_OctahedralGaussian_h
#define atlas_grid_reduced_OctahedralGaussian_h

#include "atlas/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace reduced {

class OctahedralGaussian: public ReducedGaussian {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "octahedral_gaussian"; }

    static std::vector<long> computePL(const size_t N, const size_t start=20);

    OctahedralGaussian(const util::Config& params);

    virtual eckit::Properties spec() const;

protected:

    void setup( const size_t N, const size_t start=20 );

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
