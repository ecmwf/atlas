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

		std::vector<long> computePL(const size_t N);
		
    OctahedralGaussian(): ReducedGaussian() {};
    OctahedralGaussian(const util::Config& params);
    
    virtual eckit::Properties spec() const;

  protected:

    void setup(size_t N);

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
