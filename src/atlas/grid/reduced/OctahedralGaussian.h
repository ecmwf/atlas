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
    
		std::vector<long> computePL(const size_t N);
		
    OctahedralGaussian(): ReducedGaussian() {};
    OctahedralGaussian(const util::Config& params);

    eckit::Properties spec() const;
        
  protected:

    void setup(size_t N);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);
    
};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
