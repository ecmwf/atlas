#ifndef atlas_grid_reduced_ClassicGaussian_h
#define atlas_grid_reduced_ClassicGaussian_h

#include "atlas/grid/reduced/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace reduced {

class ClassicGaussian: public ReducedGaussian {

  public:

    static std::string grid_type_str();

    static std::string className();
    
    virtual std::string shortName() const;
    virtual std::string gridType() const { return "classic_gaussian"; }

    ClassicGaussian(): ReducedGaussian() {};
    ClassicGaussian(const util::Config& params);
    ClassicGaussian(size_t N);
   
  protected:

    void setup(size_t N);

 };


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
