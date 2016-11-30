#ifndef atlas_grid_reduced_ClassicGaussian_h
#define atlas_grid_reduced_ClassicGaussian_h

#include "atlas/grid/Structured.h"

namespace atlas {
namespace grid {
namespace reduced {

class ClassicGaussian: public Structured {

  public:

    static std::string grid_type_str();

    static std::string className();
    
    ClassicGaussian(): Structured() {};
    ClassicGaussian(const util::Config& params);

    eckit::Properties spec() const;
    
    virtual const domain::Domain * domain_ptr() const { return domain_; }
    virtual const domain::Domain &domain() const { return *domain_; }
    
    
  protected:

    void setup(size_t N);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);
    
    domain::Domain * domain_;

};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
