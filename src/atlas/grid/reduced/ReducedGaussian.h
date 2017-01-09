#ifndef atlas_grid_reduced_ReducedGaussian_h
#define atlas_grid_reduced_ReducedGaussian_h

#include "atlas/grid/Structured.h"

namespace atlas {
namespace grid {
namespace reduced {

class ReducedGaussian: public Structured {

  public:

    static std::string grid_type_str();

    static std::string className();
    
    virtual std::string shortName() const;
    virtual std::string gridType() const { return "reduced_gaussian"; }

    ReducedGaussian(): Structured() {};
    ReducedGaussian(const util::Config& params);
    ReducedGaussian(const int N, const long pl[]);

    virtual eckit::Properties spec() const;
    
  protected:

    void setup(const size_t N, const long pl[]);

    //virtual void set_typeinfo() = 0;
    //static eckit::Value domain_spec(const Domain& dom);
    
};


}  // namespace reduced
}  // namespace grid
}  // namespace atlas


#endif
