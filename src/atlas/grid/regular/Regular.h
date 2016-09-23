#ifndef atlas_grid_regular_Regular_h
#define atlas_grid_regular_Regular_h

#include "atlas/grid/Structured.h"


namespace atlas {
namespace grid {
namespace regular {

class Regular: public Structured {

  public:

    static std::string grid_type_str();

    static std::string className();

    Regular(const util::Config& params);
    Regular();

    virtual const domain::Domain& domain() const { return *domain_; }
    
    virtual eckit::Properties spec() const;
    
  protected:

    void setup();

    //virtual void set_typeinfo() = 0;

    //static eckit::Value domain_spec(const Domain& dom);

  protected:

    domain::Domain * domain_;
    
    spacing::Spacing * spacing_x_;
    spacing::Spacing * spacing_y_;
    
};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
