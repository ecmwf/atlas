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
    
    virtual std::string shortName() const { return "regular"; }

    Regular(const util::Config& params);
    Regular();

    virtual eckit::Properties spec() const;
    
  protected:

    void setup();

  protected:

    spacing::Spacing * spacing_x_;
    spacing::Spacing * spacing_y_;
    
};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
