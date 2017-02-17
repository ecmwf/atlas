#ifndef atlas_grid_regular_RegularRegional_h
#define atlas_grid_regular_RegularRegional_h

#include "atlas/grid/regular/Regular.h"
#include "atlas/grid/domain/Domain.h"
#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace regular {

class RegularRegional: public Regular {

public:

    static std::string grid_type_str();

    static std::string className();

    virtual std::string shortName() const;
    virtual std::string gridType() const { return "regular_regional"; }

    RegularRegional(const util::Config& params);
    RegularRegional();

protected:

    struct Parse {
      bool valid;
      operator bool() const { return valid; }
      virtual void apply( RegularRegional& ) const =0;
    };

    struct ParseUniformCentred : Parse {
      ParseUniformCentred(const eckit::Parametrisation&);
      long nx;
      long ny;
      double dx;
      double dy;
      std::vector<double> centre_lonlat;
      bool endpoint_x;
      bool endpoint_y;
      virtual void apply( RegularRegional& ) const;
    };

    struct ParseBounds : Parse {
      ParseBounds(const eckit::Parametrisation&);
      long nx;
      long ny;
      double xmin;
      double xmax;
      double ymin;
      double ymax;
      bool endpoint_x;
      bool endpoint_y;
      virtual void apply( RegularRegional& ) const;
    };

    struct ParseLonLatBounds : Parse {
      ParseLonLatBounds(const eckit::Parametrisation&);
      long nx;
      long ny;
      double north;
      double west;
      double south;
      double east;
      bool endpoint_x;
      bool endpoint_y;
      virtual void apply( RegularRegional& ) const;
    };

    void setup(const util::Config& params);
    std::string shortName();

};


}  // namespace regular
}  // namespace grid
}  // namespace atlas


#endif
