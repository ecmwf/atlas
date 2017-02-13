#include "atlas/grid/regular/RegularRegional.h"
#include "atlas/grid/domain/RectangularDomain.h"

namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularRegional,RegularRegional::grid_type_str());

std::string RegularRegional::grid_type_str() {
    return "regular_regional";
}

std::string RegularRegional::className() {
    return "atlas.grid.regular.RegularRegional";
}

std::string RegularRegional::shortName() const {
    std::ostringstream s;
    s << "RR"<< nlonmin() << "x" << nlat();
    return s.str();
}

void RegularRegional::setup(const util::Config& config) {

    util::Config config_dom, config_proj;

    long nx, ny;
    std::vector<double> bbox(4);
    std::vector<double> sw(2), ne(2), center(2);
    double dx, dy;

    // read subconfigurations

    // projection
    if ( config.get("projection",config_proj) ) {
    } else {
      // default, error, or hardcoded default?
      config_proj.set("type","lonlat");
    }
    projection_.reset( projection::Projection::create(config_proj) );

    // dimensions
    if ( ! config.get("nx",nx) ) throw eckit::BadParameter("nx missing in Params",Here());
    if ( ! config.get("ny",ny) ) throw eckit::BadParameter("ny missing in Params",Here());

    // domain
    if ( config.get("domain",config_dom) ) {
      // domain is specified either by bbox, by sw and ne, or by center and resolution
      if ( !config_dom.get("bounding_box",bbox) ) {
        if ( config_dom.get("center",center) && config_dom.get("dx",dx) && config_dom.get("dy",dy)  ) {
          // coordinates of center
          double crd_center[] = {center[0],center[1]};
          projection_->lonlat2coords(crd_center);
          // calculate bbox
          double x = crd_center[0];
          double y = crd_center[1];
          double LX = (nx-1)*dx;
          double LY = (ny-1)*dy;
          bbox[0]=x-0.5*LX; bbox[1]=x+0.5*LX;
          bbox[2]=y-0.5*LY; bbox[3]=y+0.5*LY;
        } else {
          if ( config_dom.get("sw",sw) && config_dom.get("ne",ne) )  {
            // corner coordinates
            double crd_sw[] = {sw[0],sw[1]};
            double crd_ne[] = {ne[0],ne[1]};
            // sw corner
            projection_->lonlat2coords(crd_sw);
            // put in bbox
            bbox[0]=crd_sw[0]; bbox[2]=crd_sw[1];
            // ne corner
            projection_->lonlat2coords(crd_ne);
            // put in bbox
            bbox[1]=crd_ne[0]; bbox[3]=crd_ne[1];
          } else {
            throw eckit::BadParameter("RegularRegional grid domain should be specified by (i) bbox, (ii) center, dx and dy, or (iii) ne and sw",Here());
          }
        }
        // put bbox in config_dom
        config_dom.set("bounding_box",bbox);
      }
      // set domainType if it's missing
      std::string domainType;
      if (!config_dom.get("type",domainType)) config_dom.set("type","rectangular");
      // create domain from configuration
      domain_.reset( domain::Domain::create(config_dom) );
      // check if the domain is rectangular
      if( not dynamic_cast<domain::RectangularDomain*>(domain_.get()) ) {
        throw eckit::BadParameter("RegularRegional grid requires a RectangularDomain",Here());
      }

    } else {
      // error
      throw eckit::BadParameter("domain is required for a RegularRegional grid",Here());
    }

    // spacing_x
    {
      util::Config config_spacing;
      std::string spacingType;

      if ( !config.get("spacing_x",spacingType) ) spacingType="uniform";
      // set configuration of spacing_x
      config_spacing.set("spacingType",spacingType);
      config_spacing.set("xmin",bbox[0]);
      config_spacing.set("xmax",bbox[1]);
      config_spacing.set("N",nx);
      // create spacing
      spacing_x_=spacing::Spacing::create(config_spacing);
    }

    // spacing_y
    {
      util::Config config_spacing;
      std::string spacingType;

      if ( !config.get("spacing_y",spacingType) ) spacingType="uniform";
      // set configuration of spacing_y
      config_spacing.set("spacingType",spacingType);
      config_spacing.set("xmin",bbox[2]);
      config_spacing.set("xmax",bbox[3]);
      config_spacing.set("N",ny);
      // create spacing
      spacing_y_=spacing::Spacing::create(config_spacing);
    }

    // setup regular grid
    Regular::setup();

}

RegularRegional::RegularRegional(const util::Config& config) :
    Regular()
{
    setup(config);
}

RegularRegional::RegularRegional() :
    Regular()
{
}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

