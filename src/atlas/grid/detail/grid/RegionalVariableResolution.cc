/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "RegionalVariableResolution.h"

#include "eckit/types/FloatCompare.h"

#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/NormaliseLongitude.h"

using atlas::grid::LinearSpacing;
using XSpace = atlas::StructuredGrid::XSpace;
using YSpace = atlas::StructuredGrid::YSpace;

namespace atlas {
namespace grid {
namespace {  // anonymous

static class regional_var_resolution : public GridBuilder {
public:
    regional_var_resolution(): GridBuilder("regional_variable_resolution") {}

    void print(std::ostream&) const override {
        // os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian
        // grid";
    }


    const Grid::Implementation* create(const std::string& /*name*/, const Grid::Config&) const override {
        throw_NotImplemented("There are no named regional_var_resolution grids implemented.", Here());
        return nullptr;
    }

    //create return pointer, data type to return
    const Grid::Implementation* create(const Grid::Config& config) const override {
        // read projection subconfiguration
        double inner_xmin  = 0;
        double inner_xmax  = 0;
        double inner_ymin  = 0;
        double inner_ymax  = 0;
        double outer_xmin  = 0;
        double outer_xmax  = 0;
        double outer_ymin  = 0;
        double outer_ymax  = 0;
        double delta_inner = 0.;
        Projection projection;
        {
            util::Config config_all, config_proj, config_inner, config_outer;
            util::Config config_pr, configwx, configwy;

            config_all.set("type", "variable_resolution");

            if (config.get("projection", config_proj)) {
                config_all.set(config_proj);

                if (config_proj.getString("type") == "lonlat") {
                    config_all.set("type", "variable_resolution");
                }
                else if (config_proj.getString("type") == "rotated_lonlat") {
                    config_all.set("type", "rotated_variable_resolution");
                }
                else {
                    throw_Exception("Only \"lonlat\" or \"rotated_lonlat\" projection is supported");
                }
            }

            config.get("inner", config_inner);

            config.get("outer", config_outer);

            // merge of multiple configs use |
            config.get("progression", config_pr);
            config_all.set("progression", config_pr);

            config_all.set("inner", config_inner);
            config_all.set("outer", config_outer);
            projection = Projection(config_all);
        }

        //< Add different parts in the grid something like
        //< int inner_xmin = config.getInt("inner.xmin");

        //compute number of points in the regular grid
        util::Config config_grid;

        config.get("inner.xmin", inner_xmin);
        config.get("inner.xend", inner_xmax);
        config.get("inner.ymin", inner_ymin);
        config.get("inner.yend", inner_ymax);
        config.get("outer.xmin", outer_xmin);
        config.get("outer.xend", outer_xmax);
        config.get("outer.ymin", outer_ymin);
        config.get("outer.yend", outer_ymax);
        config.get("inner.dx", delta_inner);

        constexpr float epstest = std::numeric_limits<float>::epsilon();
        int nx_reg              = ((outer_xmax - outer_xmin + epstest) / delta_inner) + 1;
        int ny_reg              = ((outer_ymax - outer_ymin + epstest) / delta_inner) + 1;

        YSpace yspace = LinearSpacing{outer_ymin, outer_ymax, ny_reg};
        XSpace xspace = LinearSpacing{outer_xmin, outer_xmax, nx_reg};

        //< RegularGrid is a type of structuredGrid
        auto domain_ = RectangularDomain{{outer_xmin, outer_xmax}, {outer_ymin, outer_ymax}};
        //< allocate memory to make class, create an object using new "constructor"
        return new StructuredGrid::grid_t(xspace, yspace, projection, domain_);
    }

    void force_link() {}

} regional_var_resolution_;

}  // namespace

namespace detail {
namespace grid {

void force_link_Regional_var_resolution() {
    regional_var_resolution_.force_link();
}

}  // namespace grid
}  // namespace detail

}  // namespace grid
}  // namespace atlas
