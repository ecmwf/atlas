/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grid_CustomStructured_h
#define atlas_grid_CustomStructured_h

#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/grid/Structured.h"
#include "atlas/internals/Parameters.h"


namespace atlas {
namespace grid {


/**
 * @brief CustomStructured Grid
 *
 * This class is a base class for all grids that can be described by
 * constant latitudes with a uniform distribution of points per latitude
 * in zonal direction.
 * This means any full grid and reduced grid, both regular, gaussian or other
 * such distribution can be represented with this class
 */
class CustomStructured: public Structured {
  public:

    static std::string className();

    static std::string grid_type_str();

    CustomStructured(const eckit::Parametrisation&);

    CustomStructured(
        size_t nlat,
        const double lats[],
        const long pl[]
    );

    CustomStructured(
        size_t nlat,
        const double lats[],
        const long pl[],
        const double lonmin[],
        const double lonmax[]
    );

    virtual eckit::Properties spec() const;

/*
// domain is now at the level of Grid.
    virtual const domain::Domain& domain() const {
        return domain_;
    }

  private:

    domain::Domain domain_;
 */
};


}  // namespace grid
}  // namespace atlas


#endif
