/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/util/Config.h"
#include "eckit/memory/Builder.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/internals/Parameters.h"


namespace atlas {
namespace grid {
namespace detail {
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

    virtual std::string shortName() const { return "structured"; }
    virtual std::string gridType() const { return "structured"; }


    CustomStructured( const Config& );

    CustomStructured(
        size_t nlat,
        const double lats[],
        const long pl[]
    );

};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
