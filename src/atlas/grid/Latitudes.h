/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Peter Bispham
/// @author Tiago Quintino
/// @date 2013

#ifndef eckit_grid_Coordinates_H
#define eckit_grid_Coordinates_H

#include <vector>

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

class Latitudes {

public: // methods


    /// generates positive gaussian latitudes N hemi only
    static void gaussian(size_t resolution, std::vector<double>& lats);

    /// generates positive latitudes equally spaced in N hemi only
    static void uniform(size_t resolution, std::vector<double>& lats);

private:

    // Generates latitudes in N hemisphere 
    static void initialGaussianLatitudes(size_t resolution, std::vector<double>& lats);
    
    // Only tested on in N hemisphere (positive) values and we mirror
    // result to southern hemisphere
    static void refineGaussianLatitude(size_t resolution, double& value);
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
