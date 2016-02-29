/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Jan 2014

#ifndef atlas_grids_rgg_OctahedralReducedGaussianGrid_h
#define atlas_grids_rgg_OctahedralReducedGaussianGrid_h

#include "eckit/memory/Builder.h"
#include "atlas/grid/ReducedGaussianGrid.h"

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

class OctahedralReducedGaussianGrid : public ReducedGaussianGrid {
public:

    static std::string className() { return "atlas.grids.rgg.OctahedralReducedGaussianGrid"; }
    static std::string grid_type_str() { return "oct"; }

    OctahedralReducedGaussianGrid(const size_t N, const size_t octahedralPoleStart = 20);

    OctahedralReducedGaussianGrid( const eckit::Parametrisation& arg1);

    /// Computes the PL for the Octohedral distribution
    /// number of points at latitude closest to pole
    static std::vector<long> computePL(const size_t N, const size_t start);

private:

    void construct(const size_t N, const size_t start);

    void set_typeinfo();

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif // atlas_grids_rgg_OctahedralReducedGaussianGrid_h
