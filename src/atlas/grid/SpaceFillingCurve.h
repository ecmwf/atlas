/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>

#include "atlas/domain/Domain.h"
#include "atlas/library/config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {

// -------------------------------------------------------------------------------------

/// @brief Class to compute a global index given a coordinate, based on the
/// Hilbert Spacefilling Curve.
///
/// This algorithm is based on:
/// - John J. Bartholdi and Paul Goldsman "Vertex-Labeling Algorithms for the Hilbert Spacefilling Curve"\n
/// It is adapted to return contiguous numbers of the gidx_t type, instead of a double [0,1]
///
/// Given a bounding box and number of hilbert recursions, the bounding box can be divided in
/// 2^(dim*levels) equally spaced cells. A given coordinate falling inside one of these cells, is assigned
/// the 1-dimensional Hilbert-index of this cell. To make sure that 1 coordinate corresponds to only 1
/// Hilbert index, the number of levels have to be increased.
/// In 2D, the recursion cannot be higher than 15, if you want the indices to fit in "unsigned int" type of 32bit.
/// In 2D, the recursion cannot be higher than 30, if you want the indices to fit in "unsigned int" type of 64bit.
///
///
/// No attempt is made to provide the most efficient algorithm. There exist other open-source
/// libraries with more efficient algorithms, such as libhilbert, but its LGPL license
/// is not compatible with this licence.
///
/// @author Willem Deconinck
class HilbertCurve {
public:
    /// Constructor
    /// Initializes the hilbert space filling curve with a given "space" and "levels"
    HilbertCurve(const Domain& domain, idx_t levels);

    /// Compute the hilbert code for a given point in 2D
    gidx_t operator()(const PointXY& point) const;

    /// Return the maximum hilbert code possible with the initialized levels
    ///
    /// Care has to be taken that this number is not larger than the precision of the type storing
    /// the hilbert codes.
    gidx_t nb_keys() const { return nb_keys_; }

private:  // functions
    using box_t = std::array<PointXY, 4>;

    /// @brief Recursive algorithm
    gidx_t recursive_algorithm(const PointXY& p, const box_t& box, idx_t level) const;

private:  // data
    /// Vertex label type (4 vertices in 2D)
    enum VertexLabel
    {
        A = 0,
        B = 1,
        C = 2,
        D = 3
    };

    /// Bounding box, defining the space to be filled
    const RectangularDomain domain_;

    /// maximum recursion level of the Hilbert space filling curve
    idx_t max_level_;

    /// maximum number of unique codes, computed by max_level
    gidx_t nb_keys_;
    gidx_t nb_keys_2_;
};

// -------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
