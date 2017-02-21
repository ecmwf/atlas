/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_method_FiniteElement_h
#define atlas_interpolation_method_FiniteElement_h

#include "atlas/interpolation/method/Method.h"

#include <string>
#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"
#include "atlas/array/ArrayView.h"
#include "atlas/interpolation/method/PointIndex3.h"
#include "atlas/mesh/Elements.h"


namespace atlas {
namespace interpolation {
namespace method {


class FiniteElement : public Method {
public:

    FiniteElement(const Config& config) : Method(config) {}
    virtual ~FiniteElement() {}

    /**
     * @brief Create an interpolant sparse matrix relating two (pre-partitioned) meshes,
     * using elements as per the Finite Element Method and ray-tracing to calculate
     * source mesh elements intersections (and interpolation weights) with target grid
     * node-containing rays
     * @param meshSource mesh containing source elements
     * @param meshTarget mesh containing target points
     */
    void setup(mesh::Mesh& meshSource, mesh::Mesh& meshTarget);

protected:

    /**
     * Find in which element the point is contained by projecting (ray-tracing) the
     * point to the nearest element(s), returning the (normalized) interpolation weights
     */
    static Triplets projectPointToElements(
            const array::ArrayView<double, 2>& icoords,
            const array::ArrayView<double, 2>& ilonlat,
            const mesh::Connectivity& connectivity,
            const Point &p,
            size_t ip,
            ElemIndex3::NodeList::const_iterator start,
            ElemIndex3::NodeList::const_iterator finish,
            std::ostream& failures_log );

};


}  // method
}  // interpolation
}  // atlas


#endif
