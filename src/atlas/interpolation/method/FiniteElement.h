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

    FiniteElement(const Config& config) :
        Method(config),
        fallback_to_2d_(false) {
        config.get("fallback_to_2d",fallback_to_2d_);
    }

    virtual ~FiniteElement() {}

    virtual void setup(const FunctionSpace& source, const FunctionSpace& target) override;

protected:

    /**
     * @brief Create an interpolant sparse matrix relating two (pre-partitioned) meshes,
     * using elements as per the Finite Element Method and ray-tracing to calculate
     * source mesh elements intersections (and interpolation weights) with target grid
     * node-containing rays
     * @param meshSource mesh containing source elements
     * @param meshTarget mesh containing target points
     */
    //virtual void setup(Mesh& meshSource, Mesh& meshTarget) override;

    /**
     * Find in which element the point is contained by projecting (ray-tracing) the
     * point to the nearest element(s), returning the (normalized) interpolation weights
     */
    Triplets projectPointToElements(
            size_t ip,
            const ElemIndex3::NodeList& elems,
            std::ostream& failures_log ) const;


protected:

    mesh::MultiBlockConnectivity* connectivity_;
    std::unique_ptr<array::ArrayView<double,2>> icoords_;
    std::unique_ptr<array::ArrayView<double,2>> ilonlat_;
    std::unique_ptr<array::ArrayView<double,2>> ocoords_;
    std::unique_ptr<array::ArrayView<double,2>> olonlat_;
    bool fallback_to_2d_;
};


}  // method
}  // interpolation
}  // atlas


#endif
