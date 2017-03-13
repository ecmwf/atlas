
/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/ArrayView.h"
#include "atlas/interpolation/method/PointIndex3.h"
#include "atlas/mesh/HybridElements.h"


namespace atlas {
namespace interpolation {
namespace method {


ElemIndex3* create_element_centre_index(const mesh::Mesh& mesh) {

    const array::ArrayView< double, 2 > centres(mesh.cells().field( "centre" ));

    std::vector< ElemIndex3::Value > p;
    p.reserve(mesh.cells().size());

    for (size_t j = 0; j < mesh.cells().size(); ++j) {
        p.push_back(ElemIndex3::Value(
                        ElemIndex3::Point(centres(j, XX), centres(j, YY), centres(j, ZZ)),
                        ElemIndex3::Payload(j) ));
    }


    ElemIndex3* tree = new ElemIndex3();
    tree->build(p.begin(), p.end());

    return tree;
}


} // namespace method
} // namespace interpolation
} // namespace atlas

