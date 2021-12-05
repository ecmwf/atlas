/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>

#include "atlas/array/ArrayView.h"
#include "atlas/domain/detail/GlobalDomain.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildTorusXYZField.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildTorusXYZField::BuildTorusXYZField(const std::string& name): name_(name) {}

Field& BuildTorusXYZField::operator()(Mesh& mesh, const Domain& dom, double r0, double r1, idx_t nx, idx_t ny) const {
    return operator()(mesh.nodes(), dom, r0, r1, nx, ny);
}

Field& BuildTorusXYZField::operator()(mesh::Nodes& nodes, const Domain& dom, double r0, double r1, idx_t nx,
                                      idx_t ny) const {
    // fill xyz with torus coordinates. r0 and r1 are large and small radii,
    // respectively.

    auto domain = RectangularDomain(dom);
    ATLAS_ASSERT(domain);
    const double xmin = domain.xmin();
    const double xmax = domain.xmax();
    const double ymin = domain.ymin();
    const double ymax = domain.ymax();

    if (!nodes.has_field(name_)) {
        const idx_t npts                         = nodes.size();
        const array::ArrayView<double, 2> lonlat = array::make_view<double, 2>(nodes.xy());
        array::ArrayView<double, 2> xyz          = array::make_view<double, 2>(
            nodes.add(Field(name_, array::make_datatype<double>(), array::make_shape(npts, 3))));

        const double pi = M_PI;
        const double c1 = 2. * pi / double(nx) * (nx - 1) / (xmax - xmin);
        const double c2 = 2. * pi / double(ny) * (ny - 1) / (ymax - ymin);
        for (idx_t n = 0; n < npts; ++n) {
            double lon = -pi + c1 * (lonlat(n, 0) - xmin);
            double lat = -pi + c2 * (lonlat(n, 1) - ymin);

            xyz(n, 0) = std::cos(lon) * (r0 + r1 * std::cos(lat));
            xyz(n, 1) = std::sin(lon) * (r0 + r1 * std::cos(lat));
            xyz(n, 2) = r1 * std::sin(lat);
        }
    }
    return nodes.field(name_);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
