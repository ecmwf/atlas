/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/field/FieldBuilder.h"
#include "atlas/functionspace/StructuredColumns.h"

namespace atlas {
namespace field {
namespace builders {

/**
 * \brief Build a field containing a representative grid-box area in square metres.
 *
 * The source object must be a FunctionSpace backed by
 * functionspace::StructuredColumns. The underlying projection must either be
 * `lonlat` or `rotated_lonlat`, with coordinates in degrees, or report its
 * coordinate units as metres. Other function spaces and projection units are
 * rejected.
 *
 * Grid boxes are bounded by the midpoints between neighbouring row and column
 * coordinates. The algorithm assumes that the longitudinal or projected
 * x-spacing is constant within each structured row. Consequently, one area is
 * computed per row and assigned to every local point in that row, including
 * halo points. Halo rows use StructuredColumns::compute_j() to select the
 * corresponding internal row spacing and StructuredColumns::compute_xy() to
 * obtain reflected boundary coordinates.
 *
 * For a spherical lon/lat grid, the area between longitudes
 * \f$\lambda_W,\lambda_E\f$ and latitudes \f$\phi_S,\phi_N\f$ is evaluated by
 * integrating the spherical surface element:
 *
 * \f[
 * A = \int_{\lambda_W}^{\lambda_E}\int_{\phi_S}^{\phi_N}
 *     R^2\cos(\phi)\,\mathrm{d}\phi\,\mathrm{d}\lambda
 *   = R^2\left|\lambda_E-\lambda_W\right|
 *     \left|\sin(\phi_N)-\sin(\phi_S)\right|.
 * \f]
 *
 * Angular quantities are converted from degrees to radians. A strip crossing
 * either pole is split at \f$\phi=\pm\frac{\pi}{2}\f$ so mirrored latitude
 * contributions do not cancel.
 *
 * For projections expressed in metres, each grid box is treated as a planar
 * rectangle:
 *
 * \f[
 * A = \left|x_E-x_W\right|\,\left|y_N-y_S\right|
 *   = \Delta x\,\left|y_N-y_S\right|.
 * \f]
 *
 * The result is a rank-one double-precision Field with no vertical levels.
 * Its name defaults to `area` and can be changed with `field_name`. The sphere
 * radius defaults to geometry::Earth().radius() and can be set with `radius`.
 * If `geometry` is supplied, its radius takes precedence over `radius`.
 */
class GridBoxArea : public AbstractFieldBuilder {
public:
    GridBoxArea(const util::Config& config, const util::ObjectHandleBase& object);

    Field operator()() const override;

private:
    functionspace::StructuredColumns fs_;
    std::string field_name_;
    double radius_;
};

}  // namespace builders
}  // namespace field
}  // namespace atlas