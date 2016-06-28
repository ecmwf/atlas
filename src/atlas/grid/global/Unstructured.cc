/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include "eckit/memory/Builder.h"
#include "atlas/grid/global/Unstructured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/ArrayView.h"


namespace atlas {
namespace grid {
namespace global {


Unstructured::Unstructured(const mesh::Mesh& m) :
    Global(),
    points_ ( new std::vector< Grid::Point > (m.nodes().size() ) ) {

    double lat_min = std::numeric_limits<double>::max();
    double lat_max = std::numeric_limits<double>::min();
    double lon_min = lat_min;
    double lon_max = lat_max;

    array::ArrayView<double,2> lonlat (m.nodes().lonlat());
    std::vector<Point> &p = *points_;
    const size_t npts = p.size();

    for( size_t n=0; n<npts; ++n) {
        p[n].assign(lonlat(n,internals::LON),lonlat(n,internals::LAT));
        lat_min = std::min( lat_min, p[n].lat() );
        lat_max = std::max( lat_max, p[n].lat() );
        lon_min = std::min( lon_min, p[n].lon() );
        lon_max = std::max( lon_max, p[n].lon() );
    }
}


Unstructured::Unstructured(const eckit::Parametrisation& p):
    Global() {
    NOTIMP;
}


Unstructured::Unstructured(std::vector<Point>* pts) :
    Global(),
    points_(pts) {

    const std::vector<Point> &p = *points_;
    const size_t npts = p.size();

    double lat_min = std::numeric_limits<double>::max();
    double lat_max = std::numeric_limits<double>::min();
    double lon_min = lat_min;
    double lon_max = lat_max;

    for (size_t n = 0; n < npts; ++n) {
        lat_min = std::min( lat_min, p[n].lat() );
        lat_max = std::max( lat_max, p[n].lat() );
        lon_min = std::min( lon_min, p[n].lon() );
        lon_max = std::max( lon_max, p[n].lon() );
    }
}


Unstructured::~Unstructured() {
}


Grid::uid_t Unstructured::shortName() const {
    if ( shortName_.empty() ) {
        std::ostringstream s;
        s <<  "unstructured." << Grid::hash().substr(0, 7);
        shortName_ = s.str();
    }
    return shortName_;
}


void Unstructured::hash(eckit::MD5 &md5) const {
    ASSERT(points_);

    const std::vector< Point > &pts = *points_;
    md5.add(&pts[0], sizeof(Point)*pts.size());

    for (size_t i = 0; i < pts.size(); i++) {
        const Point &p = pts[i];
        md5 << p.lon() << p.lat();
    }

}


size_t Unstructured::npts() const {
    ASSERT(points_);
    return points_->size();
}


void Unstructured::lonlat(std::vector<Grid::Point>& crds) const {
    ASSERT(points_);
    crds.resize(npts());
    for (size_t i = 0; i < npts(); ++i)
        crds[i].assign(
            (*points_)[i].lon(),
            (*points_)[i].lat() );
}


eckit::Properties Unstructured::spec() const {
    if (cached_spec_)
        return *cached_spec_;

    cached_spec_.reset( new eckit::Properties );

    cached_spec_->set("grid_type",gridType());

    std::vector<double> coords;
    coords.resize(2*npts());
    Grid::copyLonLatMemory(&coords[0],2*npts());

    cached_spec_->set( "lonlat", eckit::makeVectorValue<double>(coords) );

    return *cached_spec_;
}


void Unstructured::print(std::ostream& os) const {
    os << "Unstructured(Npts:" << npts() << ")";
}


register_BuilderT1(Grid, Unstructured, Unstructured::grid_type_str());


} // namespace global
} // namespace grid
} // namespace atlas
