/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/gaussian/OctahedralGaussian.h"


namespace atlas {
namespace grid {
namespace gaussian {


std::string OctahedralGaussian::className() {
    return "atlas.grid.gaussian.OctahedralGaussian";
}


std::string OctahedralGaussian::grid_type_str() {
    return "octahedral_gaussian";
}


std::vector<long> OctahedralGaussian::computePL(const size_t N) {
    const size_t start = 20;
    std::vector<long> pl(N);
    for(size_t jlat=0; jlat < N; ++jlat) {
        pl[jlat] = start + 4*jlat;
    }
    return pl;
}


OctahedralGaussian::OctahedralGaussian(const size_t N, const Domain& dom) :
    Gaussian() {
    domain_ = dom;
    construct(N);
    set_typeinfo();
}


OctahedralGaussian::OctahedralGaussian(const eckit::Parametrisation& params) :
    Gaussian() {
    // TODO: set domain from params
    domain_ = Domain::makeGlobal();

    size_t N;
    params.get("N",N);

    construct(N);
    set_typeinfo();
}


void OctahedralGaussian::construct(const size_t N) {
    std::vector<long> pl = computePL(N);
    setup_N_hemisphere(N,pl.data());
}


void OctahedralGaussian::set_typeinfo() {
    std::ostringstream s;
    s << "O"<< N();
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


eckit::ConcreteBuilderT1<Grid,OctahedralGaussian> builder_OctahedralGaussian (OctahedralGaussian::grid_type_str());


}  // namspace gaussian
}  // namespace grid
}  // namespace atlas
