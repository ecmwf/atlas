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


eckit::ConcreteBuilderT1<Grid, OctahedralGaussian> builder_OctahedralGaussian(OctahedralGaussian::grid_type_str());


std::string OctahedralGaussian::grid_type_str() {
    return "octahedral_gaussian";
}


std::string OctahedralGaussian::className() {
    return "atlas.grid.gaussian.OctahedralGaussian";
}


std::string OctahedralGaussian::gridType() const {
    return grid_type_str();
}


OctahedralGaussian::OctahedralGaussian(size_t N, const Domain& domain) :
    Gaussian() {
    std::vector<long> pl = computePL(N);
    domain_ = domain;

    setup_N_hemisphere(N, pl.data(), domain);
}


OctahedralGaussian::OctahedralGaussian(const eckit::Parametrisation& params) :
    Gaussian() {

    size_t N;
    if (!params.get("N", N)) {
        throw eckit::BadParameter("N missing in Params", Here());
    }

    domain_ = Domain::makeGlobal();
    std::vector<double> p_domain(4);
    if (params.get("domain", p_domain)) {
      domain_ = Domain(p_domain[0], p_domain[1], p_domain[2], p_domain[3]);
    }

    std::vector<long> pl = computePL(N);

    setup_N_hemisphere(N, pl.data(), domain_);
}


std::string OctahedralGaussian::shortName() const {
    if (shortName_.empty()) {
        std::stringstream s;
        s << "O" << N();
        if (!domain_.isGlobal()) {
            s << "-local";
        }
        shortName_ = s.str();
    }
    return shortName_;
}


std::vector<long> OctahedralGaussian::computePL(const size_t N) {
    const size_t start = 20;
    std::vector<long> pl(N);
    for(size_t jlat=0; jlat < N; ++jlat) {
        pl[jlat] = start + 4*jlat;
    }
    return pl;
}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas

