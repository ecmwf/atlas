/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/gaussian/ClassicGaussian.h"

#include "eckit/memory/Builder.h"
#include "atlas/grid/gaussian/classic/PointsPerLatitude.h"


namespace atlas {
namespace grid {
namespace gaussian {


eckit::ConcreteBuilderT1<Grid, ClassicGaussian> builder_ClassicGaussian(ClassicGaussian::grid_type_str());


std::string ClassicGaussian::grid_type_str() {
    return "classic_gaussian";
}


std::string ClassicGaussian::className() {
    return "atlas.grid.gaussian.ClassicGaussian";
}


std::string ClassicGaussian::gridType() const {
    return grid_type_str();
}


ClassicGaussian::ClassicGaussian(const size_t N, const Domain& domain) :
    Gaussian() {
    std::vector<long> pl(N);
    classic::points_per_latitude_npole_equator(N, pl.data());
    domain_ = domain;

    setup_N_hemisphere(N, pl.data(), domain);
}


ClassicGaussian::ClassicGaussian(const eckit::Parametrisation& params) :
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

    std::vector<long> pl(N);
    classic::points_per_latitude_npole_equator(N, pl.data());

    setup_N_hemisphere(N, pl.data(), domain_);
}


std::string ClassicGaussian::shortName() const {
    if (shortName_.empty()) {
        std::stringstream s;
        s << "N" << N();
        if (!domain_.isGlobal()) {
            s << "-local";
        }
        shortName_ = s.str();
    }
    return shortName_;
}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas

