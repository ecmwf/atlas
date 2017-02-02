/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/gaussian/ReducedGaussian.h"

#include <typeinfo>
#include "eckit/memory/Builder.h"


namespace atlas {
namespace grid {
namespace gaussian {


eckit::ConcreteBuilderT1<Grid, ReducedGaussian> builder_ReducedGaussian(ReducedGaussian::grid_type_str());


std::string ReducedGaussian::grid_type_str() {
    return "reduced_gaussian";
}


std::string ReducedGaussian::className() {
    return "atlas.grid.gaussian.ReducedGaussian";
}


std::string ReducedGaussian::gridType() const {
    return grid_type_str();
}


ReducedGaussian::ReducedGaussian(const size_t N, const long nlons[], const Domain& domain) :
    Gaussian() {
    domain_ = domain;

    setup_N_hemisphere(N, nlons, domain_);
}


ReducedGaussian::ReducedGaussian(const eckit::Parametrisation& params) :
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

    std::vector<long> pl;
    if (!params.get("pl", pl)) {
        throw eckit::BadParameter("pl missing in Params", Here());
    }

    setup_N_hemisphere(N, pl.data(), domain_);
}


std::string ReducedGaussian::shortName() const {
    if (shortName_.empty()) {
        std::stringstream s;
        s << "reduced_gaussian.N" << N();
        if (!domain_.isGlobal()) {
            s << "-local";
        }
        shortName_ = s.str();
    }
    return shortName_;
}


extern "C" {


    Structured* atlas__grid__gaussian__ReducedGaussian_int(size_t N, int pl[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl,pl+N);
        return new ReducedGaussian(N,pl_vector.data());
    }


    Structured* atlas__grid__gaussian__ReducedGaussian_long(size_t N, long pl[]) {
        return new ReducedGaussian(N,pl);
    }


}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas

