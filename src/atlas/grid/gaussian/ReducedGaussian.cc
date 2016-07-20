/*
 * (C) Copyright 1996-2016 ECMWF.
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


register_BuilderT1(Grid,ReducedGaussian,ReducedGaussian::grid_type_str());


std::string ReducedGaussian::grid_type_str() {
    return "reduced_gaussian";
}


std::string ReducedGaussian::className() {
    return "atlas.grid.gaussian.ReducedGaussian";
}


ReducedGaussian::ReducedGaussian(const size_t N, const long nlons[], const Domain& dom) :
    Gaussian() {
    domain_ = dom;
    setup_N_hemisphere(N,nlons);
    set_typeinfo();
}


ReducedGaussian::ReducedGaussian(const eckit::Parametrisation& params) :
    Gaussian() {
    // TODO: set domain from params
    domain_ = Domain::makeGlobal();

    setup(params);
    set_typeinfo();
}


void ReducedGaussian::setup( const eckit::Parametrisation& params ) {
    if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
    size_t N;
    params.get("N",N);

    std::vector<long> pl;
    params.get("pl",pl);
    setup_N_hemisphere(N,pl.data());
}


void ReducedGaussian::set_typeinfo() {
    std::stringstream s;
    s << "reduced_gaussian.N" << N();
    shortName_ = s.str();
    grid_type_ = grid_type_str();
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


} // namespace gaussian
} // namespace grid
} // namespace atlas
