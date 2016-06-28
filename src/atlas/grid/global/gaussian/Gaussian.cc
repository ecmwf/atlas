/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/global/gaussian/Gaussian.h"

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/grid/global/gaussian/latitudes/Latitudes.h"


namespace atlas {
namespace grid {
namespace global {
namespace gaussian {


void Gaussian::LatitudesNorthPoleToEquator(const size_t N, double latitudes[]) {
    global::gaussian::latitudes::gaussian_latitudes_npole_equator(N,latitudes);
}


void Gaussian::LatitudesNorthPoleToSouthPole(const size_t N, double latitudes[]) {
    global::gaussian::latitudes::gaussian_latitudes_npole_spole(N,latitudes);
}


void Gaussian::QuadratureNorthPoleToEquator(const size_t N, double weights[]) {
    std::vector<double> lats(N);
    latitudes::gaussian_quadrature_npole_equator(N,lats.data(),weights);
}


void Gaussian::QuadratureNorthPoleToSouthPole(const size_t N, double weights[]) {
    std::vector<double> lats(2*N);
    latitudes::gaussian_quadrature_npole_equator(N,lats.data(),weights);
}


std::string Gaussian::grid_type_str() {
    return "gaussian";
}


std::string Gaussian::className() {
    return "atlas.grid.global.gaussian.Gaussian";
}


Gaussian::Gaussian(const Domain& dom) :
    Structured(dom) {
}


void Gaussian::setup_N_hemisphere(const size_t N, const long pl[]) {
    Structured::N_ = N;
    // hemisphere
    std::vector<double> lats (N);
    LatitudesNorthPoleToEquator(N,lats.data());
    Structured::setup_lat_hemisphere(N,lats.data(),pl);
}


eckit::Properties Gaussian::spec() const {
    eckit::Properties grid_spec;

    grid_spec.set("grid_type", gridType() );
    grid_spec.set("short_name",shortName());

    grid_spec.set("nlat",nlat());
    grid_spec.set("N", N() );

    grid_spec.set("pl",eckit::makeVectorValue(pl()));

    return grid_spec;
}


}  // namespace gaussian
}  // namespace global
}  // namespace grid
}  // namespace atlas
