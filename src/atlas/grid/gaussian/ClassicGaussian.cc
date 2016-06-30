/*
 * (C) Copyright 1996-2016 ECMWF.
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


register_BuilderT1(Grid,ClassicGaussian,ClassicGaussian::grid_type_str());


std::string ClassicGaussian::grid_type_str() {
    return "classic_gaussian";
}


std::string ClassicGaussian::className() {
    return "atlas.grid..gaussian.ClassicGaussian";
}


void ClassicGaussian::set_typeinfo() {
    std::stringstream s;
    s << "N" << N();
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


ClassicGaussian::ClassicGaussian(const size_t N) :
    Gaussian() {
    std::vector<long> nlon(N);
    classic::points_per_latitude_npole_equator(N,nlon.data());
    setup_N_hemisphere(N,nlon.data());
    set_typeinfo();
}


ClassicGaussian::ClassicGaussian(const eckit::Parametrisation& params) :
    Gaussian() {
    size_t N;
    if( ! params.get("N",N) )
        throw eckit::BadParameter("N missing in Params",Here());
    std::vector<long> nlon(N);
    classic::points_per_latitude_npole_equator(N,nlon.data());
    setup_N_hemisphere(N,nlon.data());
    set_typeinfo();
}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas
