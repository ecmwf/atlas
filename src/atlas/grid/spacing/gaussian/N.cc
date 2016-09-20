/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Jan 2015


#include "atlas/grid/spacing/gaussian/N.h"


namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {


std::string GaussianLatitudes::className() {
    return "GaussianLatitudes";
}


void GaussianLatitudes::assign(double lats[], const size_t size) const {
    ASSERT( size >= lats_.size() );
    for(size_t jlat=0; jlat < lats_.size(); ++jlat)
        lats[jlat] = lats_[jlat];
}


void GaussianLatitudes::assign(std::vector<double>& lats) const {
    lats = lats_;
}


template<typename CONCRETE>
void load() {
    eckit::ConcreteBuilderT0<GaussianLatitudes,CONCRETE> builder("tmp");
}


void regist() {
    load<N16>();
    load<N24>();
}


}  // namespace gaussian
}  // namespace spacing
}  // namespace grid
}  // namespace atlas

