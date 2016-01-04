/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "atlas/runtime/Log.h"
#include "eckit/log/BigNum.h"
#include "eckit/value/Value.h"
#include "eckit/geometry/RotateGrid.h"

#include "atlas/grids/LocalGrid.h"

using namespace eckit;
using namespace eckit::geometry;
using namespace std;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

LocalGrid::LocalGrid(Grid* grid, const Domain& domain)
    : Grid(domain),
      grid_(grid),
      localPts_(0) {

    ASSERT(grid);
    ASSERT( !dynamic_cast<LocalGrid*>(grid) ); // FIXME: add support of LocalGrid of LocalGrid

    cropPoints();
}

LocalGrid::~LocalGrid() {}

std::string LocalGrid::shortName() const {

    if ( shortName_.empty() ) {
        shortName_ = "local." + grid_->shortName();
    }
    return shortName_;
}

void LocalGrid::hash(eckit::MD5& md5) const {

    md5.add("local.");

    grid_->hash(md5);

    domain().hash(md5);
}

void LocalGrid::cropPoints() {

    std::vector<Grid::Point> gpts;
    grid_->lonlat(gpts);

//    Log::info() << " DOMAIN: " << domain_ << std::endl;
//    Log::info() << " ----> NESTED GRID: " << gpts.size() << std::endl;

    size_t accepted = 0;
    size_t discarded = 0;


    for (size_t i = 0; i < gpts.size(); ++i) {
        const Grid::Point& p = gpts[i];
        if ( domain_.contains(p) ) {
//            Log::info() << "  ++ POINT " << p << std::endl;
            localPts_.push_back(p);
            ++accepted;
        }
        else {
//            Log::info() << " DISCARDED POINT " << p << std::endl;
            ++discarded;
        }
    }

    Log::info() << "Local Grid contains " << eckit::BigNum(accepted) << std::endl;
    Log::info() << "Local Grid discards " << eckit::BigNum(discarded) << std::endl;
}

size_t LocalGrid::npts() const { return localPts_.size(); }

void LocalGrid::print(ostream& os) const {
    os << "LocalGrid("
       << "Domain:" << domain_
       << ",Grid:" << *grid_
       << ",npts:" << npts()
       << ")";
}

void LocalGrid::lonlat(std::vector<Point>& pts) const { pts = localPts_; }

std::string LocalGrid::gridType() const {
    std::ostringstream os;
    os << "local." << grid_->gridType();
    return os.str();
}

eckit::Properties LocalGrid::spec() const {
    NOTIMP;
}


//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
