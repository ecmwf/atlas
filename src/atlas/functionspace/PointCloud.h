/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/field/Field.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace functionspace {

//------------------------------------------------------------------------------------------------------

namespace detail { 

class PointCloud : public FunctionSpaceImpl
{
public:
    PointCloud(const std::vector<PointXY>&);
    PointCloud(const Field& lonlat);
    PointCloud(const Field& lonlat, const Field& ghost);
    virtual ~PointCloud() {}
    virtual std::string name() const { return "PointCloud"; }
    virtual operator bool() const { return true; }
    virtual size_t footprint() const { return sizeof(*this); }

    const Field& lonlat() const { return lonlat_; }
    const Field& ghost() const;
    const size_t size() const { return lonlat_.shape(0); }

  private:

    Field lonlat_;
    mutable Field ghost_;
};

} // namespace datail

//------------------------------------------------------------------------------------------------------

class PointCloud : public FunctionSpace {
public:

  PointCloud( const FunctionSpace& );
  PointCloud( const Field& points );
  PointCloud( const std::vector<PointXY>& );

  operator bool() const { return valid(); }
  bool valid() const { return functionspace_; }

  const Field& lonlat() const { return functionspace_->lonlat(); }
  const Field& ghost() const { return functionspace_->ghost(); }
  const size_t size() const { return functionspace_->size(); }

private:

  const detail::PointCloud* functionspace_;
};

} // namespace functionspace
} // namespace atlas
