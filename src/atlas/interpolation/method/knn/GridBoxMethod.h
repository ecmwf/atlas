/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"

#include "atlas/functionspace.h"
#include "atlas/grid.h"


namespace atlas {
namespace interpolation {
namespace method {


class GridBoxMethod : public KNearestNeighboursBase {
public:
    GridBoxMethod( const Config& config ) : KNearestNeighboursBase( config ) {}
    virtual ~GridBoxMethod() override;

    virtual void print( std::ostream& ) const override;

protected:
    /**
     * @brief Create an interpolant sparse matrix relating two functionspaces, using grid-box average method
     * @param source functionspace containing source points
     * @param target functionspace containing target points
     */
    virtual void setup( const FunctionSpace& source, const FunctionSpace& target ) override;

    virtual void setup( const Grid& source, const Grid& target ) override;

    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

private:
    FunctionSpace source_;
    FunctionSpace target_;
    Grid sourceGrid_;
    Grid targetGrid_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
