/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"

#include <forward_list>

#include "atlas/functionspace.h"
#include "atlas/interpolation/method/knn/GridBox.h"


namespace atlas {
namespace interpolation {
namespace method {


class GridBoxMethod : public KNearestNeighboursBase {
public:
    GridBoxMethod(const Config&);
    virtual ~GridBoxMethod() override;

protected:
    virtual void print(std::ostream&) const override;

    using KNearestNeighboursBase::do_setup;
    /**
     * @brief Create an interpolant sparse matrix relating two functionspaces, using grid-box average method
     * @param source functionspace containing source points
     * @param target functionspace containing target points
     */
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&) override;
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

    bool intersect(size_t i, const GridBox& iBox, const util::IndexKDTree::ValueList&, std::vector<Triplet>&) const;

    virtual void do_execute(const FieldSet& source, FieldSet& target, Metadata&) const override = 0;
    virtual void do_execute(const Field& source, Field& target, Metadata&) const override       = 0;

    virtual Cache createCache() const override;

protected:
    static void giveUp(const std::forward_list<size_t>&);

    FunctionSpace source_;
    FunctionSpace target_;

    GridBoxes sourceBoxes_;
    GridBoxes targetBoxes_;

    double searchRadius_;

    mutable std::forward_list<size_t> failures_;

    bool matrixFree_;
    bool failEarly_;
    bool gaussianWeightedLatitudes_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
