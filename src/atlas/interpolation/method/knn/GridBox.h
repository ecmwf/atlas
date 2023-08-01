/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <iosfwd>
#include <vector>


namespace atlas {
class Grid;
}


namespace atlas {
namespace interpolation {
namespace method {


class GridBox {
public:
    // -- Exceptions
    // None

    // -- Constructors

    GridBox(double north, double west, double south, double east);

    // -- Destructor
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods

    double area() const;
    double diagonal() const;
    bool intersects(GridBox&) const;

    double north() const { return north_; }
    double west() const { return west_; }
    double south() const { return south_; }
    double east() const { return east_; }
    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

protected:
    // -- Members
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None

private:
    // -- Members

    double north_;
    double west_;
    double south_;
    double east_;

    // -- Methods

    void print(std::ostream&) const;

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

    friend std::ostream& operator<<(std::ostream& s, const GridBox& p) {
        p.print(s);
        return s;
    }
};


struct GridBoxes : std::vector<GridBox> {
    GridBoxes(const Grid&, bool gaussianWeightedLatitudes = true);
    GridBoxes();
    double getLongestGridBoxDiagonal() const;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
