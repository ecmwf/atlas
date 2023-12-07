/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

//     Purpose.
//     --------
//           eq_regions provides the code to perform a high level
//           partitioning of the surface of a sphere into regions of
//           equal area and small diameter.
//
//     Background.
//     -----------
//     This C++ implementation is ported from the MATLAB
//     "Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox" of the
//     same name developed by Paul Leopardi at the University of New South
//     Wales.
//     This version has been coded specifically for the case of partitioning the
//     surface of a sphere or S^dim (where dim=2) as denoted in the original
//     code.
//     Only a subset of the original eq_regions package has been coded to
//     determin
//     the high level distribution of regions on a sphere, as the detailed
//     distribution of grid points to each region is left to implentatios.
//
//     The following copyright notice for the eq_regions package is included
//     from
//     the original MatLab release.
//
//     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     + Release 1.10 2005-06-26 +
//     + +
//     + Copyright (c) 2004, 2005, University of New South Wales +
//     + +
//     + Permission is hereby granted, free of charge, to any person obtaining +
//     + a copy of this software and associated documentation files (the +
//     + "Software"), to deal in the Software without restriction, including +
//     + without limitation the rights to use, copy, modify, merge, publish, +
//     + distribute, sublicense, and/or sell copies of the Software, and to +
//     + permit persons to whom the Software is furnished to do so, subject to +
//     + the following conditions: +
//     + +
//     + The above copyright notice and this permission notice shall be included
//     +
//     + in all copies or substantial portions of the Software. +
//     + +
//     + THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, +
//     + EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF +
//     + MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//     +
//     + IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY +
//     + CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, +
//     + TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE +
//     + SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. +
//     + +
//     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#pragma once

#include <vector>

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

void eq_caps(int N, std::vector<int>& n_regions, std::vector<double>& s_cap);
void eq_regions(int N, double xmin[], double xmax[], double ymin[], double ymax[]);

class EqualAreaPartitioner;
class EqualRegionsPartitioner : public Partitioner {
public:
    EqualRegionsPartitioner();

    EqualRegionsPartitioner(int N);
    EqualRegionsPartitioner(int N, const eckit::Parametrisation& config);
    EqualRegionsPartitioner(const eckit::Parametrisation& config);

    void where(int partition, int& band, int& sector) const;
    int nb_bands() const { return bands_.size(); }
    int nb_regions(int band) const { return sectors_[band]; }

    using Partitioner::partition;
    virtual void partition(const Grid&, int part[]) const;

    virtual std::string type() const { return "equal_regions"; }

public:
    // Node struct that holds the longitude and latitude in millidegrees
    // (integers)
    // This structure is used in sorting algorithms, and uses less memory than
    // if x and y were in double precision.
    struct NodeInt {
        int x, y;
        int n;
        bool operator!=(const NodeInt& other) const { return n != other.n; }
        bool operator==(const NodeInt& other) const { return n == other.n; }
        void swap(NodeInt& other) {
            auto _swap = [](int& a, int& b) {
                int tmp = a;
                a       = b;
                b       = tmp;
            };
            _swap(x, other.x);
            _swap(y, other.y);
            _swap(n, other.n);
        }
        friend void swap(NodeInt& a, NodeInt& b) { a.swap(b); }
    };

private:
    void init();
    // Doesn't matter if nodes[] is in degrees or radians, as a sorting
    // algorithm is used internally
    void partition(int nb_nodes, NodeInt nodes[], int part[]) const;

    friend class EqualAreaPartitioner;
    // y in radians
    int band(const double& y) const;

    // x in radians
    int sector(int band, const double& x) const;

    // x and y in radians
    int partition(const double& x, const double& y) const;

private:
    int N_;
    std::vector<double> bands_;
    std::vector<int> sectors_;
    enum class Coordinates
    {
        XY,
        LONLAT,
    };
    Coordinates coordinates_ = Coordinates::XY;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
