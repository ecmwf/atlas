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

#include <vector>

#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"
#include "atlas/util/detail/Debug.h"

namespace atlas {
namespace util {

class PartitionPolygon;
//------------------------------------------------------------------------------------------------------

class ConvexSphericalPolygon {
public:
    static constexpr int MAX_GRIDCELL_EDGES = 4;
    static constexpr int MAX_SIZE           = 2 * MAX_GRIDCELL_EDGES + 1;

    ConvexSphericalPolygon();
    ConvexSphericalPolygon(const std::vector<PointLonLat>& points);
    //ConvexSphericalPolygon( const PartitionPolygon& );

    ConvexSphericalPolygon(const std::vector<PointXYZ>& points, const bool debug = false);

    operator bool() const { return valid_; }

    double area() const { return area_; }

    const PointXYZ& centroid() const { return centroid_; }

    static double norm_max(const PointXYZ& p, const PointXYZ& q);

    /*
   * @brief Point-on-segment test on great circle segments
   * @param[in] P given point in (x,y,z) coordinates
   * @return 
   */
    static bool between(const PointXYZ& p, const PointXYZ& p1, const PointXYZ& p2, const int debug = 0);

    /*
   * Point left of [p1,p2]
   * @param[in] P, p1, p2 given point in xyz-coordinates
   * @return 0:P_right_of_[p1,p2], -1:overlap_of_[P,p1]_and_[P,p2], 1:P_left_of_[p1,p2]
   */
    static int leftOf(const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2, const double tol, const int debug = 0);

    /*
   * @brief Segment-sph_polygon intersection
   * @param[in] s1, s2 segment endpoints in (x,y,z) coordinates
   * @param[in] start start with polygon segments [pol[start],pol[start+1]],...
   * @param[out] ip intersection point or nullptr
   * @return 0:no_intersection, 1:
   */
    int intersect(const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, const int start,
                  const bool debug = false) const;
    static PointXYZ common(const PointXYZ& s1, const PointXYZ& s2, const PointXYZ& p1, const PointXYZ& p2,
                           const int debug = 0);

    void clip(const PointXYZ& s1, const PointXYZ& s2, const int debug = 0);
    void clip(const PointLonLat& s1, const PointLonLat& s2, const int debug = 0);
    ConvexSphericalPolygon intersect(const ConvexSphericalPolygon& pol, const int debug = 0) const;

    /*
   * @brief check if two spherical polygons area equal
   * @param[in] P given point in (x,y,z) coordinates
   * @return true if equal vertices
   */
    bool equals(const ConvexSphericalPolygon& plg, const double deg_prec = 1e-10) const;

    /*
   * @return true:polygon is convex
   */
    bool validate();

    size_t size() const { return size_; }

    double cell_radius() const { return cell_radius_; }

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& out, const ConvexSphericalPolygon& p) {
        p.print(out);
        return out;
    }

    const PointXYZ& operator[](idx_t n) const {
        ATLAS_ASSERT(n < size_);
        return sph_coords_[n];
    }

    //private:

    void compute_area(const int debug = 0);

private:
    std::array<PointXYZ, MAX_SIZE> sph_coords_;
    PointXYZ centroid_;
    size_t size_;
    bool valid_;
    double area_;
    double cell_radius_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
