/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grids_ReducedGaussianGrid_h
#define atlas_grids_ReducedGaussianGrid_h

#include "atlas/grid/global/Structured.h"


namespace atlas {
namespace grid {
namespace global {
namespace gaussian {


/**
 * @brief Abstract (Reduced) Gaussian Grid
 *
 * This grid is a special case of the class Structured, in which
 * the latitudes are distributed according to the roots of the
 * Legendre Polynomials, and a equidistant distribution in zonal
 * direction, which reduce in number going closer towards poles,
 * essentially making the grid more uniform on the sphere
 * It can be constructed with following definition:
 *   N   = number of latitudes in hemisphere
 *   npts_per_lat[] = number of points on each latitude
 */
class Gaussian: public Structured {

  public:

    /**
     * @brief Compute gaussian latitudes between North pole and equator
     * @param N         [in]  Number of latitudes between pole and equator (Gaussian N number)
     * @param latitudes [out] latitudes in degrees
     */
    static void LatitudesNorthPoleToEquator(const size_t N, double latitudes[]);

    /**
     * @brief Compute gaussian latitudes between North pole and South pole
     * @param N         [in]  Number of latitudes between pole and equator (Gaussian N number)
     * @param latitudes [out] latitudes in degrees  (size 2*N)
     */
    static void LatitudesNorthPoleToSouthPole(const size_t N, double latitudes[]);

    /**
     * @brief Compute gaussian quadrature weights between North pole and equator
     * @param N         [in]  Number of latitudes between pole and equator (Gaussian N number)
     * @param weights   [out] quadrature weights
     * @note Weights are normalized so that their sum equals 0.5 between pole and equator, or 1. from pole to pole
     */
    static void QuadratureNorthPoleToEquator (const size_t N, double weights[]);

    /**
     * @brief Compute gaussian quadrature weights between North pole and South pole
     * @param N         [in]  Number of latitudes between pole and equator (Gaussian N number)
     * @param weights   [out] quadrature weights   (size 2*N)
     * @note Weights are normalized so that their sum equals 0.5 between pole and equator, or 1. from pole to pole
     */
    static void QuadratureNorthPoleToSouthPole (const size_t N, double weights[]);

    static std::string grid_type_str();

    static std::string className();

    virtual eckit::Properties spec() const;

  protected:

    /// to be used only by derived types
    Gaussian(const Domain& dom=Domain::makeGlobal());

    void setup_N_hemisphere(const size_t N, const long pl[]);

    virtual void set_typeinfo() = 0;

};


}  // namespace gaussian
}  // namespace global
}  // namespace grid
}  // namespace atlas


#endif
