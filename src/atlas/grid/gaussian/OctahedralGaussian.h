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
/// @date Jan 2014


#ifndef atlas_grid_gaussian_OctahedralGaussian_h
#define atlas_grid_gaussian_OctahedralGaussian_h

#include "atlas/grid/gaussian/Gaussian.h"


namespace atlas {
namespace grid {
namespace gaussian {


class OctahedralGaussian : public Gaussian {

  public:

    static std::string grid_type_str();

    static std::string className();

    std::string gridType() const;

    OctahedralGaussian(const size_t N, const Domain& domain=Domain::makeGlobal());

    OctahedralGaussian(const eckit::Parametrisation&);

    virtual const Domain& domain() const {
        return domain_;
    }

  private:

    std::string shortName() const;

    /// Computes the PL for the Octohedral distribution
    /// number of points at latitude closest to pole
    static std::vector<long> computePL(const size_t N);

  private:

    mutable std::string shortName_;

    Domain domain_;

};


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas


#endif
