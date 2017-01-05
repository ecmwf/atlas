/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_NearestNeighbour_h
#define atlas_interpolation_NearestNeighbour_h

#include "Interpolation.h"

#include <string>
#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"
#include "atlas/array/ArrayView.h"
#include "atlas/interpolation/PointIndex3.h"


namespace atlas {
namespace interpolation {


class NearestNeighbour : public Interpolation {
public:

    NearestNeighbour(const Config& config) : Interpolation(config) {}
    virtual ~NearestNeighbour() {}

    /**
     * @brief Create an interpolant sparse matrix relating two (pre-partitioned) meshes,
     * using nearest neighbour method
     * @param meshSource mesh containing source elements
     * @param meshTarget mesh containing target points
     */
    void execute(Matrix& matrix, mesh::Mesh& meshSource, mesh::Mesh& meshTarget) const;

};


}  // interpolation
}  // atlas


#endif
