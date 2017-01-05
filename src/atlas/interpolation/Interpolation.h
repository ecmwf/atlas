/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_Interpolation_h
#define atlas_interpolation_Interpolation_h

#include <string>
#include <vector>
#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/memory/NonCopyable.h"
#include "eckit/memory/ScopedPtr.h"


namespace atlas {
namespace mesh { class Mesh; }
}


namespace atlas {
namespace interpolation {


class Interpolation : public eckit::NonCopyable {
public:

    typedef eckit::ScopedPtr<Interpolation> Ptr;
    typedef eckit::Configuration            Config;

    typedef eckit::linalg::SparseMatrix Matrix;
    typedef eckit::geometry::Point3     Point;
    typedef eckit::linalg::Triplet      Triplet;
    typedef std::vector< Triplet >      Triplets;

    enum { LON=0, LAT=1 };

    Interpolation(const Config& config) : config_(config) {}
    virtual ~Interpolation() {}

    /**
     * @brief Create an interpolant sparse matrix relating two (pre-partitioned) meshes
     * @param meshSource mesh containing source elements
     * @param meshTarget mesh containing target points
     */
    virtual void execute(Matrix& matrix, mesh::Mesh& meshSource, mesh::Mesh& meshTarget) const = 0;

protected:

    static void normalise(Triplets& triplets);

    const Config& config_;

};


struct InterpolationFactory {

    static Interpolation *build(const std::string& name, const Interpolation::Config&);

protected:

    std::string name_;
    virtual Interpolation *make(const Interpolation::Config&) = 0;

    InterpolationFactory(const std::string&);
    virtual ~InterpolationFactory();
};


template<class T>
struct InterpolationBuilder : public InterpolationFactory {

    InterpolationBuilder(const std::string &name) : InterpolationFactory(name) {}

private:
    virtual Interpolation *make(const Interpolation::Config& config) { return new T(config); }
};


}  // interpolation
}  // atlas


#endif
