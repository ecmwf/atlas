/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_method_Method_h
#define atlas_interpolation_method_Method_h

#include <string>
#include <vector>
#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/memory/NonCopyable.h"
#include "eckit/memory/ScopedPtr.h"


namespace atlas {
  class Field;
  class FieldSet;
namespace mesh { class Mesh; }
}


namespace atlas {
namespace interpolation {
namespace method {


class Method : public eckit::NonCopyable {
public:

    typedef eckit::ScopedPtr< Method >  Ptr;
    typedef eckit::Configuration        Config;

    typedef eckit::linalg::SparseMatrix Matrix;
    typedef eckit::geometry::Point3     Point;
    typedef eckit::linalg::Triplet      Triplet;
    typedef std::vector< Triplet >      Triplets;

    enum { LON=0, LAT=1 };

    Method(const Config& config) : config_(config) {}
    virtual ~Method() {}

    /**
     * @brief Setup the interpolant matrix relating two (pre-partitioned) meshes
     * @param meshSource mesh containing source elements
     * @param meshTarget mesh containing target points
     */
    virtual void setup(mesh::Mesh& meshSource, mesh::Mesh& meshTarget) = 0;

    virtual void execute(const FieldSet& fieldsSource, FieldSet& fieldsTarget);
    virtual void execute(const Field&    fieldSource,  Field&    fieldTarget);

    const Matrix& matrix() const { return matrix_; }
    Matrix&       matrix()       { return matrix_; }

protected:

    static void normalise(Triplets& triplets);

    const Config& config_;

    Matrix matrix_;

};


struct MethodFactory {

    static Method *build(const std::string& name, const Method::Config&);

protected:

    std::string name_;
    virtual Method *make(const Method::Config&) = 0;

    MethodFactory(const std::string&);
    virtual ~MethodFactory();
};


template<class T>
struct MethodBuilder : public MethodFactory {

    MethodBuilder(const std::string &name) : MethodFactory(name) {}

private:
    virtual Method *make(const Method::Config& config) { return new T(config); }
};


}  // method
}  // interpolation
}  // atlas


#endif
