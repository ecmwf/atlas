/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <string>
#include <vector>
#include "eckit/config/Configuration.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"


namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
  namespace functionspace {
    class FunctionSpaceImpl;
  }
  namespace field {
    class FieldImpl;
    class FieldSetImpl;
  }
}


namespace atlas {
namespace interpolation {
namespace method {


class Method : public eckit::Owned {
public:

    typedef eckit::Parametrisation Config;

    Method(const Config& config) : config_(config) {}
    virtual ~Method() {}

    /**
     * @brief Setup the interpolator relating two functionspaces
     * @param source functionspace containing source elements
     * @param target functionspace containing target points
     */
    virtual void setup(const FunctionSpace& source, const FunctionSpace& target) = 0;

    virtual void execute(const FieldSet& source, FieldSet& target) const;
    virtual void execute(const Field&    source, Field&    target) const;

protected:

    typedef eckit::linalg::Triplet      Triplet;
    typedef std::vector< Triplet >      Triplets;
    typedef eckit::linalg::SparseMatrix Matrix;

    static void normalise(Triplets& triplets);

    const Config& config_;

    // NOTE : Matrix-free or non-linear interpolation operators do not have matrices,
    //        so do not expose here, even though only linear operators are now implemented.
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

class InterpolationMethod {
public:

  using Implementation = method::Method;
  using Config = Implementation::Config;

  InterpolationMethod() {}
  InterpolationMethod( const InterpolationMethod& );
  InterpolationMethod( const Config&, const FunctionSpace& source, const FunctionSpace& target );

  void execute(const FieldSet& source, FieldSet& target) const { get()->execute(source,target); }
  void execute(const Field&    source, Field&    target) const { get()->execute(source,target); }

  const Implementation* get() const { return implementation_.get(); }

  operator bool() const { return implementation_; }

private:

  eckit::SharedPtr<const Implementation> implementation_;

};


extern "C" {

InterpolationMethod::Implementation* atlas__Interpolation__new(const eckit::Parametrisation* config, const functionspace::FunctionSpaceImpl* source, const functionspace::FunctionSpaceImpl* target);
void atlas__Interpolation__delete(InterpolationMethod::Implementation* This);
void atlas__Interpolation__execute_field(InterpolationMethod::Implementation* This, const field::FieldImpl* source, field::FieldImpl* target);
void atlas__Interpolation__execute_fieldset(InterpolationMethod::Implementation* This, const field::FieldSetImpl* source, field::FieldSetImpl* target);

}

}  // interpolation
}  // atlas
