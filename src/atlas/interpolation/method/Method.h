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

#include <iosfwd>
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
}  // namespace atlas

namespace atlas {
namespace interpolation {

class Method : public eckit::Owned {
public:
    typedef eckit::Parametrisation Config;

    Method( const Config& config ) : config_( config ) {}
    virtual ~Method() {}

    /**
   * @brief Setup the interpolator relating two functionspaces
   * @param source functionspace containing source elements
   * @param target functionspace containing target points
   */
    virtual void setup( const FunctionSpace& source, const FunctionSpace& target ) = 0;

    virtual void execute( const FieldSet& source, FieldSet& target ) const;
    virtual void execute( const Field& source, Field& target ) const;

    virtual void print( std::ostream& ) const = 0;

protected:
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;
    using Matrix   = eckit::linalg::SparseMatrix;

    static void normalise( Triplets& triplets );

    const Config& config_;

    // NOTE : Matrix-free or non-linear interpolation operators do not have
    // matrices,
    //        so do not expose here, even though only linear operators are now
    //        implemented.
    Matrix matrix_;
};

struct MethodFactory {
    static Method* build( const std::string& name, const Method::Config& );

protected:
    std::string name_;
    virtual Method* make( const Method::Config& ) = 0;

    MethodFactory( const std::string& );
    virtual ~MethodFactory();
};

template <class T>
struct MethodBuilder : public MethodFactory {
    MethodBuilder( const std::string& name ) : MethodFactory( name ) {}

private:
    virtual Method* make( const Method::Config& config ) { return new T( config ); }
};

}  // namespace interpolation
}  // namespace atlas
