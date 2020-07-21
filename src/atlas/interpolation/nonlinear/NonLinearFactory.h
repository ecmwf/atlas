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

#include <string>

#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


class NonLinear;


class NonLinearFactory : public util::Factory<NonLinearFactory> {
public:
    static std::string className() { return "NonLinearFactory"; }
    static const NonLinear* build( const std::string&, const util::Config& );
    using Factory::Factory;

private:
    virtual const NonLinear* make( const util::Config& ) = 0;
};


template <class T>
class NonLinearFactoryBuilder : public NonLinearFactory {
private:
    virtual const NonLinear* make( const util::Config& config ) override { return new T( config ); }

public:
    using NonLinearFactory::NonLinearFactory;
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
