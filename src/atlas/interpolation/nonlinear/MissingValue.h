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

#include "eckit/config/Parametrisation.h"

#include "atlas/util/Factory.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


/// @brief Missing values indicator base class
struct MissingValue {
    using Config = eckit::Parametrisation;
    virtual ~MissingValue();
    virtual bool operator()( const double& ) const = 0;
};


/// @brief Missing values indicator factory
struct MissingValueFactory : util::Factory<MissingValueFactory> {
    using Config = MissingValue::Config;
    using Factory::Factory;
    static std::string className() { return "MissingValueFactory"; }
    static const MissingValue* build( const std::string&, const Config& );

private:
    virtual const MissingValue* make( const Config& ) = 0;
};


/// @brief Missing values indicator builder for factory registration
template <class T>
class MissingValueFactoryBuilder : public MissingValueFactory {
private:
    virtual const MissingValue* make( const Config& config ) override { return new T( config ); }

public:
    using MissingValueFactory::MissingValueFactory;
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
