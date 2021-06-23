/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/redistribution/Redistribution.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"

namespace atlas {

// Use redistribution implementation factory to make object.
Redistribution::Redistribution() : Handle() {};
Redistribution::Redistribution( const FunctionSpace& sourceFunctionSpace, const FunctionSpace& targetFunctionSpace ) :
    Handle( redistribution::detail::RedistributionImplFactory::build( sourceFunctionSpace, targetFunctionSpace ) ) {}

void Redistribution::execute( const Field& sourceField, Field& targetField ) const {
    get()->execute( sourceField, targetField );
    return;
}

void Redistribution::execute( const FieldSet& sourceFieldSet, FieldSet& targetFieldSet ) const {
    get()->execute( sourceFieldSet, targetFieldSet );
    return;
}

FunctionSpace& Redistribution::source() {
    return get()->source();
}

const FunctionSpace& Redistribution::source() const {
    return get()->source();
}

FunctionSpace& Redistribution::target() {
    return get()->target();
}

const FunctionSpace& Redistribution::target() const {
    return get()->target();
}

}  // namespace atlas
