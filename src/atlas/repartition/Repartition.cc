/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/repartition/Repartition.h"
#include "atlas/repartition/detail/RepartitionImplFactory.h"

namespace atlas {

  // Use repartition implementation factory to make object.
  Repartition::Repartition (
    const FunctionSpace& sourceFunctionSpace,
    const FunctionSpace& targetFunctionSpace) :
    Handle (repartition::RepartitionImplFactory::build (
      sourceFunctionSpace, targetFunctionSpace)) {}

  void Repartition::execute(
    const Field& sourceField, Field& targetField) const {

    get()->execute(sourceField, targetField);
    return;
  }

  void Repartition::execute(
    const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const {

    get()->execute(sourceFieldSet, targetFieldSet);
    return;
  }

  FunctionSpace& Repartition::getSourceFunctionSpace () {

    return get()->getSourceFunctionSpace();
  }

  const FunctionSpace& Repartition::getSourceFunctionSpace () const {

    return get()->getSourceFunctionSpace();
  }

  FunctionSpace& Repartition::getTargetFunctionSpace () {

    return get()->getTargetFunctionSpace();
  }

  const FunctionSpace& Repartition::getTargetFunctionSpace () const {

    return get()->getTargetFunctionSpace();
  }

}
