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
#include "atlas/redistribution/detail/RedistributeGeneric.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"


namespace atlas {

// Use redistribution implementation factory to make object.
Redistribution::Redistribution(): Handle(){};
Redistribution::Redistribution(const FunctionSpace& sourceFunctionSpace, const FunctionSpace& targetFunctionSpace,
                               const util::Config& config):
    Handle([&]() -> redistribution::detail::RedistributionImpl* {
        // If no heterogeneous_redistribution config is set, or if it is not set to "true", assert that
        // the source and target function spaces are the same type (default)
        std::string isHeteroRedistString = "false";
        config.get("heterogeneous_redistribution", isHeteroRedistString);
        bool isHeteroRedist = (isHeteroRedistString == "true");

        if (!isHeteroRedist) {
            ATLAS_ASSERT(sourceFunctionSpace.type() == targetFunctionSpace.type());
        };

        std::string type = redistribution::detail::RedistributeGeneric::static_type();
        config.get("type", type);

        auto impl = redistribution::detail::RedistributionImplFactory::build(type);
        impl->setup(sourceFunctionSpace, targetFunctionSpace);
        return impl;
    }()) {}

void Redistribution::execute(const Field& source, Field& target) const {
    get()->execute(source, target);
    return;
}

void Redistribution::execute(const FieldSet& source, FieldSet& target) const {
    get()->execute(source, target);
    return;
}

const FunctionSpace& Redistribution::source() const {
    return get()->source();
}

const FunctionSpace& Redistribution::target() const {
    return get()->target();
}

}  // namespace atlas
