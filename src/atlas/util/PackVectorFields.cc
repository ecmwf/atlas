/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/util/PackVectorFields.h"

#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "atlas/array.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/functionspace.h"
#include "atlas/option.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "eckit/config/LocalConfiguration.h"

namespace atlas {
namespace util {

namespace {

using array::helpers::arrayForEachDim;
using eckit::LocalConfiguration;

void addOrReplaceField(FieldSet& fieldSet, const Field& field) {
    const auto fieldName = field.name();
    if (fieldSet.has(fieldName)) {
        fieldSet[fieldName] = field;
    }
    else {
        fieldSet.add(field);
    }
}

Field& getOrCreateField(FieldSet& fieldSet, const FunctionSpace& functionSpace, const Config& config) {
    const auto fieldName = config.getString("name");
    if (!fieldSet.has(fieldName)) {
        auto field = functionSpace.createField(config);
        field->set_dirty(false);  // We will inherit the dirty state from the component fields.
        fieldSet.add(field);
    }
    return fieldSet[fieldName];
}

void checkFieldCompatibility(const Field& componentField, const Field& vectorField) {
    ATLAS_ASSERT(componentField.functionspace().size() == vectorField.functionspace().size());
    ATLAS_ASSERT(componentField.levels() == vectorField.levels());
    ATLAS_ASSERT(componentField.variables() == 0);
    ATLAS_ASSERT(vectorField.variables() > 0);

    const auto checkStandardShape = [](const Field& field) {
        // Check for "standard" Atlas field shape.
        auto dim         = 0;
        const auto rank  = field.rank();
        const auto shape = field.shape();
        if (field.functionspace().size() != shape[dim++]) {
            return false;
        }
        if (const auto levels = field.levels(); levels && (dim >= rank || levels != shape[dim++])) {
            return false;
        }
        if (const auto variables = field.variables(); variables && (dim >= rank || variables != shape[dim++])) {
            return false;
        }
        if (dim != rank) {
            return false;
        }
        return true;
    };

    ATLAS_ASSERT(checkStandardShape(componentField));
    ATLAS_ASSERT(checkStandardShape(vectorField));
}

template <typename ComponentField, typename VectorField, typename Functor>
void copyFieldData(ComponentField& componentField, VectorField& vectorField, const Functor& copier) {
    checkFieldCompatibility(componentField, vectorField);

    auto componentViewVariant = array::make_view_variant(componentField);

    const auto componentVisitor = [&](auto componentView) {
        if constexpr (array::is_rank<1, 2>(componentView)) {
            using ComponentView          = std::decay_t<decltype(componentView)>;
            constexpr auto ComponentRank = ComponentView::rank();
            using Value                  = typename ComponentView::non_const_value_type;

            auto vectorView     = array::make_view<Value, ComponentRank + 1>(vectorField);
            constexpr auto Dims = std::make_integer_sequence<int, ComponentRank>{};
            arrayForEachDim(Dims, execution::par, std::tie(componentView, vectorView), copier);
        }
        else {
            ATLAS_THROW_EXCEPTION("Unsupported component field rank: " + std::to_string(componentView.rank()));
        }
    };

    std::visit(componentVisitor, componentViewVariant);
}

}  // namespace

FieldSet pack_vector_fields(const FieldSet& fields, FieldSet packedFields) {
    // Get the number of variables for each vector field.
    auto vectorSizeMap = std::map<std::string, idx_t>{};
    for (const auto& field : fields) {
        auto vectorFieldName = std::string{};
        if (field.metadata().get("vector_component.vector_field_name", vectorFieldName)) {
            ++vectorSizeMap[vectorFieldName];
        }
    }

    // Pack vector fields.
    for (const auto& field : fields) {
        auto vectorFieldName = std::string{};
        if (!field.metadata().get("vector_component.vector_field_name", vectorFieldName)) {
            // Not a vector component field.
            addOrReplaceField(packedFields, field);
            continue;
        }

        // Field is vector field component.
        auto vectorSize  = vectorSizeMap[vectorFieldName];
        auto vectorIndex = size_t{};
        field.metadata().get("vector_component.index", vectorIndex);
        const auto& componentField = field;

        // Get or create vector field.
        const auto vectorFieldConfig = option::name(vectorFieldName) | option::levels(componentField.levels()) |
                                       option::vector(vectorSize) | option::datatype(componentField.datatype());
        auto& vectorField = getOrCreateField(packedFields, componentField.functionspace(), vectorFieldConfig);

        // Copy field data.
        const auto copier = [&](auto&& componentElem, auto&& vectorElem) { vectorElem(vectorIndex) = componentElem; };
        copyFieldData(componentField, vectorField, copier);

        // Copy metadata.
        const auto componentFieldMetadata = componentField.metadata();
        auto componentFieldMetadataVector = std::vector<LocalConfiguration>(vectorSize);
        vectorField.metadata().get("component_field_metadata", componentFieldMetadataVector);
        componentFieldMetadataVector[vectorIndex] = componentFieldMetadata;
        vectorField.metadata().set("component_field_metadata", componentFieldMetadataVector);

        // If any component is dirty, the whole field is dirty.
        vectorField.set_dirty(vectorField.dirty() || componentField.dirty());
    }
    return packedFields;
}

FieldSet unpack_vector_fields(const FieldSet& fields, FieldSet unpackedFields) {
    for (const auto& field : fields) {
        auto componentFieldMetadataVector = std::vector<LocalConfiguration>{};
        if (!field.metadata().get("component_field_metadata", componentFieldMetadataVector)) {
            // Not a vector field.
            addOrReplaceField(unpackedFields, field);
            continue;
        }

        // Field is vector.
        const auto& vectorField = field;

        auto vectorIndex = 0;
        for (const auto& componentFieldMetadata : componentFieldMetadataVector) {
            // Get or create field.
            auto componentFieldName = std::string{};
            componentFieldMetadata.get("name", componentFieldName);
            const auto componentFieldConfig = option::name(componentFieldName) | option::levels(vectorField.levels()) |
                                              option::datatype(vectorField.datatype());
            auto& componentField = getOrCreateField(unpackedFields, vectorField.functionspace(), componentFieldConfig);

            // Copy field data.
            const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
                componentElem = vectorElem(vectorIndex);
            };
            copyFieldData(componentField, vectorField, copier);

            // Copy metadata.
            componentField.metadata() = componentFieldMetadata;
            componentField.set_dirty(vectorField.dirty());

            ++vectorIndex;
        }
    }
    return unpackedFields;
}

}  // namespace util
}  // namespace atlas
