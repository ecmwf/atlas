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
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "eckit/config/LocalConfiguration.h"

namespace atlas {
namespace util {

namespace {

using eckit::LocalConfiguration;

using array::DataType;
using array::helpers::arrayForEachDim;

void addOrReplaceField(FieldSet& fieldSet, const Field& field) {
  const auto fieldName = field.name();
  if (fieldSet.has(fieldName)) {
    fieldSet[fieldName] = field;
  } else {
    fieldSet.add(field);
  }
}

Field& getOrCreateField(FieldSet& fieldSet, const FunctionSpace& functionSpace,
                        const Config& config) {
  const auto fieldName = config.getString("name");
  if (fieldSet.has(fieldName)) {
    const auto field = fieldSet[fieldName];
  } else {
    fieldSet.add(functionSpace.createField(config));
  }
  return fieldSet[fieldName];
}

void checkFieldCompatibility(const Field& componentField,
                             const Field& vectorField) {
  ATLAS_ASSERT(componentField.functionspace().size() ==
               vectorField.functionspace().size());
  ATLAS_ASSERT(componentField.levels() == vectorField.levels());
  ATLAS_ASSERT(componentField.variables() == 0);
  ATLAS_ASSERT(vectorField.variables() > 0);

  const auto checkStandardShape = [](const Field& field) {
    // Check for "standard" Atlas field shape.
    auto dim = 0;
    const auto rank = field.rank();
    const auto shape = field.shape();
    if (field.functionspace().size() != shape[dim++]) {
      return false;
    }
    if (const auto levels = field.levels();
        levels && (dim >= rank || levels != shape[dim++])) {
      return false;
    }
    if (const auto variables = field.variables();
        variables && (dim >= rank || variables != shape[dim++])) {
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
void copyFieldData(ComponentField& componentField, VectorField& vectorField,
                const Functor& copier) {
  checkFieldCompatibility(componentField, vectorField);

  const auto copyArrayData = [&](auto value, auto rank) {
    // Resolve value-type and rank from arguments.
    using Value = decltype(value);
    constexpr auto Rank = decltype(rank)::value;

    // Iterate over fields.
    auto vectorView = array::make_view<Value, Rank>(vectorField);
    auto componentView = array::make_view<Value, Rank - 1>(componentField);
    constexpr auto Dims = std::make_integer_sequence<int, Rank - 1>{};
    arrayForEachDim(Dims, execution::par, std::tie(componentView, vectorView),
                    copier);
  };

  const auto selectRank = [&](auto value) {
    switch (vectorField.rank()) {
      case 2:
        return copyArrayData(value, std::integral_constant<int, 2>{});
      case 3:
        return copyArrayData(value, std::integral_constant<int, 3>{});
      default:
        ATLAS_THROW_EXCEPTION("Unsupported vector field rank: " +
                              std::to_string(vectorField.rank()));
    }
  };

  const auto selectType = [&]() {
    switch (vectorField.datatype().kind()) {
      case DataType::kind<double>():
        return selectRank(double{});
      case DataType::kind<float>():
        return selectRank(float{});
      case DataType::kind<long>():
        return selectRank(long{});
      case DataType::kind<int>():
        return selectRank(int{});
      default:
        ATLAS_THROW_EXCEPTION("Unknown datatype: " +
                              std::to_string(vectorField.datatype().kind()));
    }
  };

  selectType();
}

}  // namespace

namespace pack_vector_fields {

FieldSet pack(const FieldSet& fields, FieldSet packedFields) {
  // Get the number of variables for each vector field.
  auto vectorSizeMap = std::map<std::string, idx_t>{};
  for (const auto& field : fields) {
    auto vectorFieldName = std::string{};
    if (field.metadata().get("vector field name", vectorFieldName)) {
      ++vectorSizeMap[vectorFieldName];
    }
  }
  auto vectorIndexMap = std::map<std::string, idx_t>{};

  // Pack vector fields.
  for (const auto& field : fields) {
    auto vectorFieldName = std::string{};
    if (!field.metadata().get("vector field name", vectorFieldName)) {
      // Not a vector component field.
      addOrReplaceField(packedFields, field);
      continue;
    }

    // Field is vector field component.
    const auto& componentField = field;

    // Get or create vector field.
    const auto vectorFieldConfig =
        option::name(vectorFieldName) |
        option::levels(componentField.levels()) |
        option::variables(vectorSizeMap[vectorFieldName]) |
        option::datatype(componentField.datatype());
    auto& vectorField = getOrCreateField(
        packedFields, componentField.functionspace(), vectorFieldConfig);

    // Copy field data.
    const auto vectorIndex = vectorIndexMap[vectorFieldName]++;
    const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
      vectorElem(vectorIndex) = componentElem;
    };
    copyFieldData(componentField, vectorField, copier);

    // Copy metadata.
    const auto componentFieldMetadata = componentField.metadata();
    auto componentFieldMetadataVector = std::vector<LocalConfiguration>{};
    vectorField.metadata().get("component field metadata",
                               componentFieldMetadataVector);
    componentFieldMetadataVector.push_back(componentFieldMetadata);
    vectorField.metadata().set("component field metadata",
                               componentFieldMetadataVector);
  }
  return packedFields;
}

FieldSet unpack(const FieldSet& fields, FieldSet unpackedFields) {
  for (const auto& field : fields) {
    auto componentFieldMetadataVector = std::vector<LocalConfiguration>{};
    if (!field.metadata().get("component field metadata",
                              componentFieldMetadataVector)) {
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
      const auto componentFieldConfig =
          option::name(componentFieldName) |
          option::levels(vectorField.levels()) |
          option::datatype(vectorField.datatype());
      auto& componentField = getOrCreateField(
          unpackedFields, vectorField.functionspace(), componentFieldConfig);

      // Copy field data.
      const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
        componentElem = vectorElem(vectorIndex);
      };
      copyFieldData(componentField, vectorField, copier);

      // Copy metadata.
      componentField.metadata() = componentFieldMetadata;

      ++vectorIndex;
    }
  }
  return unpackedFields;
}

}  // namespace pack_vector_fields
}  // namespace util
}  // namespace atlas
