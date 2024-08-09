#include "atlas/util/PackVectorFields.h"

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

template <int I>
using Integer = std::integral_constant<int, I>;
template <int N>
using IntegerSequence = std::make_integer_sequence<int, N>;

using array::DataType;
using array::helpers::arrayForEachDim;

template <typename T>
T getT(const Config& config, const std::string& key) {
  auto value = T{};
  const auto success = config.get(key, value);
  ATLAS_ASSERT_MSG(success,
                   "\"" + key + "\" entry missing from configuration.");
  return value;
}

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
void copyFields(ComponentField& componentField, VectorField& vectorField,
                const Functor& copier) {
  checkFieldCompatibility(componentField, vectorField);

  const auto copyData = [&](auto value, auto rank) {
    // Resolve value-type and rank from arguments.
    using Value = decltype(value);
    constexpr auto Rank = decltype(rank)::value;

    // Iterate over fields.
    auto vectorView = array::make_view<Value, Rank>(vectorField);
    auto componentView = array::make_view<Value, Rank - 1>(componentField);
    constexpr auto Dims = IntegerSequence<Rank - 1>{};
    arrayForEachDim(Dims, execution::par, std::tie(componentView, vectorView),
                    copier);
  };

  const auto selectRank = [&](auto value) {
    switch (vectorField.rank()) {
      case 2:
        return copyData(value, Integer<2>{});
      case 3:
        return copyData(value, Integer<3>{});
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

PackVectorFields::PackVectorFields(const LocalConfiguration& config) {
  const auto vectorFieldList =
      getT<LocalConfigurations>(config, "vector fields");

  for (const auto& vectorField : vectorFieldList) {
    const auto vectorName = vectorField.getString("name");
    const auto componentNameList = vectorField.getStringVector("components");

    // Vector information for each component field.
    auto componentIndex = idx_t{0};
    const auto vectorSize = componentNameList.size();
    for (const auto& componentName : componentNameList) {
      const auto vectorFieldConfig =
          Config("component index", componentIndex++) |
          option::name(vectorName) | option::variables(vectorSize);
      componentToVectorMap_.set(componentName, vectorFieldConfig);
    }

    // Component information for each vector field.
    auto componentFieldConfigs = std::vector<Config>{};
    for (const auto& componentName : componentNameList) {
      componentFieldConfigs.push_back(option::name(componentName));
    }
    vectorToComponentMap_.set(vectorName, componentFieldConfigs);
  }
}

FieldSet PackVectorFields::pack(const FieldSet& fields,
                                FieldSet packedFields) const {
  for (const auto& field : fields) {
    if (!componentToVectorMap_.has(field.name())) {
      // Not a vector component field.
      addOrReplaceField(packedFields, field);
      continue;
    }

    // Field is vector field component.
    const auto& componentField = field;

    // Get or create vector field.
    const auto vectorFieldConfig =
        getT<Config>(componentToVectorMap_, componentField.name()) |
        option::levels(componentField.levels()) |
        option::datatype(componentField.datatype());
    auto& vectorField = getOrCreateField(
        packedFields, componentField.functionspace(), vectorFieldConfig);

    // Copy field data.
    const auto componentIndex =
        getT<idx_t>(vectorFieldConfig, "component index");
    const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
      vectorElem(componentIndex) = componentElem;
    };
    copyFields(componentField, vectorField, copier);

    // Copy metadata.
    const auto componentMetadata = componentField.metadata();
    vectorField.metadata().set(componentField.name() + " metadata",
                               componentMetadata);
  }
  return packedFields;
}

FieldSet PackVectorFields::unpack(const FieldSet& fields,
                                  FieldSet unpackedFields) const {
  for (const auto& field : fields) {
    if (!vectorToComponentMap_.has(field.name())) {
      // Not a vector field.
      addOrReplaceField(unpackedFields, field);
      continue;
    }

    // Field is vector.
    const auto& vectorField = field;

    const auto componentFieldConfigs =
        getT<std::vector<Config>>(vectorToComponentMap_, vectorField.name());

    for (auto componentIndex = size_t{0};
         componentIndex < componentFieldConfigs.size();
         ++componentIndex) {
      // Get or create field.
      const auto componentFieldConfig =
          componentFieldConfigs[componentIndex] |
          option::levels(vectorField.levels()) |
          option::datatype(vectorField.datatype());
      auto& componentField = getOrCreateField(
          unpackedFields, vectorField.functionspace(), componentFieldConfig);

      // Copy field data.
      const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
        componentElem = vectorElem(componentIndex);
      };
      copyFields(componentField, vectorField, copier);

      // Copy metadata.
      componentField.metadata() = vectorField.metadata().get<Config>(
          componentField.name() + " metadata");
    }
  }
  return unpackedFields;
}

}  // namespace util
}  // namespace atlas
