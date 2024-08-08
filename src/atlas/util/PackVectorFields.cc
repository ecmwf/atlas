#include "atlas/util/PackVectorFields.h"

#include <string>
#include <type_traits>

#include "atlas/array.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/functionspace.h"
#include "atlas/option.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
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
T getT(const LocalConfiguration& config, const std::string& key) {
  auto value = T{};
  const auto success = config.get(key, value);
  ATLAS_ASSERT_MSG(success,
                   "\"" + key + "\" entry missing from configuration.");
  return value;
}

bool checkShapeCompatibility(const Field& field) {
  // Check for "standard" Atlas field shape.
  auto dim = 0;
  const auto rank = field.rank();
  const auto shape = field.shape();
  if (field.functionspace().size() != shape[dim++]) {
    return false;
  }
  if (const auto levels = field.levels(); levels &&
      (dim >= rank || levels != shape[dim++])) {
    return false;
  }
  if (const auto variables = field.variables(); variables &&
      (dim >= rank || variables != shape[dim++])) {
    return false;
  }
  if (dim != rank) {
    return false;
  }
  return true;
}

Field getOrCreateField(FieldSet& fieldSet, const FunctionSpace& functionSpace,
                       const Config& config) {
  const auto fieldName = config.getString("name");
  if (fieldSet.has(fieldName)) {
    const auto field = fieldSet[fieldName];
    ATLAS_ASSERT(field.levels() == getT<idx_t>(config, "levels"));
    ATLAS_ASSERT(field.variables() == getT<idx_t>(config, "variables"));
    ATLAS_ASSERT(field.datatype() ==
                 getT<DataType::kind_t>(config, "datatype"));
    ATLAS_ASSERT(checkShapeCompatibility(field));
    return field;
  }
  return fieldSet.add(functionSpace.createField(config));
}

template <typename ComponentField, typename VectorField, typename Functor>
void copyFields(ComponentField& componentField, VectorField& vectorField,
                const Functor& copier) {
  const auto copyData = [&](auto value, auto rank) {
    // Resolve value-type and rank from arguments.
    using Value = decltype(value);
    constexpr auto Rank = decltype(rank)::value;

    // Iterate over fields.
    auto vectorView = array::make_view<Value, Rank>(vectorField);
    const auto componentView =
        array::make_view<Value, Rank - 1>(componentField);
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
    idx_t componentIndex = 0;
    const idx_t vectorSize = componentNameList.size();
    for (const auto& componentName : componentNameList) {
      componentToVectorMap_.insert(
          {componentName, {vectorName, componentIndex++, vectorSize}});
    }

    // Component information for each vector field.
    vectorToComponentMap_.insert({vectorName, {componentNameList}});
  }
}

FieldSet PackVectorFields::pack(const FieldSet& fields) const {
  auto packedFields = FieldSet{};
  for (const auto& field : fields) {
    const auto mapItr = componentToVectorMap_.find(field.name());
    if (mapItr == componentToVectorMap_.end()) {
      // Not in the list of vector fields.
      packedFields.add(field);
      continue;
    }
    const auto& componentField = field;
    ATLAS_ASSERT_MSG(!componentField.variables(),
                     "Component field must be scalar.");

    // Set vector field config.
    const auto vectorInfo = mapItr->second;
    const auto vectorFieldConf = option::name(vectorInfo.vectorName) |
                                 option::levels(componentField.levels()) |
                                 option::variables(vectorInfo.vectorSize) |
                                 option::datatype(componentField.datatype());
    auto vectorField = getOrCreateField(
        packedFields, componentField.functionspace(), vectorFieldConf);

    // Copy field data.
    const auto copier = [&](auto&& componentElem, auto&& vectorElem) {
      vectorElem(vectorInfo.componentIndex) = componentElem;
    };
    copyFields(componentField, vectorField, copier);

    // Copy metadata.
    const auto componentMetadata = componentField.metadata();
    vectorField.metadata().set(componentField.name() + " metadata",
                               componentMetadata);
  }
  return packedFields;
}

FieldSet PackVectorFields::unpack(const FieldSet& fields) const {}

}  // namespace util
}  // namespace atlas
