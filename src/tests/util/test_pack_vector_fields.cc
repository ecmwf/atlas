/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/util/PackVectorFields.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

FieldSet setFields(const FunctionSpace& functionSpace,
                   const std::vector<util::Config>& fieldConfigs) {
  auto fields = FieldSet{};

  // Set unique values to all field elements.
  auto value = 0;
  for (const auto& fieldConfig : fieldConfigs) {
    auto field = fields.add(functionSpace.createField(fieldConfig));
    for (auto arrayIdx = size_t{0}; arrayIdx < field.size(); arrayIdx++) {
      field->data<float>()[arrayIdx] = value++;
    }
    field.metadata().set("comment", "This field is made with love.");
    auto vectorFieldName = std::string{};
    if (fieldConfig.get("vector_field_name", vectorFieldName)) {
      field.metadata().set("vector_field_name", vectorFieldName);
    }
  }
  return fields;
}

FieldSet createOrderedTestFields() {
  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 are contiguous in field set.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("scalar") | option::levels(1) |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::levels(1) |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::levels(1) |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createUnorderedTestFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 are not contiguous in field set.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createInconsistentRankFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 have different ranks.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::levels(10) |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createInconsistentDatatypeFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 have different datatypes.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::datatype(DataType::kind<int>()) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createInconsistentLevelsFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 have different number of levels.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::levels(10) |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::levels(20) |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createInconsistentVariablesFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 have different number of variables.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector_component_0") |
                         option::datatype(DataType::kind<float>()) |
                         option::variables(2) |
                         util::Config{"vector_field_name", "vector"});
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector_component_1") |
                         option::datatype(DataType::kind<float>()) |
                         util::Config{"vector_field_name", "vector"});

  return setFields(functionSpace, fieldConfigs);
}

void checkTestFields(const FieldSet& fields) {
  auto value = 0;
  for (const auto& field : fields) {
    for (auto arrayIdx = size_t{0}; arrayIdx < field.size(); arrayIdx++) {
      EXPECT(field->data<float>()[arrayIdx] == value++);
    }
    EXPECT(field.metadata().get<std::string>("comment") ==
           "This field is made with love.");
  }
}

CASE("Basic pack and unpack") {
  const auto fields = createOrderedTestFields();

  const auto packedFields = util::pack_vector_fields(fields);

  EXPECT(!packedFields.has("vector_component_0"));
  EXPECT(!packedFields.has("vector_component_1"));
  EXPECT(packedFields.has("vector"));
  EXPECT(packedFields.has("scalar"));

  const auto unpackedFields = util::unpack_vector_fields(packedFields);

  EXPECT(unpackedFields.has("vector_component_0"));
  EXPECT(unpackedFields.has("vector_component_1"));
  EXPECT(!unpackedFields.has("vector"));
  EXPECT(unpackedFields.has("scalar"));

  checkTestFields(unpackedFields);
}

CASE("unpack into existing field set") {
  auto fields = createUnorderedTestFields();

  const auto packedFields = util::pack_vector_fields(fields);

  EXPECT(!packedFields.has("vector_component_0"));
  EXPECT(!packedFields.has("vector_component_1"));
  EXPECT(packedFields.has("vector"));
  EXPECT(packedFields.has("scalar"));

  // Need to unpack into existing field to guarantee field order is preserved.
  array::make_view<float, 1>(fields["vector_component_0"]).assign(0.);
  array::make_view<float, 1>(fields["vector_component_1"]).assign(0.);
  util::unpack_vector_fields(packedFields, fields);

  EXPECT(fields.has("vector_component_0"));
  EXPECT(fields.has("vector_component_1"));
  EXPECT(!fields.has("vector"));
  EXPECT(fields.has("scalar"));

  checkTestFields(fields);
}

CASE("check that bad inputs throw") {
  // Try to apply pack to inconsistent field sets.
  EXPECT_THROWS(util::pack_vector_fields(createInconsistentRankFields()));
  EXPECT_THROWS(
      util::pack_vector_fields(createInconsistentDatatypeFields()));
  EXPECT_THROWS(
      util::pack_vector_fields(createInconsistentLevelsFields()));
  EXPECT_THROWS(
      util::pack_vector_fields(createInconsistentVariablesFields()));
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
