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
  fieldConfigs.push_back(option::name("vector component 0") |
                         option::levels(1) |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector component 1") |
                         option::levels(1) |
                         option::datatype(DataType::kind<float>()));

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createUnorderedTestFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 are not contiguous in field set.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector component 0") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector component 1") |
                         option::datatype(DataType::kind<float>()));

  return setFields(functionSpace, fieldConfigs);
}

FieldSet createInconsistentTestFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  // Note: vector components 0 and 1 have different ranks.
  auto fieldConfigs = std::vector<util::Config>{};
  fieldConfigs.push_back(option::name("vector component 0") |
                         option::levels(10) |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("scalar") |
                         option::datatype(DataType::kind<float>()));
  fieldConfigs.push_back(option::name("vector component 1") |
                         option::datatype(DataType::kind<float>()));

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

util::Config createFieldPackerConfig() {
  using ConfVec = std::vector<util::Config>;
  using StrVec = std::vector<std::string>;

  // yaml equivalent is:
  //
  // vector fields:
  // - name: <name of vector field 1>
  //   components:
  //   - <name of component 1 of vector field 1>
  //   - <name of component 2 of vector field 1>
  //   - <name of component ...>
  // - name: ...
  //   components: ...

  const auto components = StrVec{"vector component 0", "vector component 1"};
  const auto vectorConfig =
      util::Config{"name", "vector"} | util::Config("components", components);
  const auto fieldPackerConfig =
      util::Config{"vector fields", ConfVec{vectorConfig}};

  return fieldPackerConfig;
}

CASE("Basic pack and unpack") {
  using util::PackVectorFields;
  const auto vectorFieldsConf = createFieldPackerConfig();

  const auto fieldPacker = PackVectorFields{vectorFieldsConf};

  const auto fields = createOrderedTestFields();

  const auto packedFields = fieldPacker.pack(fields);

  EXPECT(!packedFields.has("vector component 0"));
  EXPECT(!packedFields.has("vector component 1"));
  EXPECT(packedFields.has("vector"));
  EXPECT(packedFields.has("scalar"));

  const auto unpackedFields = fieldPacker.unpack(packedFields);

  EXPECT(unpackedFields.has("vector component 0"));
  EXPECT(unpackedFields.has("vector component 1"));
  EXPECT(!unpackedFields.has("vector"));
  EXPECT(unpackedFields.has("scalar"));

  checkTestFields(unpackedFields);
}

CASE("unpack into existing field set") {
  using util::PackVectorFields;
  const auto vectorFieldsConf = createFieldPackerConfig();

  const auto fieldPacker = PackVectorFields{vectorFieldsConf};

  auto fields = createUnorderedTestFields();

  const auto packedFields = fieldPacker.pack(fields);

  EXPECT(!packedFields.has("vector component 0"));
  EXPECT(!packedFields.has("vector component 1"));
  EXPECT(packedFields.has("vector"));
  EXPECT(packedFields.has("scalar"));

  // Need to unpack into existing field to guarantee field order is preserved.
  array::make_view<float, 1>(fields["vector component 0"]).assign(0.);
  array::make_view<float, 1>(fields["vector component 1"]).assign(0.);
  fieldPacker.unpack(packedFields, fields);

  EXPECT(fields.has("vector component 0"));
  EXPECT(fields.has("vector component 1"));
  EXPECT(!fields.has("vector"));
  EXPECT(fields.has("scalar"));

  checkTestFields(fields);
}

CASE("check that bad inputs throw") {
  using util::PackVectorFields;

  // Config must have a "vector fields" entry.
  EXPECT_THROWS(PackVectorFields{util::Config{}});
  const auto noPacking = PackVectorFields{
      util::Config{"vector fields", std::vector<util::Config>{}}};

  // These will do nothing.
  const auto fields = createInconsistentTestFields();
  const auto packedFields = noPacking.pack(fields);
  const auto unpackedFields = noPacking.unpack(packedFields);

  checkTestFields(fields);
  checkTestFields(packedFields);
  checkTestFields(unpackedFields);

  // Try to apply packer to inconsistent field set.
  const auto vectorFieldsConf = createFieldPackerConfig();
  const auto fieldPacker = PackVectorFields{vectorFieldsConf};
  EXPECT_THROWS(fieldPacker.pack(fields));

}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
