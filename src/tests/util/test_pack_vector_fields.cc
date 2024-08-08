#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/PackVectorFields.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

FieldSet createTestFields() {
  auto fields = FieldSet{};

  const auto grid = Grid("O16");
  const auto functionSpace = functionspace::StructuredColumns(grid);

  fields.add(functionSpace.createField<double>(option::name("scalar") | option::levels(1)));
  fields.add(functionSpace.createField<double>(option::name("vector component 0") | option::levels(1)));
  fields.add(functionSpace.createField<double>(option::name("vector component 1") | option::levels(1)));

  return fields;
}

CASE("Test pack and unpack (Rank 1)") {
  using util::Config;
  using util::PackVectorFields;

  using ConfVec = std::vector<Config>;
  using StrVec = std::vector<std::string>;

  const auto windComponents = StrVec{"vector component 0", "vector component 1"};
  const auto windConfig =
      Config("name", "vector") | Config("components", windComponents);
  const auto vectorFieldsConf = Config("vector fields", ConfVec{windConfig});

  const auto fieldPackager = PackVectorFields(vectorFieldsConf);

  const auto fields = createTestFields();

  const auto packedFields = fieldPackager.pack(fields);

  for (const auto& field : packedFields) {
    Log::info() << field.shape(0) << std::endl
                << field.name() << std::endl
                << field.levels() << std::endl
                << field.variables() << std::endl
                << field.rank() << std::endl;
  }
  Log::info() << packedFields["vector"].metadata() << std::endl;
}
}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
