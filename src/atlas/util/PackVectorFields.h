#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/util/Config.h"

namespace eckit {
class LocalConfiguration;
}

namespace atlas {
namespace util {

using eckit::LocalConfiguration;
using LocalConfigurations = std::vector<LocalConfiguration>;

class PackVectorFields {
 public:
  PackVectorFields(const LocalConfiguration& config);
  FieldSet pack(const FieldSet& fields,
                FieldSet packedFields = FieldSet{}) const;
  FieldSet unpack(const FieldSet& fields,
                  FieldSet unpackedFields = FieldSet{}) const;

 private:
  Config componentToVectorMap_{};
  Config vectorToComponentMap_{};
};

}  // namespace util
}  // namespace atlas
