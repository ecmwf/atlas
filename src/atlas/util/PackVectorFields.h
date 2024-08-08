#pragma once

#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"

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
  FieldSet pack(const FieldSet& fields) const;
  FieldSet unpack(const FieldSet& fields) const;

 private:
  using VectorName = std::string;
  using ComponentName = std::string;
  struct ComponentInfo {
    VectorName vectorName{};
    idx_t componentIndex{};
    idx_t vectorSize{};
  };
  struct VectorInfo {
    std::vector<ComponentName> componentNameList{};
  };

  using PackingMap = std::map<ComponentName, ComponentInfo>;
  using UnpackingMap = std::map<VectorName, VectorInfo>;

  PackingMap componentToVectorMap_{};
  UnpackingMap vectorToComponentMap_{};
};

}  // namespace util
}  // namespace atlas
