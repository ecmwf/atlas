#pragma once

#include "atlas/field.h"

namespace eckit {
class LocalConfiguration;
}

namespace atlas {
namespace util {

namespace pack_vector_fields {

/// @brief   Packs vector field components into vector fields
///
/// @details Iterates through @param fields and creates vector fields from any
///          component field with the "vector name" string metadata. These, as
///          well as any present scalar fields are added to the return-value
///          field set.
///          Note, a mutable @param packedFields field set can be supplied if
///          one needs to guarantee the order of the packed fields
FieldSet pack(const FieldSet& fields, FieldSet packedFields = FieldSet{});

/// @brief   Unpacks vector field into vector field components.
///
/// @details Undoes "pack" operation when a set of packed fields are supplied
///          as @param fields. A mutable @param unpackedFields field set can be
///          supplied if one needs to guarantee the order of the unpacked
///          fields.
FieldSet unpack(const FieldSet& fields, FieldSet unpackedFields = FieldSet{});
};  // namespace pack_vector_fields

}  // namespace util
}  // namespace atlas
