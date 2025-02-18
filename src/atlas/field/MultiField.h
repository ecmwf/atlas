/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date June 2015

#pragma once

#include <map>
#include <vector>

#include "atlas/array/Array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace field {
class MultiFieldImpl;
}
}

namespace atlas {
namespace field {

/**
 * \brief MultiField class that owns a collection of fields that are co-allocated
 *
 * Fields can only be described by parametrisation during the construction.
 * Once setup, no additional fields can be added.
 *
 * Fields have to all be of same memory layout and data type
 */

class MultiField : public util::ObjectHandle<MultiFieldImpl> {
public:  // methods
         //-- Constructors
    using Handle::Handle;

    MultiField(const eckit::Configuration&);
    MultiField(const array::DataType datatype, const array::ArrayShape& shape,
        const std::vector<std::string>& var_names);

    //-- Accessors

    const Field& field(const std::string& name) const;
    Field& field(const std::string& name);
    bool has(const std::string& name) const;
    std::vector<std::string> field_names() const;

    const Field& field(const idx_t idx) const;
    Field& field(const idx_t idx);
    idx_t size() const;

    const Field& operator[](const idx_t idx) const;
    Field& operator[](const idx_t idx);

    const Field& operator[](const std::string& name) const;
    Field& operator[](const std::string& name);

    const util::Metadata& metadata() const;
    util::Metadata& metadata();

    // -- Modifiers

    /// @brief Implicit conversion to Array
    operator const array::Array&() const;
    operator array::Array&();

    operator const FieldSet&() const;
    operator FieldSet&();

    /// @brief Access contained Array
    const array::Array& array() const;
    array::Array& array();

private:
    template<typename datatype>
    void create(const std::vector<int> shape, const std::vector<std::string> var_names);
};

/**
 * \brief MultiFieldArrayRegistry
 */

class MultiFieldArrayRegistry : public field::FieldObserver {
private:
    MultiFieldArrayRegistry() {}

public:
    static MultiFieldArrayRegistry& instance() {
        static MultiFieldArrayRegistry inst;
        return inst;
    }
    void onFieldDestruction(FieldImpl& field) override {
        std::lock_guard<std::mutex> guard(lock_);
        map_.erase(&field);
    }

    ~MultiFieldArrayRegistry() override = default;

    void add(Field& field, std::shared_ptr<array::Array> array) {
        std::lock_guard<std::mutex> guard(lock_);
        map_.emplace(field.get(), array);
        field->attachObserver(*this);
    }

public:
    std::mutex lock_;
    std::map<FieldImpl*,std::shared_ptr<array::Array>> map_;

};

// ------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
