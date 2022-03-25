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

#include "atlas/array/Array.h"
#include "atlas/field/Field.h"
#include "atlas/util/Config.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Object.h"

namespace eckit {
class Parametrisation;
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
class MultiField : public util::Object {
public:  // methods
         //-- Constructors
    MultiField();

    MultiField(const std::string& generator, const eckit::Parametrisation& = util::Config());

    //-- Accessors

    const Field& field(const std::string& name) const { return fields_[field_index_.at(name)]; }
    Field& field(const std::string& name) { return fields_[field_index_.at(name)]; }
    bool has(const std::string& name) const { return (field_index_.find(name) != field_index_.end()); }
    std::vector<std::string> field_names() const;

    const Field& field(const idx_t idx) const { return fields_[idx]; }
    Field& field(const idx_t idx) { return fields_[idx]; }
    idx_t size() const { return static_cast<idx_t>(fields_.size()); }

    const Field& operator[](const idx_t idx) const { return fields_[idx]; }
    Field& operator[](const idx_t idx) { return fields_[idx]; }

    const Field& operator[](const std::string& name) const { return field(name); }
    Field& operator[](const std::string& name) { return field(name); }

    const util::Metadata& metadata() const;
    util::Metadata& metadata();

    // -- Modifiers

    void initialize(const std::string& generator, const eckit::Parametrisation& = util::Config());

    array::Array& allocate(array::DataType datatype, array::ArraySpec&&);

    /// @brief Implicit conversion to Array
    operator const array::Array&() const { return *array_; }
    operator array::Array&() { return *array_; }

    /// @brief Access contained Array
    const array::Array& array() const { return *array_; }
    array::Array& array() { return *array_; }


public:  // temporary public for prototyping
    std::map<std::string, int> field_index_;
    std::vector<Field> fields_;
    std::unique_ptr<array::Array> array_;
    util::Metadata metadata_;
};

//------------------------------------------------------------------------------------------------------

class MultiFieldCreator : public util::Object {
public:
    MultiFieldCreator(const eckit::Parametrisation& = util::Config());

    virtual ~MultiFieldCreator();

    virtual void generate(MultiField&, const eckit::Parametrisation& = util::Config()) const = 0;
};

//------------------------------------------------------------------------------------------------------

class MultiFieldCreatorFactory {
public:
    /*!
   * \brief build FieldPoolCreator with options specified in parametrisation
   * \return mesh generator
   */
    static MultiFieldCreator* build(const std::string& FieldPool_generator,
                                    const eckit::Parametrisation& = util::Config());

    /*!
   * \brief list all registered field creators
   */
    static void list(std::ostream&);
    static bool has(const std::string& name);

private:
    virtual MultiFieldCreator* make(const eckit::Parametrisation& = util::Config()) = 0;

    std::string name_;

protected:
    MultiFieldCreatorFactory(const std::string&);
    virtual ~MultiFieldCreatorFactory();
};

template <class T>
class MultiFieldCreatorBuilder : public MultiFieldCreatorFactory {
    virtual MultiFieldCreator* make(const eckit::Parametrisation& param = util::Config()) { return new T(param); }

public:
    MultiFieldCreatorBuilder(const std::string& name): MultiFieldCreatorFactory(name) {}
};

// ------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
MultiField* atlas__FieldPool__new();
void atlas__FieldPool__initialize(MultiField* This, const char* generator, const eckit::Parametrisation* params);
void atlas__FieldPool__delete(MultiField* This);
int atlas__FieldPool__has(MultiField* This, const char* name);
FieldImpl* atlas__FieldPool__field_by_name(MultiField* This, const char* name);
FieldImpl* atlas__FieldPool__field_by_index(MultiField* This, idx_t index);
idx_t atlas__FieldPool__size(const MultiField* This);
util::Metadata* atlas__FieldPool__metadata(MultiField* This);
}

}  // namespace field
}  // namespace atlas
