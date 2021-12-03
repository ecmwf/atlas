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
 * \brief State class that owns a collection of fields
 */
class State : public util::Object {
public:  // methods
         //-- Constructors
    State();

    State(const std::string& generator, const eckit::Parametrisation& = util::Config());

    //-- Accessors

    const Field& field(const std::string& name) const;
    Field& field(const std::string& name);
    bool has(const std::string& name) const { return (fields_.find(name) != fields_.end()); }
    std::vector<std::string> field_names() const;

    const Field& field(const idx_t idx) const;
    Field& field(const idx_t idx);
    idx_t size() const { return static_cast<idx_t>(fields_.size()); }

    const Field& operator[](const idx_t idx) const { return field(idx); }
    Field& operator[](const idx_t idx) { return field(idx); }

    const Field& operator[](const std::string& name) const { return field(name); }
    Field& operator[](const std::string& name) { return field(name); }

    const util::Metadata& metadata() const;
    util::Metadata& metadata();

    // -- Modifiers

    void initialize(const std::string& generator, const eckit::Parametrisation& = util::Config());

    Field add(Field);

    void remove(const std::string& name);

private:
    typedef std::map<std::string, Field> FieldMap;

private:
    FieldMap fields_;
    util::Metadata metadata_;
};

//------------------------------------------------------------------------------------------------------

class StateGenerator : public util::Object {
public:
    StateGenerator(const eckit::Parametrisation& = util::Config());

    virtual ~StateGenerator();

    virtual void generate(State&, const eckit::Parametrisation& = util::Config()) const = 0;
};

//------------------------------------------------------------------------------------------------------

class StateGeneratorFactory {
public:
    /*!
   * \brief build StateCreator with options specified in parametrisation
   * \return mesh generator
   */
    static StateGenerator* build(const std::string& state_generator, const eckit::Parametrisation& = util::Config());

    /*!
   * \brief list all registered field creators
   */
    static void list(std::ostream&);
    static bool has(const std::string& name);

private:
    virtual StateGenerator* make(const eckit::Parametrisation& = util::Config()) = 0;

    std::string name_;

protected:
    StateGeneratorFactory(const std::string&);
    virtual ~StateGeneratorFactory();
};

template <class T>
class StateGeneratorBuilder : public StateGeneratorFactory {
    virtual StateGenerator* make(const eckit::Parametrisation& param = util::Config()) { return new T(param); }

public:
    StateGeneratorBuilder(const std::string& name): StateGeneratorFactory(name) {}
};

// ------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
State* atlas__State__new();
void atlas__State__initialize(State* This, const char* generator, const eckit::Parametrisation* params);
void atlas__State__delete(State* This);
void atlas__State__add(State* This, FieldImpl* field);
void atlas__State__remove(State* This, const char* name);
int atlas__State__has(State* This, const char* name);
FieldImpl* atlas__State__field_by_name(State* This, const char* name);
FieldImpl* atlas__State__field_by_index(State* This, idx_t index);
idx_t atlas__State__size(const State* This);
util::Metadata* atlas__State__metadata(State* This);
}

}  // namespace field
}  // namespace atlas
