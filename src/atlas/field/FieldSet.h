/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file FieldSet.h
/// @author Willem Deconinck
/// @author Pedro Maciel
/// @date Jan 2015

#pragma once

#include <cstdlib>
#include <iterator>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "eckit/deprecated.h"

#include "atlas/array_fwd.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/field/detail/FieldImpl.h"

namespace atlas {

class FieldSet;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace field {

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSetImpl : public util::Object {
public:  // types
    using iterator       = std::vector<Field>::iterator;
    using const_iterator = std::vector<Field>::const_iterator;

    template <typename T>
    static constexpr bool is_index() {
        return std::is_integral<T>::value or std::is_enum<T>::value;
    }

    template <bool pred>
    using enable_if_t = typename std::enable_if<pred, int>::type;

    template <typename T>
    using enable_if_index_t = enable_if_t<is_index<T>()>;

private:

    class FieldObserver : public field::FieldObserver {
    public:
        FieldObserver(FieldSetImpl& fieldset) : fieldset_(fieldset) {}

    private:
        void onFieldRename(FieldImpl& field) override;

    private:
        FieldSetImpl& fieldset_;
    };


public:  // methods
    /// Constructs an empty FieldSet
    FieldSetImpl(const std::string& name = "untitled");
    virtual ~FieldSetImpl();

    idx_t size() const { return static_cast<idx_t>(fields_.size()); }
    bool empty() const { return !fields_.size(); }

    void clear();

    const std::string& name() const { return name_; }
    std::string& name() { return name_; }

    template <typename Index, enable_if_index_t<Index> = 0>
    const Field& operator[](Index i) const {
        return field(i);
    }

    template <typename Index, enable_if_index_t<Index> = 0>
    Field& operator[](Index i) {
        return field(i);
    }

    const Field& operator[](const std::string& name) const { return field(name); }
    Field& operator[](const std::string& name) { return field(name); }

    template <typename Index, enable_if_index_t<Index> = 0>
    const Field& field(Index i) const {
        if (i >= size())
            throw_OutOfRange("fieldset", i, size(), Here());
        return fields_[i];
    }

    template <typename Index, enable_if_index_t<Index> = 0>
    Field& field(Index i) {
        if (i >= size())
            throw_OutOfRange("fieldset", i, size(), Here());
        return fields_[i];
    }

    const std::vector<std::string>& field_names() const;

    Field add(const Field&);
    void add(const FieldSet&);

    bool has(const std::string& name) const;

    Field& field(const std::string& name) const;

    iterator begin() { return fields_.begin(); }
    iterator end() { return fields_.end(); }
    const_iterator begin() const { return fields_.begin(); }
    const_iterator end() const { return fields_.end(); }
    const_iterator cbegin() const { return fields_.begin(); }
    const_iterator cend() const { return fields_.end(); }

    const util::Metadata& metadata() const { return metadata_; }
    util::Metadata& metadata() { return metadata_; }

    void haloExchange(bool on_device = false) const;
    void adjointHaloExchange(bool on_device = false) const;
    void set_dirty(bool = true) const;

    void updateHost() const;
    void updateHost(std::initializer_list<std::string> names) const;
    void updateHost(std::initializer_list<int> indices) const;
    void updateDevice() const;
    void updateDevice(std::initializer_list<std::string> names) const;
    void updateDevice(std::initializer_list<int> indices) const;
    void allocateDevice() const;
    void allocateDevice(std::initializer_list<std::string> names) const;
    void allocateDevice(std::initializer_list<int> indices) const;
    void deallocateDevice() const;
    void deallocateDevice(std::initializer_list<std::string> names) const;
    void deallocateDevice(std::initializer_list<int> indices) const;



protected:                                // data
    std::vector<Field> fields_;           ///< field storage
    std::string name_;                    ///< internal name
    util::Metadata metadata_;             ///< metadata associated with the FieldSet
    std::map<std::string, idx_t> index_;  ///< name-to-index map, to refer fields by name
    std::vector<std::string>     field_names_;  ///< field names
    std::map<std::string, idx_t> duplicates_;  ///< name-to-duplicates map, to refer fields by name

    friend class FieldObserver;
    FieldObserver field_observer_;
};

// C wrapper interfaces to C++ routines
extern "C" {
FieldSetImpl* atlas__FieldSet__new(char* name);
void atlas__FieldSet__delete(FieldSetImpl* This);
void atlas__FieldSet__add_field(FieldSetImpl* This, FieldImpl* field);
void atlas__FieldSet__add_fieldset(FieldSetImpl* This, FieldSetImpl* fieldset);
int atlas__FieldSet__has_field(const FieldSetImpl* This, char* name);
const char* atlas__FieldSet__name(FieldSetImpl* This);
idx_t atlas__FieldSet__size(const FieldSetImpl* This);
FieldImpl* atlas__FieldSet__field_by_name(FieldSetImpl* This, char* name);
FieldImpl* atlas__FieldSet__field_by_idx(FieldSetImpl* This, idx_t idx);
void atlas__FieldSet__data_int_specf(FieldSetImpl* This, char* name, int*& field_data, int& rank, int*& field_shapef,
                                  int*& field_stridesf);
void atlas__FieldSet__data_long_specf(FieldSetImpl* This, char* name, long*& field_data, int& rank, int*& field_shapef,
                                   int*& field_stridesf);
void atlas__FieldSet__data_float_specf(FieldSetImpl* This, char* name, float*& field_data, int& rank, int*& field_shapef,
                                    int*& field_stridesf);
void atlas__FieldSet__data_double_specf(FieldSetImpl* This, char* name, double*& field_data, int& rank, int*& field_shapef,
                                     int*& field_stridesf);
void atlas__FieldSet__data_int_specf_by_idx(FieldSetImpl* This, int& idx, int*& field_data, int& rank, int*& field_shapef,
                                  int*& field_stridesf);
void atlas__FieldSet__data_long_specf_by_idx(FieldSetImpl* This, int& idx, long*& field_data, int& rank, int*& field_shapef,
                                   int*& field_stridesf);
void atlas__FieldSet__data_float_specf_by_idx(FieldSetImpl* This, int& idx, float*& field_data, int& rank, int*& field_shapef,
                                    int*& field_stridesf);
void atlas__FieldSet__data_double_specf_by_idx(FieldSetImpl* This, int& idx, double*& field_data, int& rank, int*& field_shapef,
                                     int*& field_stridesf);
void atlas__FieldSet__set_dirty(FieldSetImpl* This, int value);
void atlas__FieldSet__halo_exchange(FieldSetImpl* This, int on_device);
}

}  // namespace field
#endif

//---------------------------------------------------------------------------------------------------------------------

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSet : DOXYGEN_HIDE(public util::ObjectHandle<field::FieldSetImpl>) {
public:  // types
    using iterator       = Implementation::iterator;
    using const_iterator = Implementation::const_iterator;

    template <typename T>
    using enable_if_index_t = Implementation::enable_if_index_t<T>;

public:  // methods
    using Handle::Handle;
    FieldSet();
    FieldSet(const std::string& name);
    FieldSet(const Field&);

    FieldSet clone(const eckit::Parametrisation& config = util::Config()) const;

    idx_t size() const { return get()->size(); }
    bool empty() const { return get()->empty(); }

    void clear() { get()->clear(); }

    const std::string& name() const { return get()->name(); }
    std::string& name() { return get()->name(); }

    template <typename Index, enable_if_index_t<Index> = 0>
    const Field& operator[](Index i) const {
        return get()->operator[](i);
    }

    template <typename Index, enable_if_index_t<Index> = 0>
    Field& operator[](Index i) {
        return get()->operator[](i);
    }

    const Field& operator[](const std::string& name) const { return get()->operator[](name); }
    Field& operator[](const std::string& name) { return get()->operator[](name); }

    const Field& operator[](const char* name) const { return get()->operator[](name); }
    Field& operator[](const char* name) { return get()->operator[](name); }

    template <typename Index, enable_if_index_t<Index> = 0>
    const Field& field(Index i) const {
        return get()->field(i);
    }

    template <typename Index, enable_if_index_t<Index> = 0>
    Field& field(Index i) {
        return get()->field(i);
    }

    const std::vector<std::string>& field_names() const { return get()->field_names(); }

    Field add(const Field& field) { return get()->add(field); }
    void add(const FieldSet& fieldset) { return get()->add(fieldset); }

    bool has(const std::string& name) const { return get()->has(name); }

    Field& field(const std::string& name) const { return get()->field(name); }

    iterator begin() { return get()->begin(); }
    iterator end() { return get()->end(); }
    const_iterator begin() const { return get()->begin(); }
    const_iterator end() const { return get()->end(); }
    const_iterator cbegin() const { return get()->begin(); }
    const_iterator cend() const { return get()->end(); }

    const util::Metadata& metadata() const;
    util::Metadata& metadata();

    void set_dirty(bool = true) const;

    void haloExchange(bool on_device = false) const { get()->haloExchange(on_device); }
    void adjointHaloExchange(bool on_device = false) const { get()->adjointHaloExchange(on_device); }

    // -- Methods related to host-device synchronisation
    void updateHost() const { get()->updateHost(); }
    void updateHost(std::initializer_list<std::string> names) const { get()->updateHost(names); }
    void updateHost(std::initializer_list<int> indices) const { get()->updateHost(indices); }
    void updateDevice() const { get()->updateDevice(); }
    void updateDevice(std::initializer_list<std::string> names) const { get()->updateDevice(names); }
    void updateDevice(std::initializer_list<int> indices) const{ get()->updateDevice(indices); }
    void allocateDevice() const { get()->allocateDevice(); }
    void allocateDevice(std::initializer_list<std::string> names) const { get()->allocateDevice(names); }
    void allocateDevice(std::initializer_list<int> indices) const { get()->allocateDevice(indices); }
    void deallocateDevice() const { get()->deallocateDevice(); }
    void deallocateDevice(std::initializer_list<std::string> names) const { get()->deallocateDevice(names); }
    void deallocateDevice(std::initializer_list<int> indices) const { get()->deallocateDevice(indices); }

    // Deprecated API
    DEPRECATED("use 'has' instead") bool has_field(const std::string& name) const { return get()->has(name); }
};

}  // namespace atlas
