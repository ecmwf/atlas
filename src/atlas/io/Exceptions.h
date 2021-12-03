/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <type_traits>
#include <typeinfo>

#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

std::string demangle(const char*);

template <typename T>
std::string demangle() {
    return demangle(typeid(T).name());
}

//---------------------------------------------------------------------------------------------------------------------

class Exception : public eckit::Exception {
public:
    using eckit::Exception::Exception;
    ~Exception() override;
};

//---------------------------------------------------------------------------------------------------------------------

class NotEncodable : Exception {
public:
    NotEncodable(const std::string& type_name);

    template <typename T>
    NotEncodable(const T&): NotEncodable{demangle<typename std::decay<T>::type>()} {}

    ~NotEncodable() override;
};

//---------------------------------------------------------------------------------------------------------------------

class NotDecodable : public Exception {
public:
    NotDecodable(const std::string& type_name);

    template <typename T>
    NotDecodable(const T&): NotDecodable{demangle<typename std::decay<T>::type>()} {}

    ~NotDecodable() override;
};

//---------------------------------------------------------------------------------------------------------------------

class InvalidRecord : public Exception {
public:
    InvalidRecord(const std::string& message): Exception("atlas::io::InvalidRecord: " + message) {}

    ~InvalidRecord() override;
};

//---------------------------------------------------------------------------------------------------------------------

class DataCorruption : public Exception {
public:
    DataCorruption(const std::string& message): Exception("atlas::io::DataCorruption: " + message) {}

    ~DataCorruption() override;
};

//---------------------------------------------------------------------------------------------------------------------

class WriteError : public Exception {
public:
    WriteError(const std::string& message): Exception("atlas::io::WriteError: " + message) {}

    ~WriteError() override;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
