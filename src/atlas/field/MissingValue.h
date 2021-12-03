/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <string>

#include "atlas/field/detail/MissingValue.h"
#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"


namespace atlas {
class Field;
namespace util {
class Config;
}
}  // namespace atlas


namespace atlas {
namespace field {


/**
 * @brief Missing values indicator
 */
struct MissingValue : DOXYGEN_HIDE(public util::ObjectHandle<detail::MissingValue>) {
    using Spec   = util::Config;
    using Config = detail::MissingValue::Config;

    /**
     * @brief ctor
     */
    using Handle::Handle;
    MissingValue();
    MissingValue(const Config&);
    MissingValue(const std::string& type, const Config&);
    MissingValue(const Field&);
    MissingValue(const std::string& type, const Field&);

    /**
     * @brief bool operator
     * @return if MissingValue object has been setup
     */
    using Handle::operator bool;  // (ensure this exists)

    /**
     * @brief operator() on user-defined values
     * @param [in] value user-defined value
     * @return if user-defined value indicates a missing value
     */
    bool operator()(const double& value) const;

    /**
     * @brief operator() on user-defined values
     * @param [in] value user-defined value
     * @return if user-defined value indicates a missing value
     */
    bool operator()(const float& value) const;

    /**
     * @brief operator() on user-defined values
     * @param [in] value user-defined value
     * @return if user-defined value indicates a missing value
     */
    bool operator()(const int& value) const;

    /**
     * @brief operator() on user-defined values
     * @param [in] value user-defined value
     * @return if user-defined value indicates a missing value
     */
    bool operator()(const long& value) const;

    /**
     * @brief operator() on user-defined values
     * @param [in] value user-defined value
     * @return if user-defined value indicates a missing value
     */
    bool operator()(const unsigned long& value) const;

    /**
     * @brief if missing value is represented by not-a-number
     * @return if missing value is represented by not-a-number
     */
    bool isnan() const;

    /**
     * @brief reference to internal implementation
     * @return reference to internal implementation
     */
    const detail::MissingValue& ref() const;

    /**
     * @brief fill missing value metadata on Field
     * @param [in] field field to set
     */
    void metadata(Field& field) const;
};


}  // namespace field
}  // namespace atlas
