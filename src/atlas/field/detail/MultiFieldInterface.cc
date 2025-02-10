/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <cstring>
#include <sstream>

#include "atlas/library/config.h"
#include "atlas/field/MultiField.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/field/detail/MultiFieldInterface.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

extern "C" {
MultiFieldImpl* atlas__MultiField__create(eckit::Configuration* config) {
    ATLAS_ASSERT(config != nullptr);
    auto multifield = new MultiField(*config);
    ATLAS_ASSERT(multifield);
    return multifield->get();
}

MultiFieldImpl* atlas__MultiField__create_shape(int kind, int rank, int shapef[], const char* var_names,
        size_t length, size_t size) {
    array::ArrayShape shape;
    shape.resize(rank);
    array::ArrayStrides strides;
    for (idx_t j = 0, jf = rank - 1; j < rank; ++j) {
        shape[j]   = shapef[jf--];
    }

    std::vector<std::string> var_names_str;
    for (size_t jj = 0; jj < size; ++jj) {
        std::unique_ptr<char[]> str(new char[length + 1]);
        ATLAS_ASSERT(snprintf(str.get(), length, "%s", var_names + jj * length ) >= 0);
        std::string sstr(str.get());
        sstr.erase(std::find_if(sstr.rbegin(), sstr.rend(), [](unsigned char ch) {
           return !std::isspace(ch);
        }).base(), sstr.end());
        var_names_str.push_back(sstr);
    }

    auto multifield = new MultiField(array::DataType{kind}, shape, var_names_str);
    ATLAS_ASSERT(multifield);
    return multifield->get();
}

void atlas__MultiField__delete(MultiFieldImpl* This) {
    delete This;
}

int atlas__MultiField__size(MultiFieldImpl* This) {
    return This->size();
}

FieldSetImpl* atlas__MultiField__fieldset(MultiFieldImpl* This) {
    return This->fieldset().get();
}

}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
