/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <vector>
#include <string>
#include <string_view>

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
    MultiFieldImpl* multifield;
    {
        MultiField f(*config);
        multifield = f.get();
        multifield->attach();
    }
    multifield->detach();
    ATLAS_ASSERT(multifield);
    return multifield;
}

MultiFieldImpl* atlas__MultiField__create_shape(int kind, int rank, int shapef[], const char* var_names,
        size_t length, size_t size) {
    array::ArrayShape shape;
    shape.resize(rank);
    for (idx_t j = 0, jf = rank - 1; j < rank; ++j) {
        shape[j] = shapef[jf--];
    }

    std::vector<std::string> var_names_str;
    for (size_t jj = 0; jj < size; ++jj) {
        std::string str(std::string_view(var_names + jj * length, length));
        str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
           return !std::isspace(ch);
        }).base(), str.end());
        var_names_str.push_back(str);
    }

    MultiFieldImpl* multifield;
    {
        MultiField f(array::DataType{kind}, shape, var_names_str);
        multifield = f.get();
        multifield->attach();
    }
    multifield->detach();
    ATLAS_ASSERT(multifield);
    return multifield;
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
