/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/FieldCreatorArraySpec.h"

#include <algorithm>
#include <sstream>

#include "eckit/config/Parametrisation.h"

#include "atlas/array/DataType.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

FieldImpl* FieldCreatorArraySpec::createField(const eckit::Parametrisation& params) const {
    std::vector<long> shape;
    if (!params.get("shape", shape)) {
        throw_Exception("Could not find parameter 'shape' in Parametrisation");
    }

    std::vector<idx_t> s(shape.size());

    bool fortran(false);
    params.get("fortran", fortran);
    if (fortran) {
        std::reverse_copy(shape.begin(), shape.end(), s.begin());
    }
    else {
        s.assign(shape.begin(), shape.end());
    }

    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (params.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        params.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    int alignment = 1;
    params.get("alignment", alignment);

    std::string name;
    params.get("name", name);
    Log::trace() << "Create field " << name << "\t shape=[";
    for (size_t i = 0; i < s.size(); ++i) {
        Log::trace() << s[i] << (i < s.size() - 1 ? "," : "");
    }
    Log::trace() << "]" << std::endl;

    auto field = FieldImpl::create(name, datatype, array::ArraySpec(std::move(s), array::ArrayAlignment(alignment)), params);
    field->callbackOnDestruction([field]() { Log::trace() << "Destroy field " << field->name() << std::endl; });
    return field;
}

namespace {
static FieldCreatorBuilder<FieldCreatorArraySpec> __ArraySpec("ArraySpec");
}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
