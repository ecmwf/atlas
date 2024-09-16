/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Slavko Brdar
/// @author Willem Deconinck
/// @date September 2024

#include <cmath>
#include <sstream>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/DataType.h"
#include "atlas/array/Range.h"
#include "atlas/field/MultiFieldCreatorArray.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

MultiFieldCreatorArray::MultiFieldCreatorArray() {}

MultiFieldCreatorArray::MultiFieldCreatorArray(const eckit::Configuration&) {}

MultiFieldCreatorArray::~MultiFieldCreatorArray() = default;

MultiFieldImpl* MultiFieldCreatorArray::create(const eckit::Configuration& config) const {
    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (config.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        config.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }
    std::vector<int> shape;
    config.get("shape", shape);
    const auto fields = config.getSubConfigurations("fields");
    int nflds = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        long nvar = 1;
        fields[i].get("nvar", nvar);
        nflds += nvar;
    }
    std::vector<std::string> var_names;
    var_names.resize(nflds);
    for (size_t i = 0, cnt = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        long nvar = 1;
        fields[i].get("nvar", nvar);
        if (nvar > 1) {
            for (int ivar = 0; ivar < nvar; ivar++) {
                std::stringstream ss;
                ss << name << "_" << ivar;
                var_names[cnt++] = ss.str();
            }

        }
        else {
                var_names[cnt++] = name;
        }
    }
    return create(datatype, shape, var_names);
}

MultiFieldImpl* MultiFieldCreatorArray::create(const array::DataType datatype, const std::vector<int>& shape, const std::vector<std::string>& var_names) const {
    const int dim = shape.size();
    const int nvar = var_names.size();
    ATLAS_ASSERT(nvar > 0 && dim > 2, "MultiField must have at least one field name.");

    int varidx = -1;
    for (int i = 0; i < dim; i++) {
        if (shape[i] == -1) {
            varidx = i;
            break;
        }
    }

    array::ArrayShape multiarray_shape = shape;
    multiarray_shape[varidx] = nvar;


    MultiFieldImpl* multifield = new MultiFieldImpl{array::ArraySpec{datatype, multiarray_shape}};
    auto& multiarray = multifield->array();
      
    array::ArrayShape field_shape;
    field_shape.resize(multiarray_shape.size() - 1);
    array::ArrayStrides field_strides;
    field_strides.resize(multiarray_shape.size() - 1);
    for (int ivar = 0, i = 0; ivar < nvar; ivar++) {
        if (ivar != varidx) {
            field_shape[i]   = multiarray.shape()[ivar];
            field_strides[i] = multiarray.strides()[ivar];
            ++i;
        } 
    }
    array::ArraySpec field_array_spec(field_shape, field_strides);

    for (size_t ifield = 0; ifield < nvar; ++ifield) {
        idx_t start_index = multiarray.strides()[varidx] * ifield;
        Field field;
        if (datatype.kind() == array::DataType::KIND_REAL64) {
            double* slice_begin = multiarray.data<double>() + start_index;
            field = Field(var_names[ifield], slice_begin, field_array_spec);
        }
        else if (datatype.kind() == array::DataType::KIND_REAL32) {
            float* slice_begin = multiarray.data<float>() + start_index;
            field      = Field(var_names[ifield], slice_begin, field_array_spec);
        }
        else if (datatype.kind() == array::DataType::KIND_INT32) {
            int* slice_begin = multiarray.data<int>() + start_index;
            field      = Field(var_names[ifield], slice_begin, field_array_spec);
        }
        else if (datatype.kind() == array::DataType::KIND_INT64) {
            long* slice_begin = multiarray.data<long>() + start_index;
            field      = Field(var_names[ifield], slice_begin, field_array_spec);
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
        multifield->add(field);
    }
    Log::debug() << "Creating multifield of " << datatype.str() << " type" << std::endl;
    return multifield;
}

// ------------------------------------------------------------------

namespace {  
static MultiFieldCreatorBuilder<MultiFieldCreatorArray> __MultiFieldCreatorArray("MultiFieldCreatorArray");
}

}  // namespace field
}  // namespace atlas
