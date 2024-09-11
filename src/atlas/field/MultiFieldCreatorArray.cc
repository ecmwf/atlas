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

MultiFieldImpl* MultiFieldCreatorArray::create(const eckit::Configuration&) const {
    ATLAS_NOTIMPLEMENTED;
    return nullptr;
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
      
    std::vector<int> field_shape_vec(multiarray_shape.size() - 1);
    std::vector<int> field_strides_vec(multiarray_shape.size() - 1);
    for (int ivar = 0, i = 0; ivar < nvar; ivar++) {
        if(ivar != varidx) {
            field_shape_vec[i]   = multiarray.shape()[ivar];
            field_strides_vec[i] = multiarray.strides()[ivar];
            ++i;
        } 
    }
    array::ArrayShape field_shape(field_shape_vec);
    array::ArrayStrides field_strides(field_strides_vec);
    array::ArraySpec field_array_spec = array::ArraySpec(field_shape, field_strides);

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
