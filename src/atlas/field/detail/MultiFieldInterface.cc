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
#include "atlas/field/detail/MultiFieldInterface.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

extern "C" {
MultiFieldImpl* atlas__MultiField__create(eckit::Configuration* config) {
    ATLAS_ASSERT(config != nullptr);

    // Register in factory
    // TODO: move __multiFieldCreatorIFS out of test_multifield_ifs.cc
    //MultiFieldCreatorBuilder<MultiFieldCreatorIFS> __MultiFieldCreatorIFS("MultiFieldCreatorIFS");

    //MultiField* multifield = new MultiField(*config);
    //return multifield->get();

    // TODO
    // repeat here the code for __multiFieldCreatorIFS out of test_multifield_ifs.cc
    long nproma = config->getLong("nproma");
    long nlev = config->getLong("nlev");
    long nblk = 0;
    if (config->has("nblk")) {
        nblk = config->getLong("nblk");
    }
    else if (config->has("ngptot")) {
        long ngptot = config->getLong("ngptot");
        nblk        = std::ceil(static_cast<double>(ngptot) / static_cast<double>(nproma));
    }
    else {
        ATLAS_THROW_EXCEPTION("Configuration not found: ngptot or nblk");
    }
    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (config->get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        config->get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    auto fields = config->getSubConfigurations("fields");
    long nfld   = 0;
    for (const auto& field_config : fields) {
        long nvar = 1;
        field_config.get("nvar", nvar);
        nfld += nvar;
    }
    auto multiarray_shape = array::make_shape(nblk, nfld, nlev, nproma);

MultiFieldImpl* multifield = new MultiFieldImpl{array::ArraySpec{datatype, multiarray_shape}};
    auto& multiarray = multifield->array();

    size_t multiarray_field_idx = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        Field field;
        size_t field_vars = 1;

        if (fields[i].get("nvar", field_vars)) {
            auto field_shape =
                array::make_shape(multiarray.shape(0), field_vars, multiarray.shape(2), multiarray.shape(3));
            auto field_strides    = multiarray.strides();
            auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            constexpr auto all = array::Range::all();
            const auto range   = array::Range(multiarray_field_idx, multiarray_field_idx + field_vars);
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 4>(multiarray).slice(all, range, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 4>(multiarray).slice(all, range, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            field.set_variables(field_vars);
        }
        else {
            auto field_shape   = array::make_shape(multiarray.shape(0), multiarray.shape(2), multiarray.shape(3));
            auto field_strides = array::make_strides(multiarray.stride(0), multiarray.stride(2), multiarray.stride(3));
            auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            constexpr auto all = array::Range::all();
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        field.set_levels(nlev);
        // field.set_blocks(nblk);

        multifield->add(field);

        multiarray_field_idx += field_vars;
    }
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
