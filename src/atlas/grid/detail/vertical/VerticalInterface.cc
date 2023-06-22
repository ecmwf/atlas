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

#include "VerticalInterface.h"
#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Vertical.h"
#include "atlas/runtime/Exception.h"

namespace atlas {

Vertical* atlas__Vertical__new(idx_t levels, const double z[]) {
    std::vector<double> zvec(z, z + levels);
    return new Vertical(levels, zvec);
}

Vertical* atlas__Vertical__new_interval(idx_t levels, const double z[], const double interval[]) {
    std::vector<double> zvec(z, z + levels);
    return new Vertical(levels, zvec, interval);
}

void atlas__Vertical__delete(Vertical* This) {
    delete This;
}

field::FieldImpl* atlas__Vertical__z(const Vertical* This) {
    ATLAS_ASSERT(This != nullptr);
    field::FieldImpl* field;
    {
        auto zfield = Field("z", array::make_datatype<double>(), array::make_shape(This->size()));
        auto zview  = array::make_view<double, 1>(zfield);
        for (idx_t k = 0; k < zview.size(); ++k) {
            zview(k) = (*This)[k];
        }
        field   = zfield.get();
        field->attach();
    }
    field->detach();
    return field;
}

int atlas__Vertical__size(const Vertical* This) {
    return This->size();
}

}  // namespace atlas
