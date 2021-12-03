/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/TableView.h"

namespace atlas {
namespace array {

// ----------------------------------------------------------------------------

template <bool ReadOnly>
TableView<ReadOnly>::TableView(const Table& table, bool host):
    host_(host),
    missing_value_(table.missing_value()),
    rows_(table.rows()),
    maxcols_(table.maxcols()),
    mincols_(table.mincols()),
    values_(host_ ? make_host_view<idx_t, 1>(*(table.data_[_values_]))
                  : make_device_view<idx_t, 1>(*(table.data_[_values_]))),
    displs_(host_ ? make_host_view<size_t, 1>(*(table.data_[_displs_]))
                  : make_device_view<size_t, 1>(*(table.data_[_displs_]))),
    counts_(host_ ? make_host_view<size_t, 1>(*(table.data_[_counts_]))
                  : make_device_view<size_t, 1>(*(table.data_[_counts_]))),
    const_access_(this),
    access_(this) {}

// ----------------------------------------------------------------------------

template <bool ReadOnly>
TableView<ReadOnly>::TableView(const TableView& other):
    host_(other.host_),
    missing_value_(other.missing_value_),
    rows_(other.rows_),
    maxcols_(other.maxcols_),
    mincols_(other.mincols_),
    values_(other.values_),
    displs_(other.displs_),
    counts_(other.counts_),
    const_access_(this),
    access_(this) {}

// ----------------------------------------------------------------------------

template <bool ReadOnly>
TableView<ReadOnly> TableView<ReadOnly>::operator=(const TableView& other) {
    host_          = other.host_;
    missing_value_ = other.missing_value_;
    rows_          = other.rows_;
    maxcols_       = other.maxcols_;
    mincols_       = other.mincols_;
    values_        = other.values_;
    displs_        = other.displs_;
    counts_        = other.counts_;
    return *this;
}

// ----------------------------------------------------------------------------

template <bool ReadOnly>
TableView<ReadOnly> make_table_view(const Table& table) {
    return TableView<ReadOnly>(table);
}

// ----------------------------------------------------------------------------

template class TableView<true>;
template class TableView<false>;
template TableView<true> make_table_view<true>(const Table&);
template TableView<false> make_table_view<false>(const Table&);

// ----------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
