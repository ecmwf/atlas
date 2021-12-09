/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "SerialPartitioner.h"

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::SerialPartitioner> __Serial(
    atlas::grid::detail::partitioner::SerialPartitioner::static_type());
}

atlas::grid::detail::partitioner::SerialPartitioner::SerialPartitioner(): Partitioner(1) {}
