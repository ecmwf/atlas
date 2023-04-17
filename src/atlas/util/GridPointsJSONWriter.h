/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <iosfwd>

#include "eckit/config/Parametrisation.h"

#include "atlas/grid.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

class GridPointsJSONWriter {
public:
    GridPointsJSONWriter(Grid grid, const eckit::Parametrisation& args);

    void write(std::ostream& out, eckit::Channel& info) const;

    void write(std::ostream& out, std::ostream* info = nullptr) const;

private:
    Grid grid_;
    grid::Distribution distribution_;
    int precision_;
    int verbose_;
    int nb_partitions_;
    int partition_;
    bool pretty_;
    std::string field_;
    std::vector<long> points_;
    long points_base_;
    long field_base_;
};

//------------------------------------------------------------------------------

} // namespace util
} // namespace atlas