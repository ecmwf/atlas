/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Healpix : public Structured {
public:
    using Structured::Structured;
    Healpix(long N, const std::string& ordering = "ring");
    Config meshgenerator() const override;
    Config partitioner() const override;
    static std::string static_type() { return "healpix"; }
    virtual std::string type() const override { return static_type(); }
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
