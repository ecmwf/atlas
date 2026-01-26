/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/linalg/sparse.h"
#include "atlas/util/Config.h"

namespace atlas {

class ScripIO {

public:
    explicit ScripIO(const util::Config& = util::NoConfig()) {}

    static atlas::linalg::SparseMatrixStorage read(const std::string&);
    static void write(const linalg::SparseMatrixStorage&, const std::string&);
};

}
