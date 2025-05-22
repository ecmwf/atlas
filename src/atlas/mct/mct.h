/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iterator>
#include <vector>
#include <map>

#include "atlas/grid.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mct/ParInter.h"
#include "atlas/mct/coupler.h"
#include "atlas/parallel/mpi/mpi.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;
using Index = eckit::linalg::Index;
using Value = eckit::linalg::Scalar;

namespace atlas::mct {

    void set_default_comm_to_local(int model_id);

    void setup_oneway_remap(int model_1, int model_2);

    void setup_coupler(int model_id, Grid grid);

    void finalise_coupler();

    void put_field(Field f, int model_2, int tstep = 0);

    void get_field(Field f, int model_1, int tstep = 0);

} // namespace atlas::mct
