/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/**
 * @file relayout_on_host_double.cc
 * @brief Explicit template instantiations of host relayout functions for double.
 *
 * All template implementations live in relayout_on_host.tcc.
 */

#include "relayout_on_host.tcc"

#define EXPLICIT_TEMPLATE_INSTANTIATION(RANK) \
    ATLAS_RELAYOUT_EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(double, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION(2)
EXPLICIT_TEMPLATE_INSTANTIATION(3)
EXPLICIT_TEMPLATE_INSTANTIATION(4)

#undef EXPLICIT_TEMPLATE_INSTANTIATION
