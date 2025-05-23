/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#if !defined(HIC_NAMESPACE)
  #define HIC_NAMESPACE
  #define HIC_NAMESPACE_BEGIN
  #define HIC_NAMESPACE_END
#else
  #define HIC_NAMESPACE_BEGIN namespace HIC_NAMESPACE {
  #define HIC_NAMESPACE_END }
#endif