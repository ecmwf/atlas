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

#include "eckit/config/Parametrisation.h"
#include "atlas/functionspace/HaloDescription.h"

namespace atlas::field {

class FieldImpl;

class Halo: public functionspace::HaloDescription {
public:
  bool updated() const;
  void updated(bool v);
  void update();
  void update(const eckit::Parametrisation&);
  void invalidate() { updated(false); }
  bool appended() const;

private:
  friend class FieldImpl;
  Halo(field::FieldImpl& f);
private:
  field::FieldImpl& field_;
};

} // namespace atlas::field
