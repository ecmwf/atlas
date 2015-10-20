/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_functionspace_EdgeBasedFiniteVolume_h
#define atlas_functionspace_EdgeBasedFiniteVolume_h

#include <string>
#include "atlas/functionspace/Nodes.h"

namespace atlas { class Mesh; }
namespace atlas { class FunctionSpace; }

namespace atlas {
namespace functionspace {

class EdgeBasedFiniteVolume : public Nodes {

public:

  EdgeBasedFiniteVolume(Mesh &, const Halo & = Halo(1) );

  virtual std::string name() const { return "EdgeBasedFiniteVolume"; }

  const Nodes& nodes_fs() const { return *this; }
        Nodes& nodes_fs()       { return *this; }

private: // data

    atlas::FunctionSpace* edges_; // non-const because functionspace may modify mesh
};


// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  EdgeBasedFiniteVolume* atlas__functionspace__EdgeBasedFiniteVolume__new (Mesh* mesh, int halo);
}

} // namespace functionspace
} // namespace atlas


#endif // atlas_numerics_nabla_EdgeBasedFiniteVolume_h
