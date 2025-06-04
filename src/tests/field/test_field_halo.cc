/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "tests/AtlasTestEnvironment.h"

#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/mesh.h"
#include "atlas/functionspace.h"
#include "atlas/field.h"
#include "atlas/option.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test ghost at end") {
  auto g = Grid("O32");
  auto mg = MeshGenerator("structured", util::Config("ghost_at_end", true));
  auto m = mg.generate(g);
  auto fs = functionspace::NodeColumns(m, option::halo(1));
  auto f = fs.createField<double>();
  auto h = f.halo();
  EXPECT(h.contiguous());
  EXPECT(h.appended());

  if (mpi::size() == 1) {
    EXPECT_EQ(h.size(), 192);
    EXPECT_EQ(h.begin(), 5248);
    EXPECT_EQ(h.end(), 5440);
  }
  if (mpi::size() == 4) {
    if (mpi::rank() == 0) {
      EXPECT_EQ(h.size(), 284);
      EXPECT_EQ(h.begin(), 1312);
      EXPECT_EQ(h.end(), 1596);
    }
  }

  EXPECT(not h.updated());
  h.updated(true);
  EXPECT(h.updated());
}

//-----------------------------------------------------------------------------

CASE("test ghost interleaved") {
  auto g = Grid("O32");
  auto mg = MeshGenerator("structured", util::Config("ghost_at_end", false));
  auto m = mg.generate(g);
  auto fs = functionspace::NodeColumns(m, option::halo(1));
  auto f = fs.createField<double>();

  auto h = f.halo();
  EXPECT(not h.contiguous());
  EXPECT(not h.appended());

  if (mpi::size() == 1) {
    EXPECT_EQ(h.size(), 192);
    EXPECT_EQ(h.begin(), 20);
    EXPECT_EQ(h.end(), 5440);
  }
  if (mpi::size() == 4) {
    if (mpi::rank() == 0) {
      EXPECT_EQ(h.size(), 284);
      EXPECT_EQ(h.begin(), 20);
      EXPECT_EQ(h.end(), 1596);
    }
  }

  EXPECT(not h.updated());
  h.updated(true);
  EXPECT(h.updated());
}

//-----------------------------------------------------------------------------

CASE("test halo exchange") {
  auto fs = functionspace::NodeColumns{Mesh{Grid{"O32"}}, option::halo(1)};
  auto f = fs.createField<double>();
  EXPECT(not f.halo().updated());
  f.halo().update();
  EXPECT(f.halo().updated());
  f.halo().invalidate();
  EXPECT(not f.halo().updated());
}

//-----------------------------------------------------------------------------

CASE("test halo exchange on device") {
  auto fs = functionspace::NodeColumns{Mesh{Grid{"O32"}}, option::halo(1)};
  auto f = fs.createField<double>();
  EXPECT(not f.halo().updated());
  f.updateDevice();
  f.halo().update(option::on_device());
  EXPECT(f.halo().updated());
  f.halo().invalidate();
  EXPECT(not f.halo().updated());
}

//-----------------------------------------------------------------------------


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
