/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/functionspace.h"
#include "atlas/field.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/output/Gmsh.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/interpolation/Interpolation.h"

#include "atlas/mesh/actions/BuildHalo.h"

using namespace atlas;

//------------------------------------------------------------------------------

class Program : public AtlasTool {
    virtual void execute( const Args& args );
public:
    Program( int argc, char** argv );
  };

//-----------------------------------------------------------------------------

Program::Program( int argc, char** argv ) : AtlasTool( argc, argv ) {
    add_option( new SimpleOption<std::string>( "gridA", "grid A" ) );
    add_option( new SimpleOption<std::string>( "gridB", "grid B" ) );
    add_option( new SimpleOption<bool>( "ghost", "Output ghost elements" ) );
    add_option( new SimpleOption<long>( "haloA", "Halo size" ) );
    add_option( new SimpleOption<long>( "haloB", "Halo size" ) );
}

//-----------------------------------------------------------------------------

void Program::execute( const Args& args ) {

  auto ghost = util::Config("ghost",args.getBool("ghost",false));
  auto haloA = option::halo( args.getLong("haloA",1) );
  auto haloB = option::halo( args.getLong("haloB",1) );

  auto gridA = Grid( args.getString("gridA") );
  auto gridB = Grid( args.getString("gridB") );

  auto meshgenerator = MeshGenerator( "structured" );

  auto distA = grid::Distribution( gridA, grid::Partitioner( "trans" ) );

  auto meshA = meshgenerator.generate( gridA, distA );

  numerics::fvm::Method fvmA(meshA,haloA);
  auto gmshA = output::Gmsh( "meshA.msh", ghost );
  gmshA.write(meshA);

  auto distB = grid::Distribution( gridB, grid::MatchingMeshPartitioner( meshA ) );

  auto meshB = meshgenerator.generate( gridB, distB );

  numerics::fvm::Method fvmB(meshB,haloB);

  Field fieldB = fvmB.node_columns().createField<double>();

  output::Gmsh gmshB( "meshB.msh", ghost );
  gmshB.write(meshB);
  gmshB.write(fieldB);

  Interpolation AtoB( option::type("finite-element"), fvmA.node_columns(), fvmB.node_columns() );
  Interpolation BtoA( option::type("finite-element"), fvmB.node_columns(), fvmA.node_columns() );

}

//------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    Program tool( argc, argv );
    return tool.start();
}
