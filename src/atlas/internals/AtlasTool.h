/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_internals_AtlasTool_h
#define atlas_internals_AtlasTool_h

#include <vector>
#include <iostream>

#include "atlas/atlas.h"
#include "atlas/runtime/Behavior.h"
#include "atlas/runtime/Log.h"
#include "eckit/runtime/Tool.h"
#include "eckit/runtime/Context.h"
#include "eckit/mpi/ParallelContextBehavior.h"
#include "eckit/option/CmdArgs.h"
#include "eckit/option/SimpleOption.h"
#include "atlas/internals/Debug.h"

//--------------------------------------------------------------------------------

namespace atlas {

static void usage( const std::string& name )
{
  Log::info() << "dummy usage" << std::endl;
}

class AtlasTool : public eckit::Tool {

protected:

  typedef std::vector<eckit::option::Option *> Options;
  typedef eckit::option::CmdArgs Args;

  virtual std::string indent() { return "      "; }
  virtual std::string briefDescription() { return ""; }
  virtual std::string longDescription() { return ""; }
  virtual std::string usage() { return name() + " [OPTION]... [--help] [--debug]"; }


  void add_option( eckit::option::Option* option )
  {
    options_.push_back( option );
  }

  virtual void help( std::ostream &out = Log::info() )
  {
    out << "NAME\n" << indent() << name();
    std::string brief = briefDescription();
    if ( brief.size() ) out << " - " << brief << '\n';

    std::string usg = usage();
    if ( usg.size() ) {
      out << '\n';
      out << "SYNOPSIS\n" << indent() << usg << '\n';
    }
    std::string desc = longDescription();
    if ( desc.size() ) {
      out << '\n';
      out << "DESCRIPTION\n" << indent() << desc << '\n';
    }
    out << '\n';
    out << "OPTIONS\n";
    for( Options::const_iterator it = options_.begin(); it!= options_.end(); ++it ) {
      out << indent() << (**it) << "\n\n";
    }
    out << std::flush;
  }

  virtual int numberOfPositionalArguments() { return 0; }
  virtual int minimumPositionalArguments() { return 0; }

  virtual void run()
  {
    if( eckit::Context::instance().argc()-1 < minimumPositionalArguments() )
    {
      Log::info() << "Usage: " << usage() << std::endl;
      return;
    }
    Options opts = options_;
    Args args(&atlas::usage,
        opts,
        numberOfPositionalArguments(),
        minimumPositionalArguments());
    bool _help(false);
    args.get("help",_help);
    if( _help )
    {
      if( eckit::mpi::rank() == 0 )
        help();
      return;
    }

    long debug(0);
    const char* env_debug = ::getenv("DEBUG");
    if( env_debug ) debug = ::atol(env_debug);
    args.get("debug",debug);

    eckit::Context::instance().behavior( new atlas::runtime::Behavior() );
    eckit::Context::instance().debug(debug);

    atlas_init();
    execute(args);
    atlas_finalize();
  }

  virtual void execute(const Args&) = 0;

public:

  AtlasTool(int argc,char **argv): eckit::Tool(argc,argv)
  {
    eckit::Context::instance().behavior( new eckit::mpi::ParallelContextBehavior() );
    add_option( new eckit::option::SimpleOption<bool>("help","Print this help") );
    add_option( new eckit::option::SimpleOption<long>("debug","Debug level") );
  }

private:
  Options options_;
};

} // namespace atlas

//----------------------------------------------------------------------------------

#endif
