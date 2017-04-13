/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include <iostream>

#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"
#include "eckit/runtime/Tool.h"
#include "eckit/option/CmdArgs.h"
#include "eckit/option/SimpleOption.h"
#include "eckit/option/VectorOption.h"
#include "eckit/option/Separator.h"
#include "atlas/parallel/mpi/mpi.h"

//--------------------------------------------------------------------------------

using eckit::option::SimpleOption;
using eckit::option::VectorOption;
using eckit::option::Separator;
using eckit::option::CmdArgs;
using eckit::option::Option;

namespace atlas {

static void usage( const std::string& name )
{
  Log::info() << "dummy usage" << std::endl;
}

class AtlasTool : public eckit::Tool {


protected:

  typedef std::vector<eckit::option::Option *> Options;
  typedef eckit::option::CmdArgs Args;

  virtual bool serial() { return false; }
  virtual std::string indent() { return "      "; }
  virtual std::string briefDescription() { return ""; }
  virtual std::string longDescription() { return ""; }
  virtual std::string usage() { return name() + " [OPTION]... [--help,-h] [--debug]"; }


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

  virtual int numberOfPositionalArguments() { return -1; }
  virtual int minimumPositionalArguments() { return 0; }

  bool handle_help()
  {
    for( int i=1; i<argc(); ++i )
    {
      if( argv(i) == "--help" ||
          argv(i) == "-h"     )
      {
        if( parallel::mpi::comm().rank() == 0 )
          help();
        return true;
      }
    }
    return false;
  }

  virtual void run()
  {
    if( handle_help() )
      return;

    if( argc()-1 < minimumPositionalArguments() )
    {
      Log::info() << "Usage: " << usage() << std::endl;
      return;
    }
    Options opts = options_;
    Args args(&atlas::usage,
        opts,
        numberOfPositionalArguments(),
        minimumPositionalArguments());

    atlas::Library::instance().initialise();
    execute(args);
    atlas::Library::instance().finalise();
  }

  virtual void execute(const Args&) = 0;

public:

  AtlasTool(int argc,char **argv): eckit::Tool(argc,argv)
  {
    add_option( new SimpleOption<bool>("help","Print this help") );
    add_option( new SimpleOption<long>("debug","Debug level") );
    taskID( eckit::mpi::comm("world").rank());
    if( taskID() != 0 )
        Log::reset();
  }

private:
  Options options_;
};

} // namespace atlas
