#include "eckit/log/CodeLocation.h"
#include "eckit/os/BackTrace.h"
#include "atlas/runtime/ErrorHandling.h"
#include "eckit/config/Resource.h"
#include "atlas/mpi/mpi.h"
#include "atlas/runtime/Log.h"

using namespace atlas;

namespace atlas {

Error::Error()
{
  clear();
  throws_    = eckit::Resource<bool>("atlas.error.throws;$ATLAS_ERROR_THROWS",false);
  aborts_    = eckit::Resource<bool>("atlas.error.aborts;$ATLAS_ERROR_ABORTS",true);
  backtrace_ = eckit::Resource<bool>("atlas.error.backtrace;$ATLAS_ERROR_BACKTRACE",true);
}

Error& Error::instance()
{
  static Error error_instance;
  return error_instance;
}

void Error::clear()
{
  code_      = atlas_err_cleared;
  msg_       = std::string("Error code was not set!");
}

void handle_error(const eckit::Exception& exception, const int err_code)
{
  std::stringstream msg;
  if( Error::instance().backtrace() || Error::instance().aborts() )
  {
    msg << "=========================================\n"
        << "ERROR\n"
        << "-----------------------------------------\n"
        << exception.what() << "\n";
    if( exception.location() )
      msg << "-----------------------------------------\n"
          << "LOCATION: " << exception.location() << "\n";


    msg << "-----------------------------------------\n"
        << "BACKTRACE\n"
        << "-----------------------------------------\n"
        << exception.callStack() << "\n"
        << "=========================================";
  }
  else
  {
    msg << exception.what();
  }
  Error::instance().set_code(err_code);
  Error::instance().set_msg(msg.str());

  if( Error::instance().aborts() )
  {
    Log::error() << msg.str() << std::endl;
    MPI_Abort(eckit::mpi::comm(), err_code);
  }
  if( Error::instance().throws() )
  {
    throw exception;
  }

}

}

using eckit::CodeLocation;
using eckit::Exception;
using eckit::BackTrace;
using eckit::Exception;
using eckit::NotImplemented;
using eckit::OutOfRange;
using eckit::UserError;
using eckit::AssertionFailed;
using eckit::SeriousBug;


void atlas__Error_set_aborts (int on_off)
{
  Error::instance().set_aborts(on_off);
}

void atlas__Error_set_throws (int on_off)
{
  Error::instance().set_throws(on_off);
}

void atlas__Error_set_backtrace (int on_off)
{
  Error::instance().set_backtrace(on_off);
}

void atlas__Error_success ()
{
  Error::instance().set_code(atlas_err_noerr);
  Error::instance().set_msg(std::string());
}

void atlas__Error_clear ()
{
  Error::instance().clear();
}

int atlas__Error_code()
{
  return Error::instance().code();
}

char* atlas__Error_msg()
{
  return const_cast<char*>(Error::instance().msg().c_str());
}

template< typename EXCEPTION>
EXCEPTION create_exception(char* msg, char* file, int line, char* function)
{
  if( file && std::string(file).size() && std::string(msg).size() )
    return EXCEPTION( std::string(msg), CodeLocation(file,line,function) );
  else if( file && std::string(file).size() )
    return EXCEPTION( std::string(), CodeLocation(file,line,function) );
  else if( std::string(msg).size() )
    return EXCEPTION( std::string(msg), CodeLocation() );
  else
    return EXCEPTION( std::string(), CodeLocation() );
}

void atlas__throw_exception(char* msg, char* file, int line, char* function)
{
  Exception exception ( create_exception<Exception>(msg,file,line,function) );
  handle_error(exception,atlas_err_exception);
}

void atlas__throw_notimplemented (char* msg, char* file, int line, char* function)
{
  NotImplemented exception ( create_exception<NotImplemented>(msg,file,line,function) );
  handle_error(exception,atlas_err_notimplemented);
}

void atlas__throw_outofrange (char* msg, char* file, int line, char* function)
{
  OutOfRange exception ( create_exception<OutOfRange>(msg,file,line,function) );
  handle_error(exception,atlas_err_outofrange);
}

void atlas__throw_usererror (char* msg, char* file, int line, char* function)
{
  UserError exception ( create_exception<UserError>(msg,file,line,function) );
  handle_error(exception,atlas_err_usererror);
}

void atlas__throw_assertionfailed (char* msg, char* file, int line, char* function)
{
  AssertionFailed exception ( create_exception<AssertionFailed>(msg,file,line,function) );
  handle_error(exception,atlas_err_assertionfailed);
}

void atlas__throw_seriousbug (char* msg, char* file, int line, char* function)
{
  SeriousBug exception ( create_exception<SeriousBug>(msg,file,line,function) );
  handle_error(exception,atlas_err_seriousbug);
}

void atlas__abort(char* msg, char* file, int line, char* function )
{
  Log::error() << "=========================================\n"
               << "ABORT\n";
  if( msg && std::string(msg).size() )
    Log::error() << "-----------------------------------------\n"
                 << msg << "\n";

  if( file && std::string(file).size() )
    Log::error() << "-----------------------------------------\n"
                 << "LOCATION: " << CodeLocation(file,line,function) << "\n";

  Log::error() << "-----------------------------------------\n"
               << "BACKTRACE\n"
               << "-----------------------------------------\n"
               << BackTrace::dump() << "\n"
               << "========================================="
               << std::endl;
  MPI_Abort(eckit::mpi::comm(), -1);
}

void atlas__error_example()
{
  ATLAS_ERROR_HANDLING(
    throw OutOfRange(12,5,Here())
  );
}

