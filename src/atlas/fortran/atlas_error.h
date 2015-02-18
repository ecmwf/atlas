#ifndef atlas_fortran_atlas_error_h
#define atlas_fortran_atlas_error_h

#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace fortran {

class Error
{
private: 
  Error();
  
public:
  static Error& instance();
    
  int code() const { return code_; }
  
  const std::string& msg() const { return msg_; }
  
  bool throws() const { return throws_; }
  
  bool backtrace() const { return backtrace_; }

  bool aborts() const { return aborts_; }
  
  void set_code( int c ) { code_ = c; }

  void set_msg( const std::string& m ) { msg_ = m; }

  void set_aborts( bool b ) { aborts_ = b; }
  
  void set_throws( bool b ) { throws_ = b; }

  void set_backtrace( bool b ) { backtrace_ = b; }
  
  void clear();
  
private:
  std::string msg_;
  int code_;
  bool aborts_;
  bool throws_;
  bool backtrace_;
};

}
}

extern "C"
{
  void atlas__abort (char* msg, char* file, int line, char* function);
  void atlas__throw_exception (char* msg, char* file, int line, char* function);
  void atlas__throw_notimplemented (char* msg, char* file, int line, char* function);
  int atlas__Error_code ();
  void atlas__Error_clear ();
  void atlas__Error_success ();
  void atlas__Error_set_aborts (int on_off);
  void atlas__Error_set_throws (int on_off);
  void atlas__Error_set_backtrace (int on_off);
  char* atlas__Error_msg ();
}

#endif
