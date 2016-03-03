
#include "eckit/config/FileReadPolicy.h"
#include "eckit/runtime/Context.h"
#include "eckit/runtime/ContextBehavior.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas_f/internals/atlas_read_file.h"

extern "C"
{
  
//-----------------------------------------------------------------------------

int atlas__read_file( const char* path, char* &content, int& size )
{
  
  ATLAS_ERROR_HANDLING(
    
    eckit::FileReadPolicy p = eckit::Context::instance().behavior().fileReadPolicy();

    std::stringstream ss;

    if( read( p, path, ss ) )
    {
      std::string s = ss.str();
      size = s.size()+1;
      content = new char[size];
      strcpy(content,s.c_str());
      return true;
    }
  );
  return false;
}

//-----------------------------------------------------------------------------

}

