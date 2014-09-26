#ifndef atlas_fortran_atlas_logging_h
#define atlas_fortran_atlas_logging_h
#include <stdio.h>
extern "C"
{
	void atlas__log_debug_set_level (int level);
	void atlas__log_debug (int lvl, char* msg, int endl, int flush);
	void atlas__log(int cat, int lvl, char* msg, int endl, int flush);
}

#endif
