#ifndef atlas_fortran_eckit_logging_h
#define atlas_fortran_eckit_logging_h

extern "C"
{
  void eckit__log_debug_set_level (int level);
	void eckit__log_debug (int level, char* name);
	void eckit__log_info (char* name);
	void eckit__log_warning (char* name);
	void eckit__log_error (char* name);
	void eckit__log_panic (char* name);
	void eckit__log_debug_endl (int level, char* name);
	void eckit__log_info_endl (char* name);
	void eckit__log_warning_endl (char* name);
	void eckit__log_error_endl (char* name);
	void eckit__log_panic_endl (char* name);
}

#endif
