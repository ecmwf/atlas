#ifndef atlas_fortran_atlas_logging_h
#define atlas_fortran_atlas_logging_h

#include "eckit/log/Channel.h"
typedef eckit::Channel AtlasChannel;

extern "C"
{
  void atlas__log_set_debug (int level);
  void atlas__log_debug (int lvl, char* msg, int endl, int flush);
  void atlas__log_cat(int cat, int lvl, char* msg, int endl, int flush);
  void atlas__logcat__connect_fortran_unit (int cat, int unit);
  void atlas__logcat__connect_stdout (int cat);
  void atlas__logcat__connect_stderr (int cat);
  void atlas__logcat__disconnect_stdout (int cat);
  void atlas__logcat__disconnect_stderr (int cat);
  void atlas__logcat__disconnect_fortran_unit (int cat, int unit);
  void atlas__logcat__set_prefix_stdout (int cat, char* prefix);
  void atlas__logcat__set_prefix_stderr (int cat, char* prefix);
  void atlas__logcat__set_prefix_fortran_unit (int cat, int unit, char* prefix);
  void atlas__logcat__indent (int cat, char* indent);
  void atlas__logcat__dedent (int cat);
  void atlas__logcat__clear_indentation (int cat);
  AtlasChannel* atlas__LogChannel_cat (int cat);
  void atlas__LogChannel__connect_fortran_unit (AtlasChannel* ch, int unit);
  void atlas__LogChannel__connect_stdout (AtlasChannel* ch);
  void atlas__LogChannel__connect_stderr (AtlasChannel* ch);
  void atlas__LogChannel__disconnect_fortran_unit (AtlasChannel* ch, int unit);
  void atlas__LogChannel__disconnect_stdout (AtlasChannel* ch);
  void atlas__LogChannel__disconnect_stderr (AtlasChannel* ch);
  void atlas__LogChannel__set_prefix (AtlasChannel* ch, char* prefix);
  void atlas__LogChannel__set_prefix_stdout (AtlasChannel* ch, char* prefix);
  void atlas__LogChannel__set_prefix_stderr (AtlasChannel* ch, char* prefix);
  void atlas__LogChannel__set_prefix_fortran_unit (AtlasChannel* ch, int unit, char* prefix);
  void atlas__LogChannel__indent (AtlasChannel* ch, char* indent);
  void atlas__LogChannel__indent_stdout (AtlasChannel* ch, char* indent);
  void atlas__LogChannel__indent_stderr (AtlasChannel* ch, char* indent);
  void atlas__LogChannel__indent_fortran_unit (AtlasChannel* ch, int unit, char* indent);
  void atlas__LogChannel__dedent (AtlasChannel* ch);
  void atlas__LogChannel__dedent_stdout (AtlasChannel* ch);
  void atlas__LogChannel__dedent_stderr (AtlasChannel* ch);
  void atlas__LogChannel__dedent_fortran_unit (AtlasChannel* ch, int unit);
  void atlas__LogChannel__clear_indentation (AtlasChannel* ch);
  void atlas__LogChannel__clear_indentation_stdout (AtlasChannel* ch);
  void atlas__LogChannel__clear_indentation_stderr (AtlasChannel* ch);
  void atlas__LogChannel__clear_indentation_fortran_unit (AtlasChannel* ch, int unit);

}

#endif
