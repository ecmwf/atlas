#ifndef atlas_fortran_atlas_logging_h
#define atlas_fortran_atlas_logging_h

#include "eckit/log/MultiChannel.h"
typedef eckit::MultiChannel AtlasChannel;

extern "C"
{
  void atlas__log_debug_set_level (int level);
  void atlas__log_debug (int lvl, char* msg, int endl, int flush);
  void atlas__log (int cat, int lvl, char* msg, int endl, int flush);
  void atlas__cat__connect_fortran_unit (int cat, int unit);
  void atlas__cat__connect_stdout (int cat);
  void atlas__cat__connect_stderr (int cat);
  void atlas__cat__disconnect_stdout (int cat);
  void atlas__cat__disconnect_stderr (int cat);
  void atlas__cat__disconnect_fortran_unit (int cat, int unit);
  void atlas__cat__set_prefix_stdout (int cat, char* prefix);
  void atlas__cat__set_prefix_stderr (int cat, char* prefix);
  void atlas__cat__set_prefix_fortran_unit (int cat, int unit, char* prefix);
  AtlasChannel* atlas__log_channel (int cat);
  void atlas__Channel__connect_fortran_unit (AtlasChannel* ch, int unit);
  void atlas__Channel__connect_stdout (AtlasChannel* ch);
  void atlas__Channel__connect_stderr (AtlasChannel* ch);
  void atlas__Channel__disconnect_fortran_unit (AtlasChannel* ch, int unit);
  void atlas__Channel__disconnect_stdout (AtlasChannel* ch);
  void atlas__Channel__disconnect_stderr (AtlasChannel* ch);
  void atlas__Channel__set_prefix_stdout (AtlasChannel* ch, char* prefix);
  void atlas__Channel__set_prefix_stderr (AtlasChannel* ch, char* prefix);
  void atlas__Channel__set_prefix_fortran_unit (AtlasChannel* ch, int unit, char* prefix);
}

#endif
