
debug_var(ATLAS_LIBRARIES)

ecbuild_add_executable(
  TARGET      atlas-sandbox-trans
  SOURCES     main.cc
  LIBS        ${ATLAS_LIBRARIES}
  CONDITION   HAVE_TRANS )
