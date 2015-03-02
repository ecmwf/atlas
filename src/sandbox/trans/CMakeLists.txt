
debug_var(ATLAS_LIBRARIES)

ecbuild_add_executable(
  TARGET      atlas-sandbox-trans
  SOURCES     main.cc
  INCLUDES    ${TRANS_API_INCLUDE_DIRS}
  LIBS        ${TRANS_API_LIBRARIES}
  CONDITION   HAVE_TRANS )
