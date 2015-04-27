
debug_var(ATLAS_LIBRARIES)

ecbuild_add_executable(
  TARGET      atlas-sandbox-trans
  SOURCES     main.cc
  LIBS        ${ATLAS_LIBRARIES} ${TRANSI_LIBRARIES}
  INCLUDES    ${TRANSI_INCLUDE_DIRS}
  CONDITION   ATLAS_HAVE_TRANS )

ecbuild_add_executable(
  TARGET      atlas-sandbox-interface
  SOURCES     interface.c
  LIBS        ${ATLAS_LIBRARIES} ${TRANSI_LIBRARIES}
  INCLUDES    ${TRANSI_INCLUDE_DIRS}
  CONDITION   ATLAS_HAVE_TRANS )