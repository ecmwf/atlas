
debug_var(ATLAS_LIBRARIES)

ecbuild_add_executable(
  TARGET      atlas-sandbox-trans
  SOURCES     main.cc
<<<<<<< Updated upstream
  INCLUDES    ${TRANS_API_INCLUDE_DIRS}
  LIBS        ${TRANS_API_LIBRARIES}
=======
<<<<<<< Updated upstream
  LIBS        ${ATLAS_LIBRARIES}
>>>>>>> Stashed changes
  CONDITION   HAVE_TRANS )
=======
  INCLUDES    ${TRANSI_INCLUDE_DIRS}
  LIBS        ${TRANSI_LIBRARIES}
  CONDITION   ATLAS_HAVE_TRANS )
>>>>>>> Stashed changes
