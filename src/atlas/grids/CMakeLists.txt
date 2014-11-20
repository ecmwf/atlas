list( APPEND atlas_grids_srcs
Grids.h
reduced_gg/Grids.h
reduced_gg/N16.cc    # TL31
reduced_gg/N24.cc    # TL47
reduced_gg/N32.cc    # TL63
reduced_gg/N48.cc    # TL95
reduced_gg/N64.cc    # TL127
reduced_gg/N80.cc    # TL159
reduced_gg/N96.cc    # TL191
reduced_gg/N128.cc   # TL255
reduced_gg/N160.cc   # TL319
reduced_gg/N200.cc   # TL399
reduced_gg/N256.cc   # TL511
reduced_gg/N320.cc   # TL639
reduced_gg/N400.cc   # TL799
reduced_gg/N512.cc   # TL1023
reduced_gg/N576.cc   # TL1151
reduced_gg/N640.cc   # TL1279
reduced_gg/N800.cc   # TL1599
reduced_gg/N1024.cc  # TL2047
reduced_gg/N1280.cc  # TL2559
reduced_gg/N1600.cc  # TL3199
reduced_gg/N2000.cc  # TL3999
reduced_gg/N4000.cc  # TL7999
reduced_gg/N8000.cc  # TL15999
)

ecbuild_add_library( TARGET atlas_grids
  AUTO_VERSION
  SOURCES ${atlas_grids_srcs}
  LIBS ${libatlas}
)

ecbuild_add_library( TARGET atlas_grids_fortran
  CONDITION CMAKE_Fortran_COMPILER_LOADED
  SOURCES Grids.cc atlas_grids_module.F90
  LIBS atlas_grids
)

