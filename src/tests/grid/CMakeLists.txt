# =========================================================================================
# Boost  
# =========================================================================================

# If this is enable then can not find boost libs using env BOOST_ROOT ?
#ecbuild_add_extra_search_paths( boost )

# To use static boost python ensure that Boost_USE_STATIC_LIBS is set on.
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_NO_SYSTEM_PATHS        ON)
set(Boost_DETAILED_FAILURE_MSG   ON)
#set(Boost_DEBUG                 ON)

find_package( Boost 1.53.0 COMPONENTS  unit_test_framework test_exec_monitor  )

# Available boost lib should be referenced as:
#
#    ${Boost_SYSTEM_LIBRARY}
#    ${Boost_SERIALIZATION_LIBRARY}
#    ${Boost_THREAD_LIBRARY}
#    ${Boost_FILESYSTEM_LIBRARY}
#    ${Boost_PROGRAM_OPTIONS_LIBRARY}
#    ${Boost_DATE_TIME_LIBRARY}
#
#    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} 
#    ${Boost_TEST_EXEC_MONITOR_LIBRARY} 

message( STATUS "Boost_LIBRARIES     : ${Boost_LIBRARIES}" )

# =====================================================================================
# Tests
# =====================================================================================
ecbuild_add_test(   TARGET      test_RegularLatLonGrid
               BOOST
               SOURCES     test_RegularLatLonGrid.cc
               CONDITION   HAVE_GRID
               INCLUDES    ${Boost_INCLUDE_DIRS} 
               LIBS        atlas_grid  
                           ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} 
                           ${Boost_TEST_EXEC_MONITOR_LIBRARY}
                )

# ======================================================================================
ecbuild_add_test(   TARGET      test_latlon
               SOURCES     test_latlon.cc
               CONDITION   HAVE_GRID
               LIBS        atlas_grid  )
            
ecbuild_add_test(   TARGET      test_gaussian
					SOURCES     test_gaussian.cc
					CONDITION	HAVE_GRID
					LIBS        atlas_grid )

ecbuild_add_test(   TARGET      test_field
					SOURCES     test_field.cc
					CONDITION	HAVE_GRID
					LIBS        atlas_grid )

# this should test the eckit::maths expressions
#ecbuild_add_test(   TARGET      test_action
#		SOURCES     test_action.cc
#		CONDITION	HAVE_GRID
#		LIBS        eckit atlas_grid )

 