# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################


####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_CXX_FLAGS_RELEASE     "-O3" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_CXX_FLAGS_DEBUG       "-O0 -g -traceback -fp-trap=common" )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

set( CMAKE_CXX_FLAGS_BIT         "-O2" )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_CXX_LINK_FLAGS        "" )

####################################################################

# Meaning of flags
# ----------------
# todo

