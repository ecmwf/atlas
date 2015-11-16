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

ecbuild_add_cxx_flags( "-O3"  BUILD RELEASE )

####################################################################
# DEBUG FLAGS
####################################################################

ecbuild_add_cxx_flags( "-O0 -g" BUILD DEBUG )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

ecbuild_add_cxx_flags( "-O2" BUILD BIT )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_CXX_LINK_FLAGS        "" )

####################################################################

# Meaning of flags
# ----------------
# todo

