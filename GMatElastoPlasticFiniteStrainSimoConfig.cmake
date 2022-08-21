# GMatElastoPlasticFiniteStrainSimo cmake module
#
# This module sets the target:
#
#     GMatElastoPlasticFiniteStrainSimo
#
# In addition, it sets the following variables:
#
#     GMatElastoPlasticFiniteStrainSimo_FOUND - true if the library is found
#     GMatElastoPlasticFiniteStrainSimo_VERSION - the library's version
#     GMatElastoPlasticFiniteStrainSimo_INCLUDE_DIRS - directory containing the library's headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastoPlasticFiniteStrainSimo::compiler_warnings - enable compiler warnings
#     GMatElastoPlasticFiniteStrainSimo::assert - enable library assertions
#     GMatElastoPlasticFiniteStrainSimo::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastoPlasticFiniteStrainSimo"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElastoPlasticFiniteStrainSimoTargets.cmake")
endif()

# Define "GMatElastoPlasticFiniteStrainSimo_INCLUDE_DIRS"

get_target_property(
    GMatElastoPlasticFiniteStrainSimo_INCLUDE_DIRS
    GMatElastoPlasticFiniteStrainSimo
    INTERFACE_INCLUDE_DIRECTORIES)

# Find dependencies

find_dependency(GMatTensor)
find_dependency(xtensor)

# Define support target "GMatElastoPlasticFiniteStrainSimo::compiler_warnings"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo::compiler_warnings)
    add_library(GMatElastoPlasticFiniteStrainSimo::compiler_warnings INTERFACE IMPORTED)
    target_link_libraries(GMatElastoPlasticFiniteStrainSimo::compiler_warnings INTERFACE
        GMatTensor::compiler_warnings)
endif()

# Define support target "GMatElastoPlasticFiniteStrainSimo::assert"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo::assert)
    add_library(GMatElastoPlasticFiniteStrainSimo::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticFiniteStrainSimo::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT)
endif()

# Define support target "GMatElastoPlasticFiniteStrainSimo::debug"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo::debug)
    add_library(GMatElastoPlasticFiniteStrainSimo::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticFiniteStrainSimo::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT
        GMATTENSOR_ENABLE_ASSERT
        XTENSOR_ENABLE_ASSERT)
endif()
