# GMatElastoPlasticFiniteStrainSimo cmake module
#
# This module sets the target:
#
#     GMatElastoPlasticFiniteStrainSimo
#
# In addition, it sets the following variables:
#
#     GMatElastoPlasticQPot_FOUND - true if GMatElastoPlasticFiniteStrainSimo found
#     GMatElastoPlasticQPot_VERSION - GMatElastoPlasticFiniteStrainSimo's version
#     GMatElastoPlasticQPot_INCLUDE_DIRS - the directory containing GMatElastoPlasticFiniteStrainSimo headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastoPlasticFiniteStrainSimo::compiler_warnings - enable compiler warnings
#     GMatElastoPlasticFiniteStrainSimo::assert - enable GMatElastoPlasticFiniteStrainSimo assertions
#     GMatElastoPlasticFiniteStrainSimo::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastoPlasticFiniteStrainSimo"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElastoPlasticQPotTargets.cmake")
    get_target_property(
        GMatElastoPlasticQPot_INCLUDE_DIRS
        GMatElastoPlasticFiniteStrainSimo
        INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies

find_dependency(xtensor)

# Define support target "GMatElastoPlasticFiniteStrainSimo::compiler_warnings"

if(NOT TARGET GMatElastoPlasticFiniteStrainSimo::compiler_warnings)
    add_library(GMatElastoPlasticFiniteStrainSimo::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GMatElastoPlasticFiniteStrainSimo::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GMatElastoPlasticFiniteStrainSimo::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
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
        XTENSOR_ENABLE_ASSERT GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT)
endif()
