cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastoPlasticFiniteStrainSimo)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastoPlasticFiniteStrainSimo/config.h"
             GMatElastoPlasticFiniteStrainSimo_version_defines
     REGEX "#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GMatElastoPlasticFiniteStrainSimo_version_defines})
    if(ver MATCHES
       "#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_${CMAKE_MATCH_1}
            "${CMAKE_MATCH_2}"
            CACHE INTERNAL "")
    endif()
endforeach()

set(GMATELASTOPLASTICFINITESTRAINSIMO_VERSION
    ${GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_MAJOR}.${GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_MINOR}.${GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_PATCH})

message(STATUS "Building GMatElastoPlasticFiniteStrainSimo v${GMATELASTOPLASTICFINITESTRAINSIMO_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(GMatElastoPlasticFiniteStrainSimo INTERFACE)

target_include_directories(GMatElastoPlasticFiniteStrainSimo INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastoPlasticFiniteStrainSimo INTERFACE xtensor)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

install(TARGETS GMatElastoPlasticFiniteStrainSimo EXPORT GMatElastoPlasticFiniteStrainSimo-targets)

install(
    EXPORT GMatElastoPlasticFiniteStrainSimo-targets
    FILE GMatElastoPlasticFiniteStrainSimoTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlasticFiniteStrainSimo")

set(_GMATELASTOPLASTICFINITESTRAINSIMO ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimoConfigVersion.cmake"
    VERSION ${GMATELASTOPLASTICFINITESTRAINSIMO_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMATELASTOPLASTICFINITESTRAINSIMO})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticFiniteStrainSimoConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimoConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlasticFiniteStrainSimo")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticFiniteStrainSimo.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimo.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimo.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GMatElastoPlasticFiniteStrainSimoConfig.cmake")

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
