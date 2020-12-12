
cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastoPlasticFiniteStrainSimo)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastoPlasticFiniteStrainSimo/config.h" _GMAT_ver
     REGEX "#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${_GMAT_ver})
    if(ver MATCHES
       "#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(_GMAT_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(_GMAT_VERSION ${_GMAT_VERSION_MAJOR}.${_GMAT_VERSION_MINOR}.${_GMAT_VERSION_PATCH})

message(STATUS "Building GMatElastoPlasticFiniteStrainSimo v${_GMAT_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)
find_package(GMatTensor REQUIRED)

add_library(GMatElastoPlasticFiniteStrainSimo INTERFACE)

target_include_directories(GMatElastoPlasticFiniteStrainSimo INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastoPlasticFiniteStrainSimo INTERFACE xtensor GMatTensor)

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

set(_GMAT ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimoConfigVersion.cmake"
    VERSION ${_GMAT_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMAT})

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

set(CTEST_TEST_TARGET_ALIAS RUN_ALL_TESTS)
include(CTest)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
