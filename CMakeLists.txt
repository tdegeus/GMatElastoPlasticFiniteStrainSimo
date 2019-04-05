
# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# project settings
# ----------------

# name
project(GMatElastoPlasticFiniteStrainSimo)

# file that contains the version information
set(parse_version include/GMatElastoPlasticFiniteStrainSimo/config.h)

# header files
set(headers
  include/GMatElastoPlasticFiniteStrainSimo/config.h
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d.hpp
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d_eigvals.hpp
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d_Elastic.hpp
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d_LinearHardening.hpp
  include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d_Matrix.hpp
)

# automatically parse the version number
file(READ "${parse_version}" version)
string(REGEX MATCH "define[ \t]+GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GMATELASTOPLASTICFINITESTRAINSIMO_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(GMATELASTOPLASTICFINITESTRAINSIMO_VERSION ${world_version}.${major_version}.${minor_version})

# paths
# -----

set(GMATELASTOPLASTICFINITESTRAINSIMO_ROOT_DIR "${CMAKE_INSTALL_PREFIX}")
set(GMATELASTOPLASTICFINITESTRAINSIMO_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/pkgconfig")
set(fcmake "GMatElastoPlasticFiniteStrainSimoConfig.cmake")
set(fpkg "GMatElastoPlasticFiniteStrainSimo.pc")

# options
# -------

# configure pkg-config (default: on)
option(PKGCONFIG "Build pkg-config ${fpkg} file" ON)

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${fpkg} file" OFF)
endif()

# C++ standard
# ------------

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# configure CMake
# ---------------

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticFiniteStrainSimoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticFiniteStrainSimoConfig.cmake
  PATH_VARS GMATELASTOPLASTICFINITESTRAINSIMO_INCLUDE_DIR GMATELASTOPLASTICFINITESTRAINSIMO_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# -------

# pkg-config
if(PKGCONFIG)
  configure_file(${fpkg}.in ${fpkg} @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fpkg} DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# CMake
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fcmake} DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

# header files
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${headers} DESTINATION ${INCLUDE_INSTALL_DIR})

# print information to screen
# ---------------------------

message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in")
message(STATUS "|   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "|")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake /path/to/${PROJECT_NAME} -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "|")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
