cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Adding vcpkg overlay ports with custom ports for dependencies: branch-and-bound, ModernDD, or-utils, docopt
set(VCPKG_OVERLAY_PORTS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg-custom-ports")
#
# Project details
#

project(
  PM
  VERSION 0.1.0
  LANGUAGES CXX C
  DESCRIPTION "Branch-and-Price algorithm for parallel machine scheduling problems"
)

#
# Set project options
#

include(cmake/ConfigSafeGuards.cmake)
include(cmake/StandardSettings.cmake)
include(cmake/StaticAnalyzers.cmake)
include(cmake/Utils.cmake)

# For Windows, it is necessary to link with the MultiThreaded library. Depending on how the rest of
# the project's dependencies are linked, it might be necessary to change the line to statically link
# with the library.
#
# This is done as follows:
#
# set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
#
# On Linux and Mac this variable is ignored. If any issues rise from it, try commenting it out and
# letting CMake decide how to link with it.
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...")

if(UNIX)
  add_compile_options(
    "$<$<CONFIG:DEBUG>:-D_DEBUG>"
  ) # this will allow to use same _DEBUG macro available in both Linux as well as Windows - MSCV
    # environment. Easy to put Debug specific code.
endif(UNIX)

#
# Setup alternative names
#

if(${PROJECT_NAME}_USE_ALT_NAMES)
  string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)
  string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)
else()
  set(PROJECT_NAME_LOWERCASE ${PROJECT_NAME})
  set(PROJECT_NAME_UPPERCASE ${PROJECT_NAME})
endif()

#
# Create library, setup header and source files
#

# Find all headers and implementation files
include(cmake/SourcesAndHeaders.cmake)

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME} ${exe_sources})

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following sources:")
    foreach(source IN LISTS exe_sources)
      verbose_message("* ${source}")
    endforeach()
  endif()

  if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
    add_library(${PROJECT_NAME}_LIB ${headers} ${sources})

    if(${PROJECT_NAME}_VERBOSE_OUTPUT)
      verbose_message("Found the following headers:")
      foreach(header IN LISTS headers)
        verbose_message("* ${header}")
      endforeach()
    endif()
  endif()
elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  add_library(${PROJECT_NAME} INTERFACE)

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following headers:")
    foreach(header IN LIST headers)
      verbose_message("* ${header}")
    endforeach()
  endif()
else()
  add_library(${PROJECT_NAME} ${headers} ${sources})

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following sources:")
    foreach(source IN LISTS sources)
      verbose_message("* ${source}")
    endforeach()
    verbose_message("Found the following headers:")
    foreach(header IN LISTS headers)
      verbose_message("* ${header}")
    endforeach()
  endif()
endif()

message(STATUS "Added all header and implementation files.")

#
# Set the project standard and warnings
#
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++20" COMPILER_SUPPORTS_CXX20)
check_cxx_compiler_flag("-std=gnu++20" COMPILER_SUPPORTS_CXX20_GNU)
check_cxx_compiler_flag("-std:c++20" COMPILER_SUPPORTS_CXX20_MSVC)
check_cxx_compiler_flag("-std=c++2a" COMPILER_SUPPORTS_CXX2A)

if(${PROJECT_NAME}_BUILD_HEADERS_ONLY
   AND (COMPILER_SUPPORTS_CXX20
        OR COMPILER_SUPPORTS_CXX20_GNU
        OR COMPILER_SUPPORTS_CXX20_MSVC
        OR COMPILER_SUPPORTS_CXX2A)
)
  target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
elseif(
  COMPILER_SUPPORTS_CXX20
  OR COMPILER_SUPPORTS_CXX20_GNU
  OR COMPILER_SUPPORTS_CXX20_MSVC
  OR COMPILER_SUPPORTS_CXX2A
)
  target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
else()
  message(
    FATAL_ERROR
      "The compiler ${CMAKE_CXX_COMPILER} has no C++20 support. Please use a different C++ compiler."
  )
endif()

include(cmake/CompilerWarnings.cmake)
set_project_warnings(${PROJECT_NAME})

message(STATUS "Applied compiler warnings. Using standard C++20.")

#
# Model project dependencies
#

# Identify and link with the specific "packages" the project uses find_package(package_name
# package_version REQUIRED package_type [other_options])
set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(fmt CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS atomic chrono timer random regex serialization)
find_package(docopt CONFIG REQUIRED)
find_package(range-v3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(or-utils CONFIG REQUIRED)
find_package(ModernDD CONFIG REQUIRED)
find_package(branch-and-bound CONFIG REQUIRED)

include(cmake/FindGurobi.cmake)
include(cmake/SetUpOsi.cmake)

target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE $<IF:$<CONFIG:DEBUG>,${GUROBI_CXX_DEBUG_LIBRARY},${GUROBI_CXX_LIBRARY}>
          $<IF:$<CONFIG:DEBUG>,${OSI_LIBRARY_DEBUG},${OSI_LIBRARY}>
          $<IF:$<CONFIG:DEBUG>,${OSIGRB_LIBRARY_DEBUG},${OSIGRB_LIBRARY}>
          $<IF:$<CONFIG:DEBUG>,${COINUTILS_LIBRARY_DEBUG},${COINUTILS_LIBRARY}>
          $<IF:$<CONFIG:DEBUG>,${CGL_LIBRARY_DEBUG},${CGL_LIBRARY}>
          $<$<LINK_LANG_AND_ID:CXX,GNU>:docopt_s>
          $<$<LINK_LANG_AND_ID:CXX,MSVC>:docopt>
          ${GUROBI_LIBRARY}
          fmt::fmt
          Boost::atomic
          Boost::chrono
          Boost::timer
          Boost::random
          range-v3::range-v3
          nlohmann_json::nlohmann_json
          or-utils::or-utils
          branch-and-bound::branch-and-bound
          ModernDD::ModernDD
)

verbose_message("Successfully added all dependencies and linked against them.")

#
# Set the build/user include directories
#

# Allow usage of header files in the `src` directory, but only for utilities
if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_include_directories(
    ${PROJECT_NAME} INTERFACE $<INSTALL_INTERFACE:include>
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
else()
  target_include_directories(
    ${PROJECT_NAME}
    PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ModernDD::ModernDD
            ${GUROBI_INCLUDE_DIR}
            ${OSI_INCLUDE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${Boost_INCLUDES}
            ${fmt_INCLUDES}
            ${docopt_INCLUDES}
            ${range-v3_INCLUDES}
            ${nlohmann_json_INCLUDES}
  )

  if(${PROJECT_NAME}_BUILD_EXECUTABLE AND ${PROJECT_NAME}_ENABLE_UNIT_TESTING)
    target_include_directories(
      ${PROJECT_NAME}_LIB
      PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      PRIVATE ModernDD::ModernDD
              ${GUROBI_INCLUDE_DIR}
              ${OSI_INCLUDE_DIR}
              ${Boost_INCLUDES}
              ${fmt_INCLUDES}
              ${docopt_INCLUDES}
              ${range-v3_INCLUDES}
              ${nlohmann_json_INCLUDES}
              ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    target_link_libraries(
      ${PROJECT_NAME}_LIB
      PRIVATE ${GUROBI_LIBRARY}
              ${GUROBI_CXX_LIBRARY}
              ${Boost_INCLUDES}
              ${fmt_INCLUDES}
              ${docopt_INCLUDES}
              ${range-v3_INCLUDES}
              ${nlohmann_json_INCLUDES}
              ORUtils::ORUtils
              branch-and-bound::branch-and-bound
              ModernDD::ModernDD
              ${OSI_LIBRARY}
              ${OSI_GRB_LIBRARY}
              ${COINUtils_LIBRARY}
              ${CGL_LIBRARY}
              -lm
              -lpthread
    )
  endif()
endif()

message(STATUS "Finished setting up include directories.")

#
# Provide alias to library for
#

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
else()
  add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
endif()

verbose_message("Project is now aliased as ${PROJECT_NAME}::${PROJECT_NAME}.")

#
# Format the project using the `clang-format` target (i.e: cmake --build build --target
# clang-format)
#


#
# Install library for easy downstream inclusion
#

include(GNUInstallDirs)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES
  DESTINATION include
  PUBLIC_HEADER DESTINATION include
)

install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

#
# Quick `ConfigVersion.cmake` creation
#

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# include Doxygen
include(cmake/Doxygen.cmake)

# Add Gurobi env to Binary directory
file(COPY ${CMAKE_SOURCE_DIR}/settings/gurobi.env DESTINATION ${CMAKE_BINARY_DIR})
