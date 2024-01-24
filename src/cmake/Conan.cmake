if(${PROJECT_NAME}_ENABLE_CONAN)
  #
  # Setup Conan requires and options here:
  #
  set(CONAN_REQUIRES boost/1.73.0 fmt/8.0.0 docopt.cpp/0.6.3 range-v3/0.11.0 nlohmann_json/3.10.0)
  set(CONAN_OPTIONS
      "
    boost:without_context=True
    boost:without_contract=True
    boost:without_coroutine=True
    boost:without_date_time=True
    boost:without_exception=True
    boost:without_fiber=True
    boost:without_filesystem=True
    boost:without_graph_parallel=True
    boost:without_iostreams=True
    boost:without_locale=True
    boost:without_log=True
    boost:without_mpi=True
    boost:without_nowide=True
    boost:without_program_options=True
    boost:without_python=True
    boost:without_stacktrace=True
    boost:without_test=True
    boost:without_thread=True
    boost:without_type_erasure=True
    boost:without_wave=True"
  )

  #
  # If `conan.cmake` (from https://github.com/conan-io/cmake-conan) does not exist, download it.
  #
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan...")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
         "${CMAKE_BINARY_DIR}/conan.cmake"
    )
    message(STATUS "Cmake-Conan downloaded successfully.")
  endif()

  include(${CMAKE_BINARY_DIR}/conan.cmake)

  conan_cmake_configure(
    REQUIRES
    ${CONAN_REQUIRES}
    OPTIONS
    ${CONAN_OPTIONS}
    GENERATORS
    cmake_find_package
    cmake_paths
    cmake
    visual_studio
  )
  conan_cmake_autodetect(settings)
  conan_cmake_install(PATH_OR_REFERENCE . BUILD missing SETTINGS ${settings})

  verbose_message("Conan is setup and all requires have been installed. ")
endif()
