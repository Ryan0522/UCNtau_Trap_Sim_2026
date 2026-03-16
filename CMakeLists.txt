cmake_minimum_required(VERSION 3.20)
project(ucntrap VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(UCNTRAP_BUILD_TESTS "Build tests" OFF)

add_library(ucntrap_core
    src/core/config.cpp
    src/core/runner.cpp
    src/experiment/production_tracker.cpp
    src/io/result_writer.cpp
    src/io/trace_loader.cpp
    src/numeric/integrator.cpp
    src/physics/dagger.cpp
    src/physics/planar_halbach_field.cpp
    src/physics/surface_model.cpp
    src/physics/trap_halbach_field.cpp
    src/source/pentrack_source.cpp
    src/source/random_source.cpp
)

target_include_directories(ucntrap_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

include(FetchContent)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.9.0
)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(yaml-cpp)

target_link_libraries(ucntrap_core PUBLIC yaml-cpp::yaml-cpp)

set(MPI_HOME "C:/Program Files (x86)/Microsoft SDKs/MPI")
set(MPI_GUESS_LIBRARY_NAME msmpi)

find_package(MPI REQUIRED)

if(MPI_FOUND)
    message(STATUS "MPI found: ${MPI_CXX_INCLUDE_DIRS}")
    target_link_libraries(ucntrap_core PUBLIC MPI::MPI_CXX)
else()
    target_include_directories(ucntrap_core PUBLIC "${MPI_HOME}/Include")
    target_link_libraries(ucntrap_core PUBLIC "${MPI_HOME}/Lib/x64/msmpi.lib")
endif()

add_executable(ucntrap
    src/cli/main.cpp
)

target_link_libraries(ucntrap
    PRIVATE
        ucntrap_core
)

if(UCNTRAP_BUILD_TESTS)
    enable_testing()

    add_executable(test_config tests/test_config.cpp)
    target_link_libraries(test_config PRIVATE ucntrap_core)
    add_test(NAME test_config COMMAND test_config)
endif()