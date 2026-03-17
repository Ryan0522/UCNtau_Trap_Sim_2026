cmake_minimum_required(VERSION 3.20)
project(ucntrap VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- 1. Directories setting ---
# Remove multi-configuration (Debug, Release, etc.)
set(BASE_BIN_DIR "${CMAKE_SOURCE_DIR}/bin")
set(BASE_LIB_DIR "${CMAKE_SOURCE_DIR}/lib")
set(BASE_PY_DIR  "${CMAKE_SOURCE_DIR}/python")

foreach(config ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${config} CONFIG_UPPER)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} "${BASE_BIN_DIR}")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG_UPPER} "${BASE_LIB_DIR}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG_UPPER} "${BASE_LIB_DIR}")
endforeach()

# All .exe to /bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
# All lib to /lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

include(FetchContent)

# --- 2. fetch (yaml-cpp & pybind11) --- 
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.9.0
)
FetchContent_MakeAvailable(yaml-cpp)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# --- 3. set MPI --- 
find_package(MPI REQUIRED)

# --- 4. (ucntrap_core) --- 
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

target_include_directories(ucntrap_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ucntrap_core PUBLIC yaml-cpp::yaml-cpp MPI::MPI_CXX)

# --- 5. Python wrapper (ucntrap_py) --- 
pybind11_add_module(ucntrap_py src/python_wrapper.cpp)
target_link_libraries(ucntrap_py PRIVATE ucntrap_core)

# Python wrapper exported to /python
set_target_properties(ucntrap_py PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python
)

# --- 6. main cli (CLI) --- [cite: 8, 108]
add_executable(ucntrap src/cli/main.cpp)
target_link_libraries(ucntrap PRIVATE ucntrap_core)

# --- 7. testing functions (bin/tests) --- 
enable_testing()
function(add_ucn_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE ucntrap_core)
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${config} CONFIG_UPPER)
        set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${CONFIG_UPPER} "${BASE_BIN_DIR}/tests")
    endforeach()
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BASE_BIN_DIR}/tests")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_ucn_test(traj_overlap_test tests/traj_overlap_test.cpp)
add_ucn_test(liouville_test tests/liouville_test.cpp)
add_ucn_test(compare_dist_test tests/compare_dist.cpp)