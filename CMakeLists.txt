cmake_minimum_required(VERSION 3.20)
project(ucntrap VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)

# --- 1. Fetch yaml-cpp ---
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.9.0
)
FetchContent_MakeAvailable(yaml-cpp)

# --- 2. Fetch pybind11 ---
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1 # 建議指定版本
)
FetchContent_MakeAvailable(pybind11)

# --- 3. setup MPI ---
find_package(MPI REQUIRED)

# --- 4. core libraries ---
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

# --- 5. Python module ---
pybind11_add_module(ucntrap_py src/python_wrapper.cpp)
target_link_libraries(ucntrap_py PRIVATE ucntrap_core)

# --- 6. CLI executables ---
add_executable(ucntrap src/cli/main.cpp)
target_link_libraries(ucntrap PRIVATE ucntrap_core)

# --- 7. testing ---
enable_testing()
add_executable(physics_overlap_test tests/physics_overlap_test.cpp)
target_link_libraries(physics_overlap_test PRIVATE ucntrap_core)

add_test(NAME PhysicsOverlap COMMAND physics_overlap_test)