cmake_minimum_required(VERSION 3.20)
project(ucntrap VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(UCNTRAP_BUILD_TESTS "Build tests" OFF)

add_library(ucntrap_core
    src/core/config.cpp
    src/experiment/production_tracker.cpp
    src/numeric/integrator.cpp
    src/physics/dagger.cpp
    src/physics/planar_halbach_field.cpp
    src/physics/surface_model.cpp
    src/physics/trap_halbach_field.cpp
    src/source/random_source.cpp
)

target_include_directories(ucntrap_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

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