cmake_minimum_required(VERSION 3.25)
project(signal_processing VERSION 0.1.0 LANGUAGES CXX)

add_library(signal_processing INTERFACE)
add_library(signal_processing::signal_processing ALIAS signal_processing)
target_compile_features(signal_processing INTERFACE cxx_std_23)
target_include_directories(signal_processing INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

option(SIGNAL_PROCESSING_BUILD_TESTS "Build correctness tests" ${PROJECT_IS_TOP_LEVEL})

if(SIGNAL_PROCESSING_BUILD_TESTS)
    include(CTest)
    add_executable(signal_processing_smoke tests/smoke.cpp)
    target_link_libraries(signal_processing_smoke PRIVATE signal_processing::signal_processing)
    if(MSVC)
        target_compile_options(signal_processing_smoke PRIVATE /W4 /permissive-)
    else()
        target_compile_options(signal_processing_smoke PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wdouble-promotion -Wold-style-cast)
    endif()
    add_test(NAME signal_processing_smoke COMMAND signal_processing_smoke)
endif()
