cmake_minimum_required(VERSION 3.25)
project(signal_processing VERSION 0.1.0 LANGUAGES CXX)

add_library(signal_processing INTERFACE)
add_library(signal_processing::signal_processing ALIAS signal_processing)
target_compile_features(signal_processing INTERFACE cxx_std_23)
target_include_directories(signal_processing INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

option(SIGNAL_PROCESSING_BUILD_TESTS "Build correctness tests" ${PROJECT_IS_TOP_LEVEL})

function(signal_processing_configure_test target)
    target_link_libraries(${target} PRIVATE signal_processing::signal_processing)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wdouble-promotion -Wold-style-cast)
    endif()
endfunction()

if(SIGNAL_PROCESSING_BUILD_TESTS)
    include(CTest)

    add_executable(signal_processing_smoke tests/smoke.cpp)
    signal_processing_configure_test(signal_processing_smoke)
    add_test(NAME signal_processing_smoke COMMAND signal_processing_smoke)

    add_executable(signal_processing_fft tests/fft.cpp)
    signal_processing_configure_test(signal_processing_fft)
    add_test(NAME signal_processing_fft COMMAND signal_processing_fft)

    add_executable(signal_processing_convolution_correlation tests/convolution_correlation.cpp)
    signal_processing_configure_test(signal_processing_convolution_correlation)
    add_test(NAME signal_processing_convolution_correlation
             COMMAND signal_processing_convolution_correlation)

    add_executable(signal_processing_transforms tests/transforms.cpp)
    signal_processing_configure_test(signal_processing_transforms)
    add_test(NAME signal_processing_transforms COMMAND signal_processing_transforms)
endif()
