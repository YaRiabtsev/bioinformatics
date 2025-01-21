cmake_minimum_required(VERSION 3.29)
project(bioinformatics)

set(CMAKE_CXX_STANDARD 20)

add_library(bioinformatics_lib
        src/findstr.cpp
        src/alignstr.cpp
)

target_include_directories(bioinformatics_lib PUBLIC include)

target_precompile_headers(bioinformatics_lib PRIVATE
        include/findstr.hpp
        include/alignstr.hpp
        include/tables.hpp
)

set_target_properties(bioinformatics_lib PROPERTIES UNITY_BUILD ON)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
        googlebenchmark
        URL https://github.com/google/benchmark/archive/refs/tags/v1.7.1.zip
)
FetchContent_MakeAvailable(googlebenchmark)

enable_testing()

add_executable(test_bioinformatics
        tests/test_findstr.cpp
        tests/test_alignstr.cpp
        tests/main.cpp
)
target_link_libraries(test_bioinformatics bioinformatics_lib gtest gtest_main)

add_executable(benchmark_bioinformatics
        benchmarks/benchmark_findstr.cpp
)
target_link_libraries(benchmark_bioinformatics bioinformatics_lib benchmark::benchmark)