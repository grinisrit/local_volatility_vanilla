cmake_minimum_required(VERSION 3.16)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

project(local_volatility_vanilla)

set( CMAKE_CXX_STANDARD          20  )
set( CMAKE_CXX_STANDARD_REQUIRED ON  )
set( CMAKE_CXX_EXTENSIONS        OFF )

set( CMAKE_BUILD_TYPE Release)

include_directories(include)


find_package(Enzyme REQUIRED CONFIG PATHS /opt/homebrew/Cellar/enzyme)
message(STATUS "Found EnzymeAD: ${Enzyme_DIR}")

add_executable(enzyme_test src/main.cc src/euler_calc.cc)
target_link_libraries(enzyme_test PRIVATE LLDEnzymeFlags ) 
target_compile_options(enzyme_test PRIVATE -O3)


target_compile_features(enzyme_test PRIVATE cxx_std_20)