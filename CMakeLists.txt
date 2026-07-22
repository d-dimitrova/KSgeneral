cmake_minimum_required(VERSION 3.18)

project(ksgeneral_native VERSION 1.0.0 LANGUAGES C CXX)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

add_library(ksgeneral_core OBJECT
    src/k2sample_stable1.cpp
    src/crossprob_new.cc
    src/two_sided_noncrossing_probability.cc
    src/fftwconvolver.cc
    src/fftw_wrappers.cc
    src/poisson_pmf.cc
    src/read_boundaries_file.cc
    src/string_utils.cc
)

set_target_properties(ksgeneral_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(ksgeneral_core PRIVATE cxx_std_11)
target_include_directories(ksgeneral_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(ksgeneral SHARED src/c_api.cpp $<TARGET_OBJECTS:ksgeneral_core>)
target_compile_definitions(ksgeneral PRIVATE KSG_BUILD_DLL)
target_include_directories(ksgeneral
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ksgeneral PRIVATE PkgConfig::FFTW3)
if(UNIX AND NOT APPLE)
    target_link_libraries(ksgeneral PRIVATE m)
endif()
set_target_properties(ksgeneral PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1 CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES)

add_executable(ksg_c_smoke tests/c_smoke.c)
target_include_directories(ksg_c_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ksg_c_smoke PRIVATE ksgeneral)
