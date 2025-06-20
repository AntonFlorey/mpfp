
cmake_minimum_required(VERSION 3.22)
project("mpfp")
set(CMAKE_CXX_STANDARD 20)

#set (PYBIND11_PYTHON_VERSION "3.11")

if (MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif (MSVC)

# OpenMP
if (NOT TARGET OpenMP::OpenMP_CXX)
    if(APPLE)
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
    endif()
    find_package(OpenMP REQUIRED)
endif()

# Eigen
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen"
)
# pybind11
add_subdirectory(extern/pybind11)
# TinyAD
add_subdirectory(extern/TinyAD)

# Include dirs
include_directories("extern/TinyAD/include")
include_directories("src")

# message(${Eigen3_DIR})

# module (pybind11)
pybind11_add_module (
    mpfp
    src/module.cpp
)

target_link_libraries(
    mpfp
    PRIVATE Eigen3::Eigen
    PRIVATE OpenMP::OpenMP_CXX
)

if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(
        mpfp
        PRIVATE OpenMP::OpenMP_CXX
    )
endif()

target_compile_definitions(mpfp
                           PRIVATE VERSION_INFO=${MPFP_VERSION_INFO})
