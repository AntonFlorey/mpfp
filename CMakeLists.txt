
cmake_minimum_required(VERSION 3.22)
project("mpfp")
set(CMAKE_CXX_STANDARD 20)

#set (PYBIND11_PYTHON_VERSION "3.11")

if (MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif (MSVC)

# Eigen
find_package (Eigen3 REQUIRED NO_MODULE)
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
)

target_compile_definitions(mpfp
                           PRIVATE VERSION_INFO=${MPFP_VERSION_INFO})
