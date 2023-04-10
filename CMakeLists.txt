cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(optimal_control_tasks)

if (WIN32)
    set(ENV{GSL_ROOT_DIR} ${GSL_INSTALL_DIR})
    if (EXISTS $ENV{GSL_ROOT_DIR})
        message(STATUS "gsl root directory $ENV{GSL_ROOT_DIR} exists")
    else()
        message(ERROR "gsl root directory $ENV{GSL_ROOT_DIR} does not exist")
    endif()
endif()

find_package(GSL REQUIRED)

subdirs(
    mathutils
    task1
)