cmake_minimum_required(VERSION 2.8)

project(flightcontrol)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

add_executable(${PROJECT_NAME} "main.cpp")
target_link_libraries (${PROJECT_NAME} ${JSONCPP_LIBRARIES} Eigen3::Eigen glog gflags ceres pthread cholmod lapack blas cxsparse)
