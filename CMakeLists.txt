# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 3.14.4)

project(MDSW VERSION 0.0.1)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")

include_directories(${OS_INCLUDE_DIRS} "include")
link_directories(${OS_LIB_DIRS})