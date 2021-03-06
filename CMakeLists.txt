cmake_minimum_required(VERSION 3.9)
project(MetaOpt)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-DHASH_VALUE -O2)

# set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
# set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
# set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,-R/usr/local/lib -Wl,-R/usr/local/lib64")
# find_package(Eigen3 REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include src)
link_directories(.)

aux_source_directory(src            SRC)
aux_source_directory(src/doe        SRC_DOE)
aux_source_directory(src/model      SRC_MOD)
aux_source_directory(src/optimizer  SRC_OPT)

add_library(metaopt SHARED ${SRC} ${SRC_DOE} ${SRC_MOD} ${SRC_OPT})
target_link_libraries(metaopt glog)

# add_executable(numerical example/numerical_case/numerical.cc )
# target_link_libraries(numerical glog metaopt)


find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
aux_source_directory(test      TEST)
enable_testing()
foreach(i ${TEST})
  STRING(REGEX REPLACE ".+/(.+)\\..*" "\\1" CASE_NAME ${i})
  message(STATUS "Test case : ${CASE_NAME} - ${i}")
  add_executable(${CASE_NAME} ${i})
  target_link_libraries(${CASE_NAME} gtest_main gtest pthread glog metaopt)
  add_test(${CASE_NAME} ${CASE_NAME})
endforeach(i)
