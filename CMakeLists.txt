#cmake_minimum_required(VERSION 3.14)
cmake_minimum_required(VERSION 3.10)
project(EMDKL)

set(CMAKE_CXX_STANDARD 14)


list(APPEND CMAKE_PREFIX_PATH "/Users/klapidus/root/obj")
find_package(ROOT COMPONENTS Hist)
include(${ROOT_USE_FILE})

find_path(FASTJET_INCLUDE_DIR fastjet/version.hh
        HINTS "/Users/klapidus/ssoft/fastjet-install/include")

find_library(FASTJET_LIBRARY NAMES fastjet
        HINTS "/Users/klapidus/ssoft/fastjet-install/lib")

add_library(EMDKL SHARED EMDKLAnalyser.cpp EMDKLAnalyser.h)

target_include_directories(EMDKL PUBLIC ${FASTJET_INCLUDE_DIR})
target_link_libraries(EMDKL ${ROOT_LIBRARIES} ${FASTJET_LIBRARY})