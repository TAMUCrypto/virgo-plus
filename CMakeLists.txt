cmake_minimum_required(VERSION 3.10)
project(virgo-plus)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-w -march=native -faligned-new -mavx")
set(CMAKE_C_FLAGS "-w -march=native -faligned-new -mavx")

include_directories(lib)

set(PC_OPTION "virgo" CACHE STRING "PC OPTION chosen by the user as the type of the polynomial commitment.")
set_property(CACHE PC_OPTION PROPERTY STRINGS virgo)
if(PC_OPTION STREQUAL "virgo")
    add_definitions(-DUSE_VIRGO)
    link_directories(lib/virgo/lib)
    add_subdirectory(lib/virgo)
endif()

add_subdirectory(src)