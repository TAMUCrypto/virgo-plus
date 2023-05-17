aux_source_directory(. SRC)
include_directories(.)
add_executable(virgo_plus_run ${SRC})

if(PC_OPTION STREQUAL "virgo")
    target_link_libraries(virgo_plus_run virgo_lib XKCP)
endif()
