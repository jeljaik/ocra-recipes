project(taskManagerTests)
file(GLOB source_files *.cpp)
file(GLOB header_files *.h*)
include_directories(${OCRA_INCLUDE_DIRS})
add_executable(${PROJECT_NAME} ${source_files} ${header_files})
target_link_libraries(ocra-control)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
