project(ocra CXX)

include_directories(BEFORE ${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE control_srcs src/control/*.cpp)
file(GLOB_RECURSE optim_srcs src/optim/*.cpp)
file(GLOB_RECURSE util_srcs src/util/*.cpp)

file(GLOB_RECURSE control_headers include/ocra/control/*.h)
file(GLOB_RECURSE optim_headers include/ocra/optim/*.h)
file(GLOB_RECURSE util_headers include/ocra/util/*.h)

# Define groupings for source/headers files in IDE project generation
source_group("ocra-control-src" FILES ${control_srcs})
source_group("ocra-control-headers" FILES ${control_headers})
source_group("ocra-optim-src" FILES ${optim_srcs})
source_group("ocra-optim-headers" FILES ${optim_headers})
source_group("ocra-util-src" FILES ${util_srcs})
source_group("ocra-util-headers" FILES ${util_headers})


set(OCRA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE PATH "")

# Check for SMLT
find_package(SMLT QUIET)
IF(${SMLT_FOUND})
    message(STATUS "SMLT lib found. Compiling gaussian process trajectory.")
    include_directories(${SMLT_INCLUDE_DIRS})
    add_definitions(-DUSING_SMLT=1)
ELSE()
    message(STATUS "SMLT lib NOT found. Removing gaussian process trajectory from source list.")
    list(REMOVE_ITEM control_srcs ${CMAKE_CURRENT_SOURCE_DIR}/src/control/Trajectory/GaussianProcessTrajectory.cpp)
ENDIF()

# Check for GTTraj
find_package(GTTraj QUIET)
IF(${GTTraj_FOUND})
    message(STATUS "GTTraj lib found. Compiling time optimal trajectory.")
    include_directories(${GTTraj_INCLUDE_DIRS})
    add_definitions(-DUSING_GTTRAJ=1)
ELSE()
    message(STATUS "GTTraj lib NOT found. Removing time optimal trajectory from source list.")
    list(REMOVE_ITEM control_srcs ${CMAKE_CURRENT_SOURCE_DIR}/src/control/Trajectory/TimeOptimalTrajectory.cpp)
ENDIF()



add_library(ocra-control ${control_srcs} ${control_headers})
add_library(ocra-optim ${optim_srcs} ${optim_headers})
add_library(ocra-util ${util_srcs} ${util_headers})

LIST(APPEND link_libs   ${Boost_LIBRARIES}
                        ${YARP_LIBRARIES}
                        ${qpOASES_LIBRARIES}
                        )
IF(${SMLT_FOUND})
 LIST(APPEND link_libs ${SMLT_LIBRARIES})
ENDIF()

IF(${GTTraj_FOUND})
 LIST(APPEND link_libs ${GTTraj_LIBRARIES})
ENDIF()

target_link_libraries(ocra-control ${link_libs})
target_link_libraries(ocra-optim ${link_libs})
target_link_libraries(ocra-util ${link_libs})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        FILES_MATCHING PATTERN *.h*)

install(TARGETS ocra-control ocra-optim ocra-util
        EXPORT OcraRecipes
        DESTINATION lib)

set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS ocra-control ocra-optim ocra-util)
# set_property(GLOBAL APPEND PROPERTY PACKAGE_${VARS_PREFIX}_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
