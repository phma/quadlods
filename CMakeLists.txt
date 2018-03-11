project(quadlods)
cmake_minimum_required(VERSION 2.8)

add_executable(quadlods main.cpp circletest.cpp contfrac.cpp discrepancy.cpp histogram.cpp
               ldecimal.cpp manysum.cpp polyline.cpp ps.cpp xy.cpp)
add_library(quadlib0 STATIC quadlods.cpp)
add_library(quadlib1 SHARED quadlods.cpp)
add_custom_command(OUTPUT primes.dat COMMAND quadlods ARGS sortprimes MAIN_DEPENDENCY quadlods)
set(CMAKE_CXX_FLAGS "-std=c++11")
install(TARGETS quadlods DESTINATION bin)
install(TARGETS quadlib0 quadlib1 DESTINATION lib)
install(FILES quadlods.h DESTINATION include)
install(FILES primes.dat DESTINATION share/quadlods)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(GMPXX REQUIRED)
set(LIBS ${LIBS} ${GMP_LIBRARY} ${GMPXX_LIBRARIES})
target_link_libraries(quadlods ${LIBS} quadlib1)

set(QUADLODS_MAJOR_VERSION 0)
set(QUADLODS_MINOR_VERSION 1)
set(QUADLODS_PATCH_VERSION 0)
set(QUADLODS_VERSION ${QUADLODS_MAJOR_VERSION}.${QUADLODS_MINOR_VERSION}.${QUADLODS_PATCH_VERSION})

set_target_properties(quadlib0 PROPERTIES OUTPUT_NAME "quadlods" VERSION ${QUADLODS_VERSION})
set_target_properties(quadlib1 PROPERTIES OUTPUT_NAME "quadlods" VERSION ${QUADLODS_VERSION})

set(CPACK_PACKAGE_VERSION_MAJOR ${QUADLODS_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${QUADLODS_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${QUADLODS_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)
