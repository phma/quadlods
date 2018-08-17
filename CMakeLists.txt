project(quadlods)
cmake_minimum_required(VERSION 2.8)

add_executable(quadlods main.cpp circletest.cpp contfrac.cpp filltest.cpp
               histogram.cpp hstep.cpp ldecimal.cpp manysum.cpp matrix.cpp polyline.cpp ps.cpp
               random.cpp xy.cpp)
add_library(quadlib0 STATIC quadlods.cpp)
add_library(quadlib1 SHARED quadlods.cpp)
add_custom_command(OUTPUT primes.dat COMMAND quadlods sortprimes)
#primes.dat does not depend on quadlods because it's just math and doesn't change.
add_custom_target(primesdat ALL DEPENDS primes.dat)
install(TARGETS quadlods DESTINATION bin)
install(TARGETS quadlib0 quadlib1 DESTINATION lib)
install(FILES quadlods.h DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/primes.dat DESTINATION share/quadlods)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/quadlods)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(GMPXX REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)
set(LIBS ${LIBS} ${GMP_LIBRARY} ${GMPXX_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(quadlods ${LIBS} quadlib1)

set(QUADLODS_MAJOR_VERSION 0)
set(QUADLODS_MINOR_VERSION 1)
set(QUADLODS_PATCH_VERSION 0)
set(QUADLODS_VERSION ${QUADLODS_MAJOR_VERSION}.${QUADLODS_MINOR_VERSION}.${QUADLODS_PATCH_VERSION})
set(COPY_YEAR 2018)

include_directories(${PROJECT_BINARY_DIR})
configure_file (config.h.in config.h)

set_target_properties(quadlib0 PROPERTIES OUTPUT_NAME "quadlods" VERSION ${QUADLODS_VERSION})
set_target_properties(quadlib1 PROPERTIES OUTPUT_NAME "quadlods" VERSION ${QUADLODS_VERSION})

set(CPACK_PACKAGE_VERSION_MAJOR ${QUADLODS_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${QUADLODS_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${QUADLODS_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)
