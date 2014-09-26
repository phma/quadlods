project(quadlods)
add_executable(quadlods quadlods.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(GMPXX REQUIRED)
set(LIBS ${LIBS} ${GMP_LIBRARY} ${GMPXX_LIBRARIES})
target_link_libraries(quadlods ${LIBS})
