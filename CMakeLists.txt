project(irafl)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s -Os -static-libgcc -static-libstdc++")

set(IRAFL_MAJOR_VERSION 1)
set(IRAFL_MINOR_VERSION 0)
set(IRAFL_PATCH_VERSION 0)
set(IRAFL_VERSION ${IRAFL_MAJOR_VERSION}.${IRAFL_MINOR_VERSION}.${IRAFL_PATCH_VERSION})

file(
    GLOB srcs
    _main.cpp
    src/irafl.cpp
    src/functions/*.cpp
)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/../glm
)

add_library(${PROJECT_NAME} SHARED ${srcs})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${IRAFL_VERSION} SOVERSION ${IRAFL_MAJOR_VERSION})

install(DIRECTORY include DESTINATION ".")
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
