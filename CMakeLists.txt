# Minimum version of CMake required to build this project
cmake_minimum_required(VERSION 3.16.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME MatrixMain)

set(PROJECT_SOURCE_DIR
	src
)

set(PROJECT_SOURCE
    ${PROJECT_SOURCE_DIR}/matrix_main.cc
)

set(PROJECT_LIBRARIES_DIR
    libs
)

set(PROJECT_LIBRARIES 
    Matrix
)

# Set name of the project
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE})

foreach(LIBRARY ${PROJECT_LIBRARIES})
    add_subdirectory("${PROJECT_LIBRARIES_DIR}/${LIBRARY}")
endforeach(LIBRARY)

target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBRARIES})
