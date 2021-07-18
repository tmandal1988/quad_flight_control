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

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Set name of the project
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE})

foreach(LIBRARY ${PROJECT_LIBRARIES})
    add_subdirectory("${PROJECT_LIBRARIES_DIR}/${LIBRARY}")
endforeach(LIBRARY)

target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBRARIES})
