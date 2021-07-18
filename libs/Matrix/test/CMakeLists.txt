cmake_minimum_required(VERSION 3.16.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TEST_PROJECT_NAME
    MatrixLibTests
)
set(LIBRARY_TESTS_SOURCE
    matrix_lib_unittests.cc
)

project(${TEST_PROJECT_NAME})

add_executable(${TEST_PROJECT_NAME} ${LIBRARY_TESTS_SOURCE})
target_link_libraries(${TEST_PROJECT_NAME} 
    gtest_main
    ${LIBRARY_NAME} # NOTE: This is defined from project above
)
add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})