cmake_minimum_required(VERSION 3.16.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBRARY_NAME
    Navio
)

add_library(${LIBRARY_NAME} INTERFACE)

set(NAVIO_LIBRARIES 
    Common
    Navio2
)

foreach(NAVIO_LIBRARY ${NAVIO_LIBRARIES})
    add_subdirectory("${NAVIO_LIBRARY}")
endforeach(NAVIO_LIBRARY)