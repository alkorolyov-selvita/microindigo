cmake_minimum_required(VERSION 3.6)
project(main LANGUAGES CXX)

add_executable(${PROJECT_NAME} main.cpp
        main.h)

target_link_libraries(${PROJECT_NAME} indigo-core object_threadsafe)
target_include_directories(${PROJECT_NAME}
        PRIVATE ${CMAKE_SOURCE_DIR}/indigo-core
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
