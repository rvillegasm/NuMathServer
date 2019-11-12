
cmake_minimum_required(VERSION 3.9.1)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "ERROR: In-source builds are not allowed!. 
  You should create a separate directory for build files.")
endif()

project(test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(SOURCES main.cpp
            external/fparser.cc
            external/fpoptimizer.cc)

find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(test "${SOURCES}")

target_include_directories(test PUBLIC "/usr/local/include/pistache")
target_link_libraries(test PUBLIC "/usr/local/lib/libpistache.a" pthread numath)
# target_link_libraries(test PRIVATE nlohmann_json::nlohmann_json)