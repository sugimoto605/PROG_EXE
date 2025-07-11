cmake_minimum_required(VERSION 3.10)
set(CMAKE_OSX_SYSROOT macosx)
project(MultiTargets)

# Build-Modeの設定
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Library Search
cmake_policy(SET CMP0167 NEW)
find_package(Boost REQUIRED COMPONENTS timer)

# P1ターゲット
add_executable(P1 P1/main.cpp)

# P2ターゲット
add_executable(P2 P2/main.cpp)

# P3ターゲット
add_executable(P3 P3/main.cpp)

# P4ターゲット
add_executable(P4 P4/main.cpp)

# BOOSTターゲット
add_executable(BOOST BOOST/main.cpp)
target_link_libraries(BOOST PRIVATE Boost::timer)