if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -DNDEBUG -Wfatal-errors -Ofast -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -DNDEBUG -Wfatal-errors -O3 -mllvm -vectorize -ffast-math -march=native") 
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -DNDEBUG -Wfatal-errors -fast -axAVX")
endif()

# Release mode by default
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_CXX_FLAGS_DEBUG   "-std=c++11 -DBZ_DEBUG")

find_package(Boost COMPONENTS iostreams system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_path(BLITZ_INCLUDE_DIR NAMES blitz/array.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${BLITZ_INCLUDE_DIR}")

if(CMAKE_COMPILER_IS_GNUCXX)
  add_executable(test-cpp-gnucxx test.cpp)
else()
  add_executable(test-cpp-nongnu test.cpp)
endif()
add_executable(inout-cpp inout.cpp)

find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES NAMES gnuplot-iostream/gnuplot-iostream.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${GNUPLOT-IOSTREAM_INCLUDE_DIR}")
add_executable(plot-cpp plot.cpp)
target_link_libraries(plot-cpp ${Boost_LIBRARIES})

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  find_library(BLITZ_LIBRARY NAMES blitz)
  target_link_libraries(plot-cpp "blitz")
  target_link_libraries(test-cpp "blitz")
  target_link_libraries(inout-cpp "blitz")
endif()
