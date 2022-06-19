cmake_minimum_required(VERSION 3.21)
project(Test_Spanner)

set(CMAKE_CXX_STANDARD 17)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

include_directories(/opt/homebrew/Cellar)

link_directories(/opt/homebrew/lib/)
link_directories(/opt/homebrew/Cellar/)

set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -std=c++17 -fexceptions -Wall -Wextra -O3 -stdlib=libc++ -Wno-unused-parameter -Wno-macro-redefined")
set(CMAKE_BUILD_TYPE "Debug")

find_package(CGAL)

add_executable(Test_Spanner main.cpp Utilities.h GraphPrinter.h LatexPrinter.h PointGenerators.h StretchFactorExact.h StretchFactorExperimental.h DelaunayL2.h WSPD.h WSPD.cpp point_partition.h neighbor.h)

target_link_libraries(Test_Spanner gmp)
target_link_libraries(Test_Spanner mpfr)
target_link_libraries(Test_Spanner CGAL)
