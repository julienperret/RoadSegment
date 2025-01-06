# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(RoadSegment)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)
find_package(GDAL CONFIG REQUIRED)

create_single_source_cgal_program("sdg.cpp")
target_link_libraries(sdg PRIVATE GDAL::GDAL)
