CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

INCLUDE(CMake.in.config)

FIND_PACKAGE(OpenFLUIDHelpers REQUIRED)
 
find_package( Boost REQUIRED COMPONENTS program_options regex )
include_directories( ${Boost_INCLUDE_DIRS} )

OPENFLUID_ADD_SIMULATOR(SIM)
