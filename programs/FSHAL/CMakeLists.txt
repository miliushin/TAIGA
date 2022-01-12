#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(shall)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Setup include directory for this project


include_directories(${PROJECT_SOURCE_DIR}/include /home/ilushin/programs/root/include)

set( ROOTPREFIX /home/ilushin/programs/root/lib/lib)
set( ROOT_LIBRARIES
   ${ROOTPREFIX}Gpad.so 
   ${ROOTPREFIX}Hist.so 
   ${ROOTPREFIX}Graf.so 
   ${ROOTPREFIX}Graf3d.so 
   ${ROOTPREFIX}Tree.so  
   ${ROOTPREFIX}Rint.so 
   ${ROOTPREFIX}Postscript.so  
   ${ROOTPREFIX}Matrix.so  
   ${ROOTPREFIX}Physics.so  
   ${ROOTPREFIX}MathCore.so 
   ${ROOTPREFIX}RIO.so 
   ${ROOTPREFIX}Net.so 
   ${ROOTPREFIX}Thread.so 
   ${ROOTPREFIX}Core.so 
   ) 


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(FSHAL main.cpp ${sources} ${headers})
target_link_libraries(FSHAL   ${ROOT_LIBRARIES})
#target_link_libraries(main)
#----------------------------------------------------------------------------
add_custom_target(shall DEPENDS FSHAL)
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FSHAL DESTINATION bin)
