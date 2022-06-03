#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Setup the project
#
cmake_minimum_required(VERSION 2.8...3.23)

project(muon)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14")
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Setup include directory for this project


include_directories(${PROJECT_SOURCE_DIR}/include )

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Locate sources and headers for this project
#NB : headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)
#----------------------------------------------------------------------------
# Add the executable
add_executable(muon_rs main.cpp ${sources} ${headers})
#target_link_libraries(muon_rs   ${ROOT_LIBRARIES})
#target_link_libraries(main)
#----------------------------------------------------------------------------
add_custom_target(muon DEPENDS muon_rs)
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
file(COPY delay DESTINATION ./)
install(TARGETS muon_rs DESTINATION bin)
