#--------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(General)

#--------------------------------------------------------------------------
# Find Geant4 package, activating all UI and Vis drivers by default
# You can set WITH_GEATH4_UIVIS to OFF via the command line or 
# ccmake/cmake-#gui to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Building with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()
#--------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#--------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#--------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(run main.cc ${sources} ${headers}) 
target_link_libraries(run ${Geant4_LIBRARIES})

#--------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build General.  This is so that we can run the executable directory 
# because it relies on these scripts being in the current working directory.
#

set(GENERAL_SCRIPTS
  vis.mac
  input
	)

foreach(_script ${GENERAL_SCRIPTS})
  configure_file(
  ${PROJECT_SOURCE_DIR}/${_script}
  ${PROJECT_BINARY_DIR}/${_script}
  COPYONLY
  )
endforeach()

#--------------------------------------------------------------------------
# Install the executable to the 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS run DESTINATION bin)

