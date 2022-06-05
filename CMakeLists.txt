cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project( FEBioRemap )

# Set a default build type if none was specified
set(default_build_type "Release")

add_compile_options(-fPIC)

# FEBIO
if(NOT EXISTS ${FEBIO_LIB_DIR}/libfecore.a)
  message(SEND_ERROR "Could not find FEBio library (libfecore.a). Check FEBIO_LIB_DIR.")
  return()
endif()
ADD_LIBRARY(fecore STATIC IMPORTED)
SET_TARGET_PROPERTIES(fecore PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libfecore.a)

if(NOT EXISTS ${FEBIO_LIB_DIR}/libfebiomech.a)
  message(SEND_ERROR "Could not find FEBio library (libfebiomech.a). Check FEBIO_LIB_DIR.")
  return()
endif()
ADD_LIBRARY(febiomech STATIC IMPORTED)
SET_TARGET_PROPERTIES(febiomech PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libfebiomech.a)


# OPENMP
find_package(OpenMP)
if ( NOT OPENMP_FOUND )
  message(STATUS "This project requires OpenMP, and will not be compiled.")
  return()  
endif()
add_compile_options(-fopenmp)


# MMG
if(NOT EXISTS ${MMG_LIB_DIR})
  message(SEND_ERROR "Could not find MMG library. Check MMG_LIB.")
endif()
if(NOT EXISTS ${MMG_INC_DIR}/mmg/libmmg.h)
  message(SEND_ERROR "Could not find libmmg.h. Check MMG_INC.")
endif()
add_definitions(-DHAS_MMG)
ADD_LIBRARY(mmg STATIC IMPORTED)
SET_TARGET_PROPERTIES(mmg PROPERTIES IMPORTED_LOCATION ${MMG_LIB_DIR}/libmmg3d.a)


# The plugin library
file(GLOB HDR_FEBioRemap "src/*.h")
file(GLOB SRC_FEBioRemap "src/*.cpp")
add_library(FEBioRemap SHARED ${HDR_FEBioRemap} ${SRC_FEBioRemap})

target_link_libraries(FEBioRemap febiomech fecore)
target_include_directories(FEBioRemap PRIVATE ${FEBIO_ROOT} ${MMG_INC_DIR})