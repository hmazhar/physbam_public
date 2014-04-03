cmake_minimum_required(VERSION 2.8)

project(PhysBAM)

# For non-multi-configuration generators (eg, make, Eclipse)
# The Visual Studio and XCode generators create a single project with all these
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "For single-configuration generators (e.g. make) set the type of build: Release, Debug, RelWithDebugInfo, MinSizeRel")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" "RelWithDebugInfo" "MinSizeRel")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PhysBAM_BINARY_DIR}/lib")
# PhysBAM has issues when build as shared libs. Keep off to build static.
option(BUILD_SHARED_LIBS "Build as shared libraries by default" OFF)

if (CMAKE_COMPILER_IS_GNUCC)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

set(PhysBAM_PL_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Path to the PhysBAM Public_Library directory")







option(PHYSBAM_USE_DOUBLE   "Turn ON to enable double support"  OFF )
option(PHYSBAM_USE_DYADIC   "Turn ON to enable dyadic support"  OFF )
option(PHYSBAM_USE_BINTREE  "Turn ON to enable bintree support" ON  )
option(PHYSBAM_USE_RLE      "Turn ON to enable RLE support"     OFF )
option(PHYSBAM_USE_ZLIB     "Turn ON to enable zlib support"    OFF )
option(PHYSBAM_USE_FFMPEG   "Turn ON to enable ffmpeg support"  OFF )
option(PHYSBAM_USE_JPEG     "Turn ON to enable libjpeg support" OFF )
option(PHYSBAM_USE_PNG      "Turn ON to enable png support"     OFF )
option(PHYSBAM_USE_OPENEXR  "Turn ON to enable OpenEXR support" OFF )
option(PHYSBAM_USE_FFTW     "Turn ON to enable FFTW support"    OFF )
option(PHYSBAM_USE_MPI      "Turn ON to enable MPI support"     OFF )

IF (NOT PHYSBAM_USE_DOUBLE_SUPPORT)
		ADD_DEFINITIONS("-DCOMPILE_WITHOUT_DOUBLE_SUPPORT")
ENDIF()
IF (NOT PHYSBAM_USE_DYADIC_SUPPORT)
		ADD_DEFINITIONS("-DCOMPILE_WITHOUT_DYADIC_SUPPORT")
ENDIF()
IF (PHYSBAM_USE_BINTREE_SUPPORT)
		ADD_DEFINITIONS("-DCOMPILE_WITH_BINTREE_SUPPORT")
ENDIF()
IF (NOT PHYSBAM_USE_RLE)
		ADD_DEFINITIONS("-DCOMPILE_WITHOUT_RLE_SUPPORT")
ENDIF()

IF (NOT PHYSBAM_USE_ZLIB)
		ADD_DEFINITIONS("-DCOMPILE_WITHOUT_ZLIB_SUPPORT")
ELSEIF()
		FIND_PACKAGE(ZLIB)
ENDIF()
IF (PHYSBAM_USE_FFMPEG)
		FIND_PACKAGE(FFMPEG)
ENDIF()
IF (PHYSBAM_USE_JPEG)
		FIND_PACKAGE(JPEG)
ENDIF()
IF (PHYSBAM_USE_PNG)
		FIND_PACKAGE(PNG)
ENDIF()
IF (PHYSBAM_USE_OPENEXR)
ENDIF()
IF (PHYSBAM_USE_FFTW)
ENDIF()
IF (PHYSBAM_USE_MPI)
		FIND_PACKAGE(MPI)
ENDIF()


ADD_SUBDIRECTORY(PhysBAM_Tools)
ADD_SUBDIRECTORY(PhysBAM_Geometry)
ADD_SUBDIRECTORY(PhysBAM_Rendering)


file(COPY ${PhysBAM_Tools_DIR} ${PhysBAM_Geometry_DIR} ${PhysBAM_OpenGL_DIR}
		DESTINATION ${PhysBAM_BINARY_DIR}/include
		FILES_MATCHING PATTERN "*.h" )
