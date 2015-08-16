cmake_minimum_required(VERSION 2.8)

project(PhysBAM)

#set path for local find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

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


option(PHYSBAM_USE_DOUBLE   		"Turn ON to enable double support"  		OFF )
option(PHYSBAM_USE_DYADIC   		"Turn ON to enable dyadic support"  		OFF )
option(PHYSBAM_USE_BINTREE  		"Turn ON to enable bintree support" 		OFF )
option(PHYSBAM_USE_RLE      		"Turn ON to enable RLE support"     		OFF )
option(PHYSBAM_USE_ZLIB     		"Turn ON to enable zlib support"    		OFF )
option(PHYSBAM_USE_ID_TYPES   		"Treat ID types as int to avoid possible performance consequences"    		OFF )
option(PHYSBAM_USE_FFMPEG   		"Turn ON to enable ffmpeg support"  		OFF )
option(PHYSBAM_USE_JPEG     		"Turn ON to enable libjpeg support" 		OFF )
option(PHYSBAM_USE_PNG      		"Turn ON to enable png support"     		OFF )
option(PHYSBAM_USE_OPENEXR  		"Turn ON to enable OpenEXR support" 		OFF )
option(PHYSBAM_USE_FFTW     		"Turn ON to enable FFTW support"    		OFF )
option(PHYSBAM_USE_MPI      		"Turn ON to enable MPI support"     		OFF )
option(PHYSBAM_USE_PTHREADS      	"Turn ON to enable PThreads support"     	OFF )
option(PHYSBAM_USE_COMPRESSION      "Turn ON to enable Compression support"     OFF )
option(PHYSBAM_BUILD_PROJECTS       "Turn ON to build example projects"         OFF )



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
	ADD_DEFINITIONS("-DUSE_LIBZ")
ENDIF()

IF (PHYSBAM_USE_ID_TYPES)			
	ADD_DEFINITIONS("-DCOMPILE_ID_TYPES_AS_INT")	
ENDIF()


IF (PHYSBAM_USE_FFMPEG)
	FIND_PACKAGE(FFMPEG)
	ADD_DEFINITIONS("-DUSE_FFMPEG")
ENDIF()

IF (PHYSBAM_USE_JPEG)
	FIND_PACKAGE(JPEG)
	ADD_DEFINITIONS("-DUSE_LIBJPEG")
ENDIF()

IF (PHYSBAM_USE_PNG)
	FIND_PACKAGE(PNG)
	ADD_DEFINITIONS("-DUSE_LIBPNG")
ENDIF()

IF (PHYSBAM_USE_OPENEXR)
	FIND_PACKAGE(OPENEXR)
	ADD_DEFINITIONS("-DUSE_USE_OPENEXR")
ENDIF()

IF (PHYSBAM_USE_FFTW)
	FIND_PACKAGE(FFTW)
	ADD_DEFINITIONS("-DUSE_FFTW")
ENDIF()

IF (PHYSBAM_USE_MPI)
	FIND_PACKAGE(MPI)
	ADD_DEFINITIONS("-DUSE_MPI")
ENDIF()

IF (PHYSBAM_USE_PTHREADS)
	ADD_DEFINITIONS("-DUSE_PTHREADS")
	find_package (Threads)
ENDIF()

IF (PHYSBAM_USE_COMPRESSION)
	ADD_DEFINITIONS("-DUSE_USE_COMPRESSION")
ENDIF()

SET (PHYSBAM_LIBS 
	${ZLIB_LIBRARY} 
	${JPEG_LIBRARY} 
	${PNG_LIBRARIES} 
	${FFMPEG_LIBRARIES} 
	${OPENEXR_LIBRARIES}
	${FFTW_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
)


ADD_SUBDIRECTORY(PhysBAM_Tools)
ADD_SUBDIRECTORY(PhysBAM_Geometry)
ADD_SUBDIRECTORY(PhysBAM_Rendering)

IF (PHYSBAM_BUILD_PROJECTS)
	ADD_SUBDIRECTORY(Projects)
ENDIF()

file(COPY ${PhysBAM_Tools_DIR} ${PhysBAM_Geometry_DIR} ${PhysBAM_OpenGL_DIR}
		DESTINATION ${PhysBAM_BINARY_DIR}/include
		FILES_MATCHING PATTERN "*.h" )
