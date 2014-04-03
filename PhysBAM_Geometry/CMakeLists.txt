set(PhysBAM_Geometry_DIR "${PhysBAM_PL_DIR}/PhysBAM_Geometry" CACHE PATH "Path to the PhysBAM Public_Library/PhysBAM_Geometry")

file(GLOB_RECURSE PhysBAM_Geometry_SOURCE "${PhysBAM_Geometry_DIR}/*.cpp")
file(GLOB_RECURSE PhysBAM_Geometry_HEADER "${PhysBAM_Geometry_DIR}/*.h")

include_directories(${PhysBAM_PL_DIR})

add_library(PhysBAM_Geometry ${PhysBAM_Geometry_SOURCE} ${PhysBAM_Geometry_HEADER})

target_link_libraries(PhysBAM_Geometry PhysBAM_Tools)