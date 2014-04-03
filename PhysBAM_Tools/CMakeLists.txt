set(PhysBAM_Tools_DIR "${PhysBAM_PL_DIR}/PhysBAM_Tools" CACHE PATH "Path to the PhysBAM Public_Library/PhysBAM_Tools")

file(GLOB_RECURSE PhysBAM_Tools_SOURCE "${PhysBAM_Tools_DIR}/*.cpp")
file(GLOB_RECURSE PhysBAM_Tools_HEADER "${PhysBAM_Tools_DIR}/*.h")

include_directories(${PhysBAM_PL_DIR})

add_library(PhysBAM_Tools ${PhysBAM_Tools_SOURCE} ${PhysBAM_Tools_HEADER})

