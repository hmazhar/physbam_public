set(PhysBAM_OpenGL_DIR "${PhysBAM_PL_DIR}/PhysBAM_Rendering/PhysBAM_OpenGL" CACHE PATH "Path to the PhysBAM Public_Library/PhysBAM_OpenGL")
set(PhysBAM_RayTracing_DIR "${PhysBAM_PL_DIR}/PhysBAM_Rendering/PhysBAM_Ray_Tracing" CACHE PATH "Path to the PhysBAM Public_Library/PhysBAM_Ray_Tracing")


file(GLOB_RECURSE PhysBAM_OpenGL_SOURCE "${PhysBAM_OpenGL_DIR}/*.cpp")
file(GLOB_RECURSE PhysBAM_OpenGL_HEADER "${PhysBAM_OpenGL_DIR}/*.h")

file(GLOB_RECURSE PhysBAM_Ray_Tracing_SOURCE "${PhysBAM_RayTracing_DIR}/*.cpp")
file(GLOB_RECURSE PhysBAM_Ray_Tracing_HEADER "${PhysBAM_RayTracing_DIR}/*.h")


include_directories(${PhysBAM_PL_DIR})

add_library(PhysBAM_OpenGL ${PhysBAM_OpenGL_SOURCE} ${PhysBAM_OpenGL_HEADER})
add_library(PhysBAM_RayTracing ${PhysBAM_Ray_Tracing_SOURCE} ${PhysBAM_Ray_Tracing_HEADER})

FIND_PACKAGE(OpenGL)
FIND_PACKAGE(GLUT)

target_link_libraries(PhysBAM_OpenGL PhysBAM_Tools PhysBAM_Geometry ${PHYSBAM_LIBS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

target_link_libraries(PhysBAM_RayTracing PhysBAM_Tools PhysBAM_Geometry ${PHYSBAM_LIBS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})