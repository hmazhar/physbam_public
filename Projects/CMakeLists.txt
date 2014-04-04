
include_directories(${PhysBAM_PL_DIR})

add_executable(smoke 
${PhysBAM_PL_DIR}/Projects/Smoke/SMOKE_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/Smoke/SMOKE_EXAMPLE.cpp
${PhysBAM_PL_DIR}/Projects/Smoke/main.cpp
)

target_link_libraries(smoke PhysBAM_Tools PhysBAM_Geometry)

add_executable(levelset 
${PhysBAM_PL_DIR}/Projects/Levelset/WATER_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/Levelset/WATER_EXAMPLE.cpp
${PhysBAM_PL_DIR}/Projects/Levelset/main.cpp
)

target_link_libraries(levelset PhysBAM_Tools PhysBAM_Geometry)


add_executable(opengl_1d  ${PhysBAM_PL_DIR}/Projects/opengl_1d/main.cpp )
target_link_libraries(opengl_1d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(opengl_2d  ${PhysBAM_PL_DIR}/Projects/opengl_2d/main.cpp )
target_link_libraries(opengl_2d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(opengl_3d  ${PhysBAM_PL_DIR}/Projects/opengl_3d/main.cpp )
target_link_libraries(opengl_3d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(strand_forward_euler 
${PhysBAM_PL_DIR}/Projects/strand_forward_euler/SIMULATION_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/strand_forward_euler/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/strand_forward_euler/main.cpp
)

target_link_libraries(strand_forward_euler PhysBAM_Tools PhysBAM_Geometry)


add_executable(strand_semi_implicit
${PhysBAM_PL_DIR}/Projects/strand_semi_implicit/SIMULATION_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/strand_semi_implicit/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/strand_semi_implicit/main.cpp
)

target_link_libraries(strand_semi_implicit PhysBAM_Tools PhysBAM_Geometry)
