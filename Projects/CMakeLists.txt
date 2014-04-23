
include_directories(${PhysBAM_PL_DIR})

add_executable(smoke 
${PhysBAM_PL_DIR}/Projects/user_projects/Smoke/SMOKE_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/Smoke/SMOKE_EXAMPLE.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/Smoke/main.cpp
)

target_link_libraries(smoke PhysBAM_Tools PhysBAM_Geometry)

add_executable(levelset 
${PhysBAM_PL_DIR}/Projects/user_projects/Levelset/WATER_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/Levelset/WATER_EXAMPLE.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/Levelset/main.cpp
)

target_link_libraries(levelset PhysBAM_Tools PhysBAM_Geometry)


add_executable(opengl_1d  ${PhysBAM_PL_DIR}/Projects/opengl_1d/main.cpp )
target_link_libraries(opengl_1d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(opengl_2d  ${PhysBAM_PL_DIR}/Projects/opengl_2d/main.cpp )
target_link_libraries(opengl_2d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(opengl_3d  ${PhysBAM_PL_DIR}/Projects/opengl_3d/main.cpp )
target_link_libraries(opengl_3d PhysBAM_Tools PhysBAM_Geometry PhysBAM_OpenGL)


add_executable(strand_forward_euler 
${PhysBAM_PL_DIR}/Projects/user_projects/strand_forward_euler/SIMULATION_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/strand_forward_euler/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/strand_forward_euler/main.cpp
)

target_link_libraries(strand_forward_euler PhysBAM_Tools PhysBAM_Geometry)


add_executable(strand_semi_implicit
${PhysBAM_PL_DIR}/Projects/user_projects/strand_semi_implicit/SIMULATION_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/strand_semi_implicit/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/strand_semi_implicit/main.cpp
)

target_link_libraries(strand_semi_implicit PhysBAM_Tools PhysBAM_Geometry)

# ball_semi_implicit
add_executable(ball_semi_implicit
${PhysBAM_PL_DIR}/Projects/user_projects/ball_semi_implicit/SIMULATION_DRIVER.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/ball_semi_implicit/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/ball_semi_implicit/main.cpp
)

target_link_libraries(ball_semi_implicit PhysBAM_Tools PhysBAM_Geometry)

# collision_backward_euler
add_executable(collision_backward_euler
${PhysBAM_PL_DIR}/Projects/user_projects/collision_backward_euler/main.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/collision_backward_euler/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/collision_backward_euler/SIMULATION_DRIVER.cpp
)

target_link_libraries(collision_backward_euler PhysBAM_Tools PhysBAM_Geometry)

# collision_semi_implicit
add_executable(collision_semi_implicit
${PhysBAM_PL_DIR}/Projects/user_projects/collision_semi_implicit/main.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/collision_semi_implicit/SIMULATION_LAYOUT.cpp
${PhysBAM_PL_DIR}/Projects/user_projects/collision_semi_implicit/SIMULATION_DRIVER.cpp
)

target_link_libraries(collision_semi_implicit PhysBAM_Tools PhysBAM_Geometry)
