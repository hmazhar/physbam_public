//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
using namespace PhysBAM;
using namespace std;
//#####################################################################
template<class T> OPENGL_MAC_VELOCITY_FIELD_2D<T>::
OPENGL_MAC_VELOCITY_FIELD_2D(GRID<TV> &grid, ARRAY<T,FACE_INDEX<2> > &face_velocities_input,T_ARRAYS_BOOL *active_cells_input,T_FACE_ARRAYS_BOOL *active_faces_input)
    : OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(vector_field,vector_locations),
     grid(grid),face_velocities(face_velocities_input),u(face_velocities.Component(1)),v(face_velocities.Component(2)),active_cells(active_cells_input),active_faces(active_faces_input)
{
    PHYSBAM_ASSERT(grid.MAC_offset == 0.5);
    Set_Velocity_Mode(CELL_CENTERED);
}

template<class T> OPENGL_MAC_VELOCITY_FIELD_2D<T>::
~OPENGL_MAC_VELOCITY_FIELD_2D()
{}

template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input)
{
    velocity_mode = velocity_mode_input;
    Update();
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0);
}

template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Update()
{
    PHYSBAM_ASSERT(grid.MAC_offset == 0.5);

    vector_field.Resize(0);
    vector_locations.Resize(0);

    if(u.counts.x-1!=v.counts.x || u.counts.y!=v.counts.y-1) return;

    int idx=1;
    if(velocity_mode == FACE_CENTERED){
        vector_field.Resize(u.counts.x*u.counts.y + v.counts.x*v.counts.y);
        vector_locations.Resize(u.counts.x*u.counts.y + v.counts.x*v.counts.y);
        for(int i=u.domain.min_corner.x;i<=u.domain.max_corner.x;i++) for(int j=u.domain.min_corner.y;j<=u.domain.max_corner.y;j++) if(!active_faces||(active_faces->Component(1))(i,j)){
            vector_field(idx)=VECTOR<T,2>(u(i,j),0);vector_locations(idx)=grid.X_Face(i,j);idx++;}
        for(int i=v.domain.min_corner.x;i<=v.domain.max_corner.x;i++) for(int j=v.domain.min_corner.y;j<=v.domain.max_corner.y;j++) if(!active_faces||(active_faces->Component(2))(i,j)){
            vector_field(idx)=VECTOR<T,2>(0,v(i,j));vector_locations(idx)=grid.Y_Face(i,j);idx++;}
        vector_field.Resize(idx-1);
        vector_locations.Resize(idx-1);}
    else{
        int number_of_cells=(u.counts.x-1)*(v.counts.y-1);
        vector_field.Resize(number_of_cells);
        vector_locations.Resize(number_of_cells);
        for(int i=u.domain.min_corner.x;i<=u.domain.max_corner.x-1;i++) for(int j=v.domain.min_corner.y;j<=v.domain.max_corner.y-1;j++) if(!active_cells||(*active_cells)(i,j)){
            vector_field(idx)=VECTOR<T,2>((T).5*(u(i,j)+u(i+1,j)),(T).5*(v(i,j)+v(i,j+1)));vector_locations(idx)=grid.Center(i,j);idx++;}
        vector_field.Resize(idx-1);
        vector_locations.Resize(idx-1);}
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_2D && grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)->index;
        T u_left,u_right,v_bottom,v_top;
        u_left=u(index.x,index.y);
        u_right=u(index.x+1,index.y);
        v_bottom=v(index.x,index.y);
        v_top=v(index.x,index.y+1);
        TV center_velocity(0.5*(u_left+u_right),0.5*(v_bottom+v_top));
        output_stream<<"    u left = "<<u_left<<",right = "<<u_right<<" avg="<<center_velocity.x<<std::endl;
        output_stream<<"    v bottom = "<<v_bottom<<",top = "<<v_top<<" avg="<<center_velocity.y<<std::endl;
        T ux=(u_right-u_left)*grid.one_over_dX.x,vy=(v_top-v_bottom)*grid.one_over_dX.y;
        output_stream<<"    divergence = "<<ux+vy<<" (ux="<<ux<<", vy="<<vy<<")"<<std::endl;}
    if(selection && selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *particle_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
        VECTOR<T,2> interp(interpolation.Clamped_To_Array(grid.Get_X_Face_Grid(),u,particle_selection->location),interpolation.Clamped_To_Array(grid.Get_Y_Face_Grid(),v,particle_selection->location));
        output_stream<<"    @ particle = "<<interp<<std::endl;}
}
//#####################################################################
// Convenience functions
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Velocity_Mode()
{
    VELOCITY_MODE new_velocity_mode = (VELOCITY_MODE)(((int)velocity_mode+1)%2);
    Set_Velocity_Mode(new_velocity_mode);
}
//#####################################################################
template class OPENGL_MAC_VELOCITY_FIELD_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_MAC_VELOCITY_FIELD_2D<double>;
#endif
