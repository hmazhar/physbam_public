//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
using namespace PhysBAM;
//#####################################################################
template<class T> OPENGL_MAC_VELOCITY_FIELD_3D<T>::
OPENGL_MAC_VELOCITY_FIELD_3D(GRID<TV> &grid)
    :OPENGL_VECTOR_FIELD_3D<T>(*(new ARRAY<TV>),*(new ARRAY<TV>),OPENGL_COLOR::Gray(0.8f),.25,false,false,true),scale(1),grid(grid),
    face_velocities(*new ARRAY<T,FACE_INDEX<3> >),u(face_velocities.Component(1)),v(face_velocities.Component(2)),w(face_velocities.Component(3))
{
    max_vectors_3d = 10000;
    PHYSBAM_ASSERT(grid.MAC_offset == 0.5);
    Set_Velocity_Mode(CELL_CENTERED);
}

template<class T> OPENGL_MAC_VELOCITY_FIELD_3D<T>::
~OPENGL_MAC_VELOCITY_FIELD_3D()
{
    delete &vector_field;delete &vector_locations;
}

template<class T> void OPENGL_MAC_VELOCITY_FIELD_3D<T>::
Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input)
{
    velocity_mode = velocity_mode_input;
    Update();
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_MAC_VELOCITY_FIELD_3D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain);
}

template<class T> void OPENGL_MAC_VELOCITY_FIELD_3D<T>::
Update()
{
    PHYSBAM_ASSERT(grid.MAC_offset == 0.5);

    vector_field.Resize(0);
    vector_locations.Resize(0);

    if(u.counts+TV_INT(-1,1,0)!=v.counts) return; // not a complete check that they match but good enough

    // Make sure the slice is good
    VECTOR<int,3> domain_start(1,1,1),domain_end(grid.counts);
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    if(slice && ((slice->mode == OPENGL_SLICE::CELL_SLICE && ((slice->index-1)/scale+1 < domain_start[slice->axis] || (slice->index-1)/scale+1 > domain_end[slice->axis])) || (slice->mode == OPENGL_SLICE::NODE_SLICE)))
        return; // Currently we don't draw anything if the slice doesn't match where the velocity field lives

    VECTOR<int,3> cell_start(face_velocities.Component(1).domain.min_corner),cell_end(face_velocities.Component(1).domain.max_corner-1);
    int m=cell_end.x-cell_start.x+1,n=cell_end.y-cell_start.y+1,mn=cell_end.z-cell_start.z+1;

    if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE)
        switch (slice->axis)
        {
            case 1: cell_start.x=cell_end.x=(slice->index-1)/scale+1;m=1; break;
            case 2: cell_start.y=cell_end.y=(slice->index-1)/scale+1;n=1; break;
            case 3: cell_start.z=cell_end.z=(slice->index-1)/scale+1;mn=1; break;
        }
    
    if(velocity_mode == FACE_CENTERED){
        int num_vectors=3*m*n*mn+n*mn+m*mn+m*n;
        vector_field.Resize(num_vectors);
        vector_locations.Resize(num_vectors);

        int idx=0;

        // u velocities
        for(int i=cell_start.x;i<=cell_end.x+1;i++)for(int j=cell_start.y;j<=cell_end.y;j++)for(int k=cell_start.z;k<=cell_end.z;k++){
            T vel=u(VECTOR<int,3>(i,j,k));
            if(vel != 0){idx++;vector_field(idx)=VECTOR<T,3>(vel,0,0);vector_locations(idx)=grid.X_Face(i,j,k);}}

        // v velocities
        for (int i=cell_start.x;i<=cell_end.x;i++) for (int j=cell_start.y;j<=cell_end.y+1;j++) for (int k=cell_start.z;k<=cell_end.z;k++){
            T vel = v(VECTOR<int,3>(i,j,k));
            if (vel != 0){idx++;vector_field(idx)=VECTOR<T,3>(0,vel,0);vector_locations(idx)=grid.Y_Face(i,j,k);}}

        // w velocities
        for (int i=cell_start.x;i<=cell_end.x;i++) for (int j=cell_start.y;j<=cell_end.y;j++) for (int k=cell_start.z;k<=cell_end.z+1;k++){
            T vel = w(VECTOR<int,3>(i,j,k));
            if (vel != 0){idx++;vector_field(idx)=VECTOR<T,3>(0,0,vel);vector_locations(idx)=grid.Z_Face(i,j,k);}}

        vector_field.Resize(idx);
        vector_locations.Resize(idx);
    }
    else
    {
        int inc=1;
        if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE){
            int num_vectors=m*n*mn;
            inc=(int)pow((double)PhysBAM::max(1, num_vectors/max_vectors_3d),(double)one_third);
            m=(int)ceil((T)m/inc);n=(int)ceil((T)n/inc);mn=(int)ceil((T)mn/inc);}

        int num_vectors=m*n*mn;
        vector_field.Resize(num_vectors);
        vector_locations.Resize(num_vectors);

        int idx=0;VECTOR<int,3> index;
        for(index.x=cell_start.x;index.x<=cell_end.x;index.x+=inc)for(index.y=cell_start.y;index.y<=cell_end.y;index.y+=inc)for(index.z=cell_start.z;index.z<=cell_end.z;index.z+=inc){
            idx++;
            vector_field(idx) = (T)0.5*VECTOR<T,3>(u(index)+u(index+VECTOR<int,3>(1,0,0)),
                                                    v(index)+v(index+VECTOR<int,3>(0,1,0)),
                                                    w(index)+w(index+VECTOR<int,3>(0,0,1)));
            vector_locations(idx) = grid.X(index);}
    }
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_3D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_3D && grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)selection)->index;
        T u_left=u(index.x,index.y,index.z),u_right=u(index.x+1,index.y,index.z),
            v_bottom=v(index.x,index.y,index.z),v_top=v(index.x,index.y+1,index.z),
            w_back=w(index.x,index.y,index.z),w_front=w(index.x,index.y,index.z+1);

        TV center_velocity(0.5*(u_left+u_right),0.5*(v_bottom+v_top),0.5*(w_back+w_front));
        output_stream<<"    u left = "<<u_left<<",right = "<<u_right<<" avg="<<center_velocity.x<<std::endl;
        output_stream<<"    v bottom = "<<v_bottom<<",top = "<<v_top<<" avg="<<center_velocity.y<<std::endl;
        output_stream<<"    w back = "<<w_back<<",front = "<<w_front<<" avg="<<center_velocity.z<<std::endl;
        T ux=(u_right-u_left)*grid.one_over_dX.x,vy=(v_top-v_bottom)*grid.one_over_dX.y,wz=(w_front-w_back)*grid.one_over_dX.z;
        output_stream<<"    divergence = "<<ux+vy+wz<<" (ux="<<ux<<", vy="<<vy<<", wz="<<wz<<")"<<std::endl;}
    if(selection && selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_3D){
        OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T> *particle_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>*)selection;
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
        TV interp(interpolation.Clamped_To_Array(grid.Get_X_Face_Grid(),u,particle_selection->location),
            interpolation.Clamped_To_Array(grid.Get_Y_Face_Grid(),v,particle_selection->location),
            interpolation.Clamped_To_Array(grid.Get_Z_Face_Grid(),w,particle_selection->location));
        output_stream<<"    @ particle = "<<interp<<std::endl;}
}
//#####################################################################
// Convenience functions
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_3D<T>::
Toggle_Velocity_Mode()
{
    VELOCITY_MODE new_velocity_mode = (VELOCITY_MODE)(((int)velocity_mode+1)%2);
    Set_Velocity_Mode(new_velocity_mode);
}
//#####################################################################
template class OPENGL_MAC_VELOCITY_FIELD_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_MAC_VELOCITY_FIELD_3D<double>;
#endif
