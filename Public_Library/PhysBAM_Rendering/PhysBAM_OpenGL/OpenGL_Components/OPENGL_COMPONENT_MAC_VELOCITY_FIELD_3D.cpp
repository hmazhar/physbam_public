//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D(const GRID<TV> &grid, const std::string &velocity_filename)
    : OPENGL_COMPONENT("MAC Velocity Field"),
     opengl_mac_velocity_field(*new GRID<TV>(grid)),
     opengl_vorticity_magnitude(opengl_mac_velocity_field.grid,*(new ARRAY<T,VECTOR<int,3> >),OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1)),
     draw_vorticity(false),velocity_filename(velocity_filename),valid(false),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    is_animation = FILE_UTILITIES::Is_Animated(velocity_filename);
    opengl_vorticity_magnitude.Set_Scale_Range(0,100);
    frame_loaded = -1;
}

template<class T,class RW> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D()
{
    delete &opengl_mac_velocity_field.grid;
}

template<class T,class RW> bool OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(velocity_filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Display(const int in_color) const
{
    if (valid && draw){
        opengl_mac_velocity_field.Display(in_color);
        if(draw_vorticity) opengl_vorticity_magnitude.Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_mac_velocity_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Reinitialize()
{
    if (draw)
    {
        if (!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(velocity_filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_mac_velocity_field.face_velocities);//u,opengl_mac_velocity_field.v,opengl_mac_velocity_field.w);
            else
                return;

            opengl_mac_velocity_field.Update();

            if(draw_vorticity){
                Update_Vorticity();
                opengl_vorticity_magnitude.Update();}

            frame_loaded = frame;
            valid = true;
        }
    }
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": ";
        opengl_mac_velocity_field.Print_Selection_Info(stream,selection);
        if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_3D && opengl_mac_velocity_field.grid.Is_MAC_Grid()){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)selection)->index;
            if(draw_vorticity) stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude.values(index)<<std::endl;}}
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Set_Vector_Size(double size)
{
    opengl_mac_velocity_field.size = size;
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Velocity_Mode()
{
    opengl_mac_velocity_field.Toggle_Velocity_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Velocity_Mode_And_Draw()
{
    if (draw)
    {
        Toggle_Velocity_Mode();
        if ((int)opengl_mac_velocity_field.velocity_mode == 0) Toggle_Draw();
    }
    else Toggle_Draw();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Increase_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Decrease_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1/1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Arrowhead()
{
    opengl_mac_velocity_field.Toggle_Arrowhead_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Toggle_Draw_Vorticity()
{
    draw_vorticity=!draw_vorticity;
    if(draw_vorticity) valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude.Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T,RW>::
Update_Vorticity()
{
    GRID<TV>& grid=opengl_mac_velocity_field.grid;
    RANGE<VECTOR<int,3> > domain_indices(grid.Domain_Indices());domain_indices.Change_Size(-VECTOR<int,3>::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup(opengl_mac_velocity_field.face_velocities);
    opengl_vorticity_magnitude.values.Resize(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid,domain_indices);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        T vorticity_magnitude=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index).Magnitude();
        opengl_vorticity_magnitude.values(index)=vorticity_magnitude;
        min_vorticity=min(min_vorticity,vorticity_magnitude);
        max_vorticity=max(max_vorticity,vorticity_magnitude);}
}
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<double,double>;
#endif
