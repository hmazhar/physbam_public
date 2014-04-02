#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D.h>
namespace PhysBAM{

template<class T,class T2,class RW> OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D(RLE_GRID_3D<T> &grid,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT("RLE Based Scalar Field"),opengl_rle_cell_scalar_field(grid,*(new ARRAY<T2>),color_map_input),opengl_triangulated_surface(0),triangulated_surface(0),
    scalar_field_filename(scalar_field_filename_input),valid(false),draw_surface(false)
{
    is_animation=FILE_UTILITIES::Is_Animated(scalar_field_filename);
    frame_loaded=-1;
}

template<class T,class T2,class RW> OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
~OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D()
{
}

template<class T,class T2,class RW> bool OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T> static void Display_Surface(const OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T,T>& component)
{
    if(!component.opengl_triangulated_surface){
        PHYSBAM_ASSERT(!component.triangulated_surface);
        LEVELSET_RLE<RLE_GRID_3D<T> > levelset(component.opengl_rle_cell_scalar_field.grid,component.opengl_rle_cell_scalar_field.value);
        component.triangulated_surface=DUALCONTOUR_RLE_3D<T>::Create_Triangulated_Surface_From_Levelset(levelset);
        component.opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*component.triangulated_surface,true,
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(.8f,0,0)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR(0,0,0.8)));
        *component.opengl_triangulated_surface->vertex_normals=*component.triangulated_surface->vertex_normals;}
    component.opengl_triangulated_surface->Display();
}
    
template<> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<float,float,float>::
Display(const int in_color) const
{
    if(valid && draw){
        if(!draw_surface) opengl_rle_cell_scalar_field.Display(in_color);
        else Display_Surface(*this);}
}

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<double,double,double>::
Display(const int in_color) const
{
    if(valid && draw){
        if(!draw_surface) opengl_rle_cell_scalar_field.Display(in_color);
        else Display_Surface(*this);}
}
#endif

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_rle_cell_scalar_field.Display(in_color);
}

template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_rle_cell_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_rle_cell_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Reinitialize()
{
    if(draw && (!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))){
        valid=false;
        delete opengl_triangulated_surface;opengl_triangulated_surface=0;
        delete triangulated_surface;triangulated_surface=0;
        std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename,frame);
        if(!FILE_UTILITIES::File_Exists(tmp_filename)) return;
        FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_rle_cell_scalar_field.value);
        opengl_rle_cell_scalar_field.Update();
        frame_loaded=frame;
        valid=true;}
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Increase_Point_Size()
{
    opengl_rle_cell_scalar_field.point_size++;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Decrease_Point_Size()
{
    opengl_rle_cell_scalar_field.point_size--;
}

template<class T,class T2,class RW> void OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,T2,RW>::
Draw_Surface()
{
    draw_surface=!draw_surface;
}

template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<float,int,float>;
template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<float,bool,float>;
template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<double,int,double>;
template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<double,bool,double>;
template class OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<double,double,double>;
#endif
}
#endif
