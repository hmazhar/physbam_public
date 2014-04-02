#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_GRID.h>
using namespace PhysBAM;

//#####################################################################
// OPENGL_COMPONENT_OCTREE_GRID
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
OPENGL_COMPONENT_OCTREE_GRID(const std::string& filename_input):OPENGL_COMPONENT("Octree Area"),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation=true;//FILE_UTILITIES::Is_Animated(filename);
    Reinitialize();
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Display(const int in_color) const
{
    if(valid&&draw)opengl_grid.Display(in_color);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Bounding_Box() const
{
    if(valid&&draw)return opengl_grid.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    opengl_grid.Print_Selection_Info(output_stream,selection);
}
//#####################################################################
// Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    return opengl_grid.Create_Or_Destroy_Selection_After_Frame_Change(old_selection,delete_selection);
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE_GRID<T,RW>::
Reinitialize(bool force)
{
    if (force || ((is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename_string=FILE_UTILITIES::Get_Frame_Filename(filename, frame);
        if (FILE_UTILITIES::File_Exists(filename_string)){
            FILE_UTILITIES::Read_From_File<RW>(filename_string,opengl_grid.grid);
            opengl_grid.grid.Node_Iterator_Data();opengl_grid.grid.Face_Iterator_Data();}
        else return;
        opengl_grid.Update();
        frame_loaded=frame;valid=true;}
}
template class OPENGL_COMPONENT_OCTREE_GRID<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_OCTREE_GRID<double,double>;
#endif
#endif
