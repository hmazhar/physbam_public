#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD.h>

namespace PhysBAM{
//#####################################################################
// Function OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD(OCTREE_GRID<T> &grid,const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
     :OPENGL_COMPONENT("Octree Cell Scalar Field"), opengl_octree_cell_scalar_field(grid,*(new ARRAY<T2>),color_map_input), 
      opengl_triangulated_surface(0),triangulated_surface(0),
      scalar_field_filename(scalar_field_filename_input),valid(false),draw_surface(false)
{
    is_animation=FILE_UTILITIES::Is_Animated(scalar_field_filename);frame_loaded=-1;
    Reinitialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
~OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD()
{}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
    if (draw_input) opengl_octree_cell_scalar_field.Set_Slice(slice);
}
//#####################################################################
// Function Display
//#####################################################################
template<> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<float,float,float>::
Display(const int in_color) const
{
    if(valid && draw){
        if(!draw_surface) opengl_octree_cell_scalar_field.Display(in_color);
        else{
            if(!opengl_triangulated_surface){
                OCTREE_GRID<float>& grid=opengl_octree_cell_scalar_field.grid;ARRAY<float>& phi=opengl_octree_cell_scalar_field.value;
                LEVELSET_OCTREE<float> levelset(grid,phi);
                DUALCONTOUR_OCTREE<float> contour(&levelset);
                if(triangulated_surface){delete &triangulated_surface->mesh;delete &triangulated_surface->particles;}
                delete triangulated_surface;delete opengl_triangulated_surface;
                triangulated_surface=contour.Get_Triangulated_Surface();
                opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<float>(*triangulated_surface,true,
                    OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0.8,(float)0,0)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0,0,(float)0.8)));
                *opengl_triangulated_surface->vertex_normals=*triangulated_surface->vertex_normals;}
            
            opengl_triangulated_surface->Display(in_color);
        }
    }
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<double,double,double>::
Display(const int in_color) const
{
    if(valid && draw){
        if(!draw_surface) opengl_octree_cell_scalar_field.Display(in_color);
        else{
            if(!opengl_triangulated_surface){
                OCTREE_GRID<double>& grid=opengl_octree_cell_scalar_field.grid;ARRAY<double>& phi=opengl_octree_cell_scalar_field.value;
                LEVELSET_OCTREE<double> levelset(grid,phi);
                DUALCONTOUR_OCTREE<double> contour(&levelset);
                if(triangulated_surface){delete &triangulated_surface->mesh;delete &triangulated_surface->particles;}
                delete triangulated_surface;delete opengl_triangulated_surface;
                triangulated_surface=contour.Get_Triangulated_Surface();
                opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<double>(*triangulated_surface,true,
                    OPENGL_MATERIAL::Plastic(OPENGL_COLOR((double)0.8,(double)0,0)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((double)0,0,(double)0.8)));
                *opengl_triangulated_surface->vertex_normals=*triangulated_surface->vertex_normals;}
            
            opengl_triangulated_surface->Display(in_color);
        }
    }
}
#endif
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Display(const int in_color) const
{
    if(valid && draw && !draw_surface) opengl_octree_cell_scalar_field.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_octree_cell_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Reinitialize(const bool force)
{
    if(draw||force){
        if(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            valid=false;
            delete opengl_triangulated_surface;opengl_triangulated_surface=0;
            if(triangulated_surface){delete &triangulated_surface->mesh;delete &triangulated_surface->particles;}delete triangulated_surface;triangulated_surface=0;
            std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename,frame);
            if(!FILE_UTILITIES::File_Exists(tmp_filename)) return;
            FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_octree_cell_scalar_field.value);
            opengl_octree_cell_scalar_field.Update();
            frame_loaded=frame;valid=true;}}
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_octree_cell_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Increase_Point_Size
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Increase_Point_Size()
{
    opengl_octree_cell_scalar_field.point_size++;
}
//#####################################################################
// Function Decrease_Point_Size
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Decrease_Point_Size()
{
    opengl_octree_cell_scalar_field.point_size--;
}
//#####################################################################
// Function Draw_Surface
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,T2,RW>::
Draw_Surface()
{
    draw_surface=!draw_surface;
}
//#####################################################################
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<float,int,float>;
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<float,bool,float>;
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<double,int,double>;
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<double,bool,double>;
template class OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<double,double,double>;
#endif
}
#endif
