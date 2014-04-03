//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>

using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
OPENGL_COMPONENT_LEVELSET_3D(const std::string& levelset_filename_input,
                             const std::string& triangulated_surface_filename_input,
                             const std::string& filename_set_input,
                             const std::string& filename_triangulated_surface_set_input,
                             bool write_generated_triangulated_surface_input,
                             bool check_triangulated_surface_file_time_input)
    : OPENGL_COMPONENT("Levelset 3D"), opengl_levelset_multiview(0), 
      levelset_filename(levelset_filename_input),triangulated_surface_filename(triangulated_surface_filename_input),
      filename_set(filename_set_input),filename_triangulated_surface_set(filename_triangulated_surface_set_input),
      write_generated_triangulated_surface(write_generated_triangulated_surface_input),
      frame_loaded(-1),check_triangulated_surface_file_time(check_triangulated_surface_file_time_input),
      set(1),set_loaded(-1),use_sets(true),draw_multiple_levelsets(true),ghost_cells(3)
{
    int number_of_sets=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,(number_of_sets+1));
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(FILE_UTILITIES::File_Exists(filename)) number_of_sets++;
        else break;}
    LOG::cout<<"Found "<<number_of_sets<<" levelsets for multiphase"<<std::endl;

    if(number_of_sets==0){use_sets=false;draw_multiple_levelsets=false;number_of_sets=1;}

    opengl_levelset_multiviews.Resize(number_of_sets);
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Levelset_Multiple_Color_Map();
    for(int i=1;i<=opengl_levelset_multiviews.m;i++){
        opengl_levelset_multiviews(i)=new OPENGL_LEVELSET_MULTIVIEW<T>();
        if(use_sets){
            OPENGL_COLOR color=color_map->Lookup(i);
            opengl_levelset_multiviews(i)->Set_Slice_Color(color,OPENGL_COLOR::Transparent());
            opengl_levelset_multiviews(i)->Set_Surface_Material(color,color);
            opengl_levelset_multiviews(i)->Set_Two_Sided(false);}}
    opengl_levelset_multiview=opengl_levelset_multiviews(1);

    if (triangulated_surface_filename.length()==0) triangulated_surface_filename="";

    is_animation=levelset_filename.find("%d")!=std::string::npos;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,
                     const OPENGL_MATERIAL &back_surface_mat)
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Surface_Material(front_surface_mat, back_surface_mat);
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat)
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Overlayed_Surface_Material(overlayed_surface_mat);
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Set_Slice_Color(const OPENGL_COLOR &inside_slice_color, const OPENGL_COLOR &outside_slice_color)
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Set_Slice_Color(inside_slice_color, outside_slice_color);
}

template<class T,class RW> bool OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    if(use_sets) return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,set));
    else return FILE_UTILITIES::File_Exists(is_animation?STRING_UTILITIES::string_sprintf(levelset_filename.c_str(),frame_input):levelset_filename);
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Display(const int in_color) const
{
    if(draw){
        if(draw_multiple_levelsets) for(int i=1;i<=opengl_levelset_multiviews.m;i++) opengl_levelset_multiviews(i)->Display(in_color);
        else opengl_levelset_multiview->Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Bounding_Box() const
{
    if (draw) return opengl_levelset_multiview->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        bool is_MAC=true;
        if(opengl_levelset_multiviews.Size() && opengl_levelset_multiviews(1)->Levelset() && !opengl_levelset_multiviews(1)->Levelset()->grid.Is_MAC_Grid()) is_MAC=false;
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_3D && is_MAC){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
            opengl_levelset_multiviews(1)->Levelset()->grid.Clamp(index,ghost_cells);
            for(int i=1;i<=opengl_levelset_multiviews.m;i++){
                const LEVELSET_3D<GRID<TV> >& levelset=*opengl_levelset_multiviews(i)->Levelset();
                output_stream<<component_name<<": phi["<<i<<"]="<<levelset.phi(index)
                             <<" curvature["<<i<<"]="<<levelset.Compute_Curvature(levelset.grid.Center(index))<<std::endl;}}
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_3D && !is_MAC){
            VECTOR<int,3> index=((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
            opengl_levelset_multiviews(1)->Levelset()->grid.Clamp(index,ghost_cells);
            for(int i=1;i<=opengl_levelset_multiviews.m;i++)  if(opengl_levelset_multiviews(i)->Levelset())
                output_stream<<component_name<<": phi["<<i<<"]="<<(*opengl_levelset_multiviews(i)->Levelset()).phi(index)<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_3D){
            OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>*)current_selection;
            VECTOR<T,3> location=selection->location;
            for(int i=1;i<=opengl_levelset_multiviews.m;i++) if(opengl_levelset_multiviews(i)->Levelset())
                output_stream<<component_name<<": phi["<<i<<"] @ particle="<<opengl_levelset_multiviews(i)->Levelset()->Phi(location)<<std::endl;}}
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Turn_Smooth_Shading_On()
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Turn_Smooth_Shading_On();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Turn_Smooth_Shading_Off()
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Turn_Smooth_Shading_Off();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && (frame_loaded != frame || set_loaded != set)) || (!is_animation && frame_loaded < 0)){
            if(use_sets){
                for(int i=1;i<=opengl_levelset_multiviews.m;i++)
                    Reinitialize_Levelset(STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,i),STRING_UTILITIES::string_sprintf(filename_triangulated_surface_set.c_str(),frame,i),opengl_levelset_multiviews(i));
                set_loaded=set;}
            else Reinitialize_Levelset(FILE_UTILITIES::Get_Frame_Filename(levelset_filename.c_str(),frame), FILE_UTILITIES::Get_Frame_Filename(triangulated_surface_filename.c_str(),frame), opengl_levelset_multiview);
            frame_loaded=frame;
        }
    }
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Reinitialize_Levelset(const std::string& levelset_filename, const std::string& triangulated_surface_filename, OPENGL_LEVELSET_MULTIVIEW<T,RW>* levelset_multiview)
{
    if(FILE_UTILITIES::File_Exists(levelset_filename)) levelset_multiview->Read_Levelset(levelset_filename);
    else return;
    if(!triangulated_surface_filename.empty() && FILE_UTILITIES::File_Exists(triangulated_surface_filename) &&
       (!check_triangulated_surface_file_time || FILE_UTILITIES::Compare_File_Times(triangulated_surface_filename, levelset_filename)>=0)){
        levelset_multiview->Read_Triangulated_Surface(triangulated_surface_filename);}
    else levelset_multiview->Generate_Triangulated_Surface(write_generated_triangulated_surface,triangulated_surface_filename);
    levelset_multiview->Initialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Toggle_Display_Overlay()
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Display_Overlay();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Toggle_Slice_Color_Mode()
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Slice_Color_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Toggle_Smooth_Slice()
{
    for(int i=1;i<=opengl_levelset_multiviews.m;i++)
        opengl_levelset_multiviews(i)->Toggle_Smooth_Slice_Texture();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Next_Set()
{
    if(!use_sets) return;
    set=min(set+1,opengl_levelset_multiviews.m);
    opengl_levelset_multiview=opengl_levelset_multiviews(set);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Previous_Set()
{
    if(!use_sets) return;
    set=max(set-1,1);
    opengl_levelset_multiview=opengl_levelset_multiviews(set);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_3D<T,RW>::
Toggle_Draw_Multiple_Levelsets()
{
    draw_multiple_levelsets=!draw_multiple_levelsets;
}

template class OPENGL_COMPONENT_LEVELSET_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_LEVELSET_3D<double,double>;
#endif
