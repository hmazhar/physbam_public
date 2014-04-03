//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
using namespace PhysBAM;
//#####################################################################

template<class T,class RW> OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
OPENGL_COMPONENT_LEVELSET_2D(const std::string& levelset_filename_input,const std::string filename_set_input)
    :OPENGL_COMPONENT("Levelset 2D"),opengl_levelset(0),levelset_filename(levelset_filename_input),filename_set(filename_set_input),
    frame_loaded(-1),set(1),use_sets(true),set_loaded(-1),valid(false),draw_multiple_levelsets(false)
{
    int number_of_sets=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_sets+1);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(FILE_UTILITIES::File_Exists(filename)) number_of_sets++;else break;}
    LOG::cout<<"Found "<<number_of_sets<<" levelsets for multiphase"<<std::endl;
    if(number_of_sets==0){use_sets=false;number_of_sets=1;}else draw_multiple_levelsets=true;

    opengl_levelsets.Resize(number_of_sets);
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Levelset_Multiple_Color_Map();
    for(int j=1;j<=opengl_levelsets.m;j++)
        opengl_levelsets(j)=new OPENGL_LEVELSET_2D<T>(*(new LEVELSET_2D<GRID<TV> >(*(new GRID<TV>),*(new ARRAY<T,VECTOR<int,2> >))),color_map->Lookup(j),OPENGL_COLOR::Transparent());
    opengl_levelset=opengl_levelsets(1);

    is_animation=FILE_UTILITIES::Is_Animated(levelset_filename);
    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
~OPENGL_COMPONENT_LEVELSET_2D()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        delete &opengl_levelsets(j)->levelset.grid;
        delete &opengl_levelsets(j)->levelset.phi;
        delete &opengl_levelsets(j)->levelset;}
}

template<class T,class RW> bool OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    if(use_sets) return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(filename_set.c_str(),set,frame_input));
    else return FILE_UTILITIES::Frame_File_Exists(levelset_filename,frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw){
        if(draw_multiple_levelsets) for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Display(in_color);
        else opengl_levelsets(set)->Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_levelset->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D && opengl_levelsets(1)->levelset.grid.Is_MAC_Grid()){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
            for(int i=1;i<=opengl_levelsets.m;i++) 
                output_stream<<component_name<<": phi["<<i<<"]="<<opengl_levelsets(i)->levelset.phi(index)
                             <<" curvature["<<i<<"]="<<opengl_levelsets(i)->levelset.Compute_Curvature(opengl_levelsets(i)->levelset.grid.Center(index))<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_2D && !opengl_levelsets(1)->levelset.grid.Is_MAC_Grid()){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;
            for(int i=1;i<=opengl_levelsets.m;i++) 
                output_stream<<component_name<<": phi["<<i<<"]="<<opengl_levelsets(i)->levelset.phi(index)<<std::endl;}
        if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
            OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)current_selection;
            VECTOR<T,2> location=selection->location;
            for(int i=1;i<=opengl_levelsets.m;i++) 
                output_stream<<component_name<<": phi["<<i<<"] @ particle="<<opengl_levelsets(i)->levelset.Phi(location)<<std::endl;}}
}


template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Reinitialize(const bool force_even_if_not_drawn)
{
    if (draw||force_even_if_not_drawn){
        if ((is_animation && (frame_loaded!=frame || set_loaded!=set)) || (!is_animation && frame_loaded<0)){
            valid=false;std::string filename;
            if(use_sets) for(int i=1;i<=opengl_levelsets.m;i++){
                filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),i,frame);
                if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename.c_str(),opengl_levelsets(i)->levelset);
                else return;
                opengl_levelsets(i)->Update();}
            else{
                filename=FILE_UTILITIES::Get_Frame_Filename(levelset_filename,frame);
                if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File<RW>(filename.c_str(),opengl_levelset->levelset);
                else return;
                opengl_levelset->Update();}
            frame_loaded=frame;set_loaded=set;valid=true;}}
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Color_Mode()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Color_Map();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Smooth()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Smooth_Texture();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Normals()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Normals();
}

template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Draw_Mode()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        int mask=((int)opengl_levelsets(j)->draw_area<<2) | ((int)opengl_levelsets(j)->draw_curve<<1) | ((int)opengl_levelsets(j)->draw_cells);
        int newmask=(mask%8)+1;
        opengl_levelsets(j)->draw_area=(newmask&4)!=0;
        opengl_levelsets(j)->draw_curve=(newmask&2)!=0;
        opengl_levelsets(j)->draw_cells=(newmask&1)!=0;
        opengl_levelsets(j)->Update();}
}
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Draw_Sign()
{
    for(int j=1;j<=opengl_levelsets.m;j++){
        opengl_levelsets(j)->dominant_sign=(opengl_levelsets(j)->dominant_sign==1)?-1:1;
        opengl_levelsets(j)->Update();}
}
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Next_Set()
{
    set=min(set+1,opengl_levelsets.m);
    LOG::cout<<"viewing levelset set "<<set<<std::endl;
    opengl_levelset=opengl_levelsets(set);
    Reinitialize();
}
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Previous_Set()
{
    set=max(set-1,1);
    LOG::cout<<"viewing levelset set "<<set<<std::endl;
    opengl_levelset=opengl_levelsets(set);
    Reinitialize();
}
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Draw_Multiple_Levelsets()
{
    draw_multiple_levelsets=!draw_multiple_levelsets;
}
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_2D<T,RW>::
Toggle_Draw_Ghost_Values()
{
    for(int j=1;j<=opengl_levelsets.m;j++) opengl_levelsets(j)->Toggle_Draw_Ghost_Values();
}

template class OPENGL_COMPONENT_LEVELSET_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_LEVELSET_2D<double,double>;
#endif
