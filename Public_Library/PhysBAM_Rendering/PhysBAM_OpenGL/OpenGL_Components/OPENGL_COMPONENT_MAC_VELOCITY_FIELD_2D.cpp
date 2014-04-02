//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D(const GRID<TV> &grid,const std::string &velocity_filename_input,const std::string directory_adaptive_input,const std::string filename_active_cells_input,const std::string filename_active_faces_input)
    : OPENGL_COMPONENT("MAC Velocity Field 2D"),draw_vorticity(false),
    velocity_filename(velocity_filename_input),directory_adaptive(directory_adaptive_input),filename_active_cells(filename_active_cells_input),filename_active_faces(filename_active_faces_input),level(1),
    use_levels(true),level_loaded(-1),valid(false),draw_divergence(false),draw_all_levels(true),draw_streamlines(false),use_seed_for_streamlines(false),opengl_divergence_field(0),
    streamlines(*new SEGMENT_MESH(),*new GEOMETRY_PARTICLES<TV>()),opengl_streamlines(streamlines),psi_N_psi_D_basedir(""),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    ARRAY<GRID<TV>*> grid_array;
    grid_array.Resize(1);
    grid_array(1)=new GRID<TV>(grid);
    Initialize(grid_array);
}

template<class T,class RW> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D(const GRID<TV> &grid,const std::string &velocity_filename_input,const std::string directory_adaptive_input,const std::string filename_active_cells_input,const std::string filename_active_faces_input,const ARRAY<GRID<TV>*> grid_array_input)
    : OPENGL_COMPONENT("MAC Velocity Field 2D"),draw_vorticity(false),velocity_filename(velocity_filename_input),directory_adaptive(directory_adaptive_input),filename_active_cells(filename_active_cells_input),filename_active_faces(filename_active_faces_input),level(1),use_levels(true),level_loaded(-1),
    valid(false),draw_divergence(false),draw_all_levels(true),draw_streamlines(false),use_seed_for_streamlines(false),opengl_divergence_field(0),streamlines(*new SEGMENT_MESH(),*new GEOMETRY_PARTICLES<TV>()),opengl_streamlines(streamlines),psi_N_psi_D_basedir(""),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    Initialize(grid_array_input);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Initialize(const ARRAY<GRID<TV>*> &grid_array_input)
{
    is_animation=FILE_UTILITIES::Is_Animated(velocity_filename);
    frame_loaded=-1;

    int number_of_levels=0;
    while(directory_adaptive!=""){
        std::string dirname=STRING_UTILITIES::string_sprintf(directory_adaptive.c_str(),(number_of_levels+1));
        LOG::cout<<"Checking "<<dirname<<std::endl;
        if(FILE_UTILITIES::Directory_Exists(dirname)) number_of_levels++;
        else break;}
    LOG::cout<<"Found "<<number_of_levels<<" levels for adaptive"<<std::endl;
    
    if(number_of_levels==0){use_levels=false;number_of_levels=1;}
    else velocity_filename=FILE_UTILITIES::Get_Short_Name(velocity_filename);

    opengl_adaptive_mac_velocity_fields.Resize(number_of_levels);
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)=new OPENGL_MAC_VELOCITY_FIELD_2D<T>(*(new GRID<TV>(*grid_array_input(i))),*(new ARRAY<T,FACE_INDEX<2> >));
    opengl_mac_velocity_field=opengl_adaptive_mac_velocity_fields(1);
    number_of_steps=2*opengl_mac_velocity_field->grid.counts.x;
    opengl_vorticity_magnitude=new OPENGL_SCALAR_FIELD_2D<T>(opengl_mac_velocity_field->grid,*(new ARRAY<T,VECTOR<int,2> >),OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));

    OPENGL_COLOR_RAMP<T>* ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(-1e+2,OPENGL_COLOR::Red());
    ramp->Add_Color(-1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(-1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(0,OPENGL_COLOR::Black());
    ramp->Add_Color(1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(1e+2,OPENGL_COLOR::Red());
    opengl_divergence_field=new OPENGL_SCALAR_FIELD_2D<T>(opengl_mac_velocity_field->grid,divergence,ramp);

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D()
{
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++){
        delete &opengl_adaptive_mac_velocity_fields(i)->grid;
        delete &opengl_adaptive_mac_velocity_fields(i)->u;
        delete &opengl_adaptive_mac_velocity_fields(i)->v;}
    delete opengl_divergence_field;
}

template<class T,class RW> bool OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    if(use_levels) return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(directory_adaptive.c_str(),level)+FILE_UTILITIES::Get_Frame_Filename(velocity_filename.c_str(), frame_input));
    else return FILE_UTILITIES::Frame_File_Exists(velocity_filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_mac_velocity_field->Print_Selection_Info(stream,selection);
        if(draw_vorticity && ((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)->index;
            stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude->values(index)<<std::endl;}}
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid){
        if(draw){
            if(draw_all_levels) for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)->Display(in_color);
            else opengl_adaptive_mac_velocity_fields(level)->Display(in_color);}
        if(draw_divergence) opengl_divergence_field->Display(in_color);
        if(draw_vorticity) opengl_vorticity_magnitude->Display(in_color);
        if(draw_streamlines) opengl_streamlines.Display(in_color);}
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_mac_velocity_field->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Reinitialize()
{
    if (draw || draw_divergence){
        if ((is_animation && (frame_loaded!=frame || level_loaded!=level)) || (!is_animation && frame_loaded < 0)){
            valid = false;
            if(use_levels) for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++){
                std::string tmp_filename=STRING_UTILITIES::string_sprintf(directory_adaptive.c_str(),i)+FILE_UTILITIES::Get_Frame_Filename(velocity_filename.c_str(),frame);
                if (FILE_UTILITIES::File_Exists(tmp_filename)) FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_adaptive_mac_velocity_fields(i)->u,opengl_adaptive_mac_velocity_fields(i)->v);
                else return;
                tmp_filename=STRING_UTILITIES::string_sprintf(directory_adaptive.c_str(),i)+FILE_UTILITIES::Get_Frame_Filename(filename_active_cells.c_str(),frame);
                LOG::cout<<"Reading active cells from"<<tmp_filename<<std::endl;
                if(FILE_UTILITIES::File_Exists(tmp_filename)){if(!opengl_adaptive_mac_velocity_fields(i)->active_cells) opengl_adaptive_mac_velocity_fields(i)->active_cells=new T_ARRAYS_BOOL();
                    FILE_UTILITIES::Read_From_File<bool>(tmp_filename,*opengl_adaptive_mac_velocity_fields(i)->active_cells);}
                else return;
                tmp_filename=STRING_UTILITIES::string_sprintf(directory_adaptive.c_str(),i)+FILE_UTILITIES::Get_Frame_Filename(filename_active_faces.c_str(),frame);
                LOG::cout<<"Reading active faces from"<<tmp_filename<<std::endl;
                if(FILE_UTILITIES::File_Exists(tmp_filename)){if(!opengl_adaptive_mac_velocity_fields(i)->active_faces) opengl_adaptive_mac_velocity_fields(i)->active_faces=new T_FACE_ARRAYS_BOOL();
                    FILE_UTILITIES::Read_From_File<bool>(tmp_filename,*opengl_adaptive_mac_velocity_fields(i)->active_faces);}
                opengl_adaptive_mac_velocity_fields(i)->Update();}
            else{
                std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(velocity_filename.c_str(), frame);
                if (FILE_UTILITIES::File_Exists(tmp_filename)) FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_mac_velocity_field->face_velocities);
                else return;
                opengl_mac_velocity_field->Update();}
            frame_loaded=frame;level_loaded=level;valid=true;
            Update_Divergence();
            Update_Streamlines();
            if(draw_vorticity){
                Update_Vorticity();
                opengl_vorticity_magnitude->Update();}
        }
    }
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Divergence()
{
    if(draw_divergence && valid){
        GRID<TV>& grid=opengl_mac_velocity_field->grid;
        ARRAY_VIEW<T,VECTOR<int,2> > &u=opengl_mac_velocity_field->u,&v=opengl_mac_velocity_field->v;
        static T_FACE_ARRAYS_BOOL psi_N;
        static T_ARRAYS_BOOL psi_D;
        bool got_all_psi=true;
        if(!psi_N_psi_D_basedir.empty()){
            std::string psi_N_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_N",psi_N_psi_D_basedir.c_str(),frame);
            std::string psi_D_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_D",psi_N_psi_D_basedir.c_str(),frame);
            if(FILE_UTILITIES::File_Exists(psi_N_filename)) FILE_UTILITIES::Read_From_File<RW>(psi_N_filename,psi_N);
            else got_all_psi=false;
            if(FILE_UTILITIES::File_Exists(psi_D_filename)) FILE_UTILITIES::Read_From_File<RW>(psi_D_filename,psi_D);
            else got_all_psi=false;}
        else got_all_psi=false;
        if(!got_all_psi){psi_N.Clean_Memory();psi_D.Clean_Memory();}
        divergence.Resize(1,grid.counts.x,1,grid.counts.y);
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
            if(got_all_psi && (psi_D(i,j) || (psi_N.Component(1)(i,j) && psi_N.Component(1)(i+1,j) && psi_N.Component(2)(i,j) && psi_N.Component(2)(i,j+1)))) divergence(i,j)=0;
            else divergence(i,j)=grid.one_over_dX.x*(u(i+1,j)-u(i,j))+grid.one_over_dX.y*(v(i,j+1)-v(i,j));}
        opengl_divergence_field->Update();}
    else divergence.Clean_Memory();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Streamlines()
{
    streamlines.Clean_Memory();streamlines.particles.array_collection->Clean_Memory();streamlines.mesh.Clean_Memory();
    if(!draw_streamlines || !valid) return;
    
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    int number_of_streamlines=100;
    T step_length=(T).5*grid.Minimum_Edge_Length();

    RANDOM_NUMBERS<T> random;
    if(use_seed_for_streamlines) random.Set_Seed(streamline_seed);
    T_LINEAR_INTERPOLATION_VECTOR linear_interpolation;
    T_FACE_ARRAYS_SCALAR mac_velocity_field(grid);
    mac_velocity_field.Component(1)=opengl_mac_velocity_field->u;
    mac_velocity_field.Component(2)=opengl_mac_velocity_field->v;
    FACE_LOOKUP_UNIFORM<GRID<TV> > V_lookup(mac_velocity_field);

    for(int i=1;i<=number_of_streamlines;i++){
        int p=streamlines.particles.array_collection->Add_Element();
        TV X=streamlines.particles.X(p)=random.Get_Uniform_Vector(grid.domain);
        for(int step=1;step<=number_of_steps;step++){
            TV velocity=linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X);
            TV X_new=X+step_length*velocity;
            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X_new));
            X_new=grid.Clamp(X+step_length*velocity);
            int new_particle=streamlines.particles.array_collection->Add_Element();
            streamlines.particles.X(new_particle)=X_new;
            streamlines.mesh.elements.Append(VECTOR<int,2>(p,new_particle));
            p=new_particle;
            X=X_new;}}
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode()
{
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)->Toggle_Velocity_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode_And_Draw()
{
    if (draw)
    {
        Toggle_Velocity_Mode();
        if ((int)opengl_mac_velocity_field->velocity_mode==0) Toggle_Draw();
    }
    else Toggle_Draw();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Increase_Vector_Size()
{
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)->Scale_Vector_Size(1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Decrease_Vector_Size()
{
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)->Scale_Vector_Size(1/1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Arrowhead()
{
    for(int i=1;i<=opengl_adaptive_mac_velocity_fields.m;i++) opengl_adaptive_mac_velocity_fields(i)->draw_arrowhead = !opengl_adaptive_mac_velocity_fields(i)->draw_arrowhead;
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Divergence()
{
    draw_divergence=!draw_divergence;
    Update_Divergence();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_All_Levels()
{
    draw_all_levels=!draw_all_levels;
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Streamlines()
{
    draw_streamlines=!draw_streamlines;
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Use_Streamline_Seed()
{
    use_seed_for_streamlines=!use_seed_for_streamlines;
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Streamline_Seed(const unsigned int seed)
{
    streamline_seed=seed;
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Lengthen_Streamlines()
{
    number_of_steps+=10;
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Shorten_Streamlines()
{
    number_of_steps=max(number_of_steps-10,0);
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Next_Level(){
    level=min(level+1,opengl_adaptive_mac_velocity_fields.m);
    opengl_mac_velocity_field=opengl_adaptive_mac_velocity_fields(level);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Previous_Level(){
    level=max(level-1,1);
    opengl_mac_velocity_field=opengl_adaptive_mac_velocity_fields(level);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Vorticity()
{
    draw_vorticity=!draw_vorticity;
    if(draw_vorticity) valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude->Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Vorticity()
{
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    RANGE<VECTOR<int,2> > domain_indices(grid.Domain_Indices());domain_indices.Change_Size(-VECTOR<int,2>::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup(opengl_mac_velocity_field->face_velocities);
    opengl_vorticity_magnitude->values.Resize(grid.Domain_Indices());
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,domain_indices);iterator.Valid();iterator.Next()){VECTOR<int,2> index=iterator.Cell_Index();
        T vorticity_magnitude=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index).Magnitude();
        opengl_vorticity_magnitude->values(index)=vorticity_magnitude;
        min_vorticity=min(min_vorticity,vorticity_magnitude);
        max_vorticity=max(max_vorticity,vorticity_magnitude);}
}
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<double,double>;
#endif
