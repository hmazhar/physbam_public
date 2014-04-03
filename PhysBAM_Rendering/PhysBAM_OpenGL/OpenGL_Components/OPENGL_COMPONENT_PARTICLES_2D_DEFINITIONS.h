//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;

template<class T,class T_PARTICLES,class RW> OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
OPENGL_COMPONENT_PARTICLES_2D(const std::string &filename_input, const std::string &filename_set_input, bool use_ids_input, bool particles_stored_per_cell_uniform_input)
    : OPENGL_COMPONENT("Particles 2D"), particles(new T_PARTICLES),opengl_points(new OPENGL_POINTS_2D<T>(*(new ARRAY<VECTOR<T,2> >))),
      opengl_vector_field(*(new ARRAY<VECTOR<T,2> >), opengl_points->points, OPENGL_COLOR::Cyan()), 
      filename(filename_input), filename_set(filename_set_input), frame_loaded(-1), set(1), set_loaded(-1),number_of_sets(0),use_sets(false),valid(false),
      draw_velocities(false), have_velocities(false), use_ids(use_ids_input), 
      particles_stored_per_cell_uniform(particles_stored_per_cell_uniform_input),
      draw_multiple_particle_sets(false)
{
    number_of_sets=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),(number_of_sets+1),frame);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(FILE_UTILITIES::File_Exists(filename)) number_of_sets++;
        else break;}
    if(number_of_sets>0){use_sets=true;draw_multiple_particle_sets=true;}
    else number_of_sets=1;

    particles_multiple.Resize(number_of_sets);opengl_points_multiple.Resize(number_of_sets);selected_ids.Resize(number_of_sets);
    particles_multiple(1)=particles;opengl_points_multiple(1)=opengl_points;
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Particle_Multiple_Color_Map();
    opengl_points->color=color_map->Lookup(1);
    for(int i=2;i<=number_of_sets;i++){
        particles_multiple(i)=new T_PARTICLES;
        opengl_points_multiple(i)=new OPENGL_POINTS_2D<T>(*(new ARRAY<VECTOR<T,2> >),color_map->Lookup(i));}
        
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}

template<class T,class T_PARTICLES,class RW> OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
~OPENGL_COMPONENT_PARTICLES_2D()
{
    for(int i=1;i<=number_of_sets;i++){delete particles_multiple(i);delete &opengl_points_multiple(i)->points;delete opengl_points_multiple(i);}
    delete &opengl_vector_field.vector_field;
}

template<class T,class T_PARTICLES,class RW> bool OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Display(const int in_color) const
{
    if (valid && draw){
        GLint mode;
        glGetIntegerv(GL_RENDER_MODE, &mode);
        if(mode == GL_SELECT){
            if(draw_multiple_particle_sets){
                glPushName(0);
                for(int i=1;i<=number_of_sets;i++){
                    glLoadName(i);
                    opengl_points_multiple(i)->Display(in_color);}
                glPopName();}
            else opengl_points->Display(in_color);
            if (draw_velocities && have_velocities) opengl_vector_field.Display(in_color);}
        else{
            if(draw_multiple_particle_sets) for(int i=1;i<=number_of_sets;i++) opengl_points_multiple(i)->Display(in_color);
            else opengl_points->Display(in_color);
            if (draw_velocities && have_velocities) opengl_vector_field.Display(in_color);}}
}

template<class T,class T_PARTICLES,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_points->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T_PARTICLES,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    int particle_set;int point_index;
    if(buffer_size==1){particle_set=set;point_index=buffer[0];}
    else if(buffer_size==2){particle_set=buffer[0];point_index=buffer[1];}
    else return 0;

    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection = new OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>(this);

    // We have the OPENGL_POINTS_2D index but need to find the particle index
    int particle_index=0;
    int active_count=0;
    for(particle_index=1;particle_index<=particles_multiple(particle_set)->array_collection->Size();particle_index++){
        active_count++;
        if(active_count==point_index) break;}
    selection->index=particle_index;
    selection->particle_set=particle_set;
    selection->location=particles_multiple(particle_set)->X(particle_index);
    ARRAY_VIEW<int>* ids=0;if(use_ids) ids=Get_Particles_Id_Array(particle_set);
    if(ids){
        selection->has_id=true;
        selection->id=(*ids)(particle_index);}
    else selection->has_id=false;
    return selection;
}

template<class T,class T_PARTICLES,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Get_Selection_By_Id(int id,int particle_set)
{
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection = 0;
    ARRAY_VIEW<int>* ids=0;if(use_ids) ids=Get_Particles_Id_Array(particle_set);
    if(ids){
        int index;
        if (ids->Find(id,index)){
            selection=new OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>(this);
            selection->index=index;
            selection->has_id=true;
            selection->id=id;
            selection->particle_set=particle_set;
            selection->location=particles_multiple(particle_set)->X(index);}}
    return selection;
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if (selection->type != OPENGL_SELECTION::COMPONENT_PARTICLES_2D) return;
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
    int particle_index = real_selection->index;
    ARRAY_VIEW<int>* ids=0;
    if (use_ids && (ids=Get_Particles_Id_Array(real_selection->particle_set))){
        Select_Particle_By_Id((*ids)(particle_index),real_selection->particle_set);}
    else{
        int point_index = 0;
        for (int i = 1; i <= particle_index; i++)
            point_index++;
        opengl_points_multiple(real_selection->particle_set)->Select_Point(point_index);}
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Clear_Highlight()
{
    if(use_ids) Clear_Id_Selection();
    else for(int i=1;i<=opengl_points_multiple.m;i++) opengl_points_multiple(i)->Clear_Selection();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    if(selection && selection->type == OPENGL_SELECTION::COMPONENT_PARTICLES_2D && selection->object == this){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
        
        if(!draw_multiple_particle_sets && real_selection->particle_set!=set) return;
    
        output_stream<<"Selected particle in ["<<component_name<<"("<<real_selection->particle_set<<")] (total number = "<<particles->array_collection->Size()<<")"<<std::endl;
    
        int current_index = -1;
        if(real_selection->has_id){
            output_stream<<"  Selected by id "<<real_selection->id<<std::endl;
            ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(real_selection->particle_set);
            if(!ids || !ids->Find(real_selection->id,current_index))
                output_stream<<"  Doesn't exist"<<std::endl;}
        else{
            if(real_selection->index <= particles_multiple(real_selection->particle_set)->array_collection->Size()){
                output_stream<<"  Selected by index "<<real_selection->index<<std::endl;
                current_index = real_selection->index;}}

        if(current_index > 0){
            // real_selection->index is index into particles array at time of selection.  Not very useful. current_index is more useful
            output_stream<<"current index = "<<current_index<<std::endl;
            Read_Write<T_PARTICLES,RW>::Print(output_stream,*particles_multiple(real_selection->particle_set),current_index);}}
}

template<class T,class T_PARTICLES,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    // TODO: reimplement transfering particles between objects.
    if(old_selection && old_selection->object == this){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)old_selection;
        if(real_selection->has_id){
            OPENGL_SELECTION* new_selection=Get_Selection_By_Id(real_selection->id,real_selection->particle_set);
            return new_selection;}
        else delete_selection=true;
    }
    return 0;
}

template<class T,class T_PARTICLES,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Selection_Bounding_Box(OPENGL_SELECTION *selection) const
{
    int current_index = Get_Current_Index_Of_Selection(selection);
    if(current_index != -1) return World_Space_Box(RANGE<VECTOR<float,2> >(VECTOR<float,2>(particles->X(current_index))));
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class T_PARTICLES,class RW> int OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Get_Current_Index_Of_Selection(OPENGL_SELECTION *selection) const
{
    PHYSBAM_ASSERT(selection->type == OPENGL_SELECTION::COMPONENT_PARTICLES_2D && selection->object == this);
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;

    int current_index = -1;

    if (real_selection->has_id){
        ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(real_selection->particle_set);
        if(ids) ids->Find(real_selection->id,current_index);}
    else if (real_selection->index <= particles_multiple(real_selection->particle_set)->array_collection->Size()) current_index = real_selection->index;

    return current_index;
}

template<class T,class T_PARTICLES,class RW> T_PARTICLES* OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Get_Particle_Set_Of_Selection(OPENGL_SELECTION *selection) const
{
    PHYSBAM_ASSERT(selection->type == OPENGL_SELECTION::COMPONENT_PARTICLES_2D && selection->object == this);
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection = (OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
    PHYSBAM_ASSERT(real_selection->particle_set<=particles_multiple.m);
    return particles_multiple(real_selection->particle_set);
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))) return;
    valid=true;have_velocities=false;

    for(int i=1;i<=number_of_sets;i++){
        std::string frame_filename;
        if(use_sets) frame_filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),i,frame);
        else frame_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
        
        try{
            std::istream* input_file=FILE_UTILITIES::Safe_Open_Input(frame_filename);
            if(particles_stored_per_cell_uniform){
                ARRAY<T_PARTICLES*,VECTOR<int,2> > particles_per_cell;
                Read_Binary<RW>(*input_file,particles_per_cell);
                ARRAY<ARRAY_COLLECTION*> initialization_array(particles_per_cell.array.Size());
                for(int j=1;j<=particles_per_cell.array.Size();j++){
                    if(particles_per_cell.array(j)) initialization_array(j)=particles_per_cell.array(j)->array_collection;
                    else initialization_array(j)=0;}
                ARRAY_VIEW<const ARRAY_COLLECTION* const> initialization_array_view(initialization_array.Size(),initialization_array.Get_Array_Pointer());
                particles_multiple(i)->array_collection->Initialize(initialization_array_view);
                particles_per_cell.Delete_Pointers_And_Clean_Memory();}
            else{
                Read_Binary<RW>(*input_file,*particles_multiple(i));}
            delete input_file;
            opengl_points_multiple(i)->Set_Points_From_Particles(*particles_multiple(i),true,Get_Particles_Id_Array(i)!=0);}
        catch(FILESYSTEM_ERROR&){valid=false;}}
    frame_loaded=frame;
#if 0
    if (draw_velocities && particles->store_velocity){
        have_velocities = true;
        opengl_vector_field.vector_field.Resize(particles->Size());
        for(int i=1;i<=particles->Size();i++)opengl_vector_field.vector_field(i)=particles->V(i);}
#endif
    if(use_ids) Apply_Id_Selection();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Toggle_Draw_Point_Numbers()
{
    for(int i=1;i<=opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_point_numbers=!opengl_points_multiple(i)->draw_point_numbers;
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Toggle_Draw_Radii()
{
    for(int i=1;i<=opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_radii=!opengl_points_multiple(i)->draw_radii;
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    if (draw_velocities) Reinitialize(true);
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Command_Prompt_Response()
{
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        std::string command;
        std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);
        sstream >> command;
        if (command == "s")
        {
            ARRAY<int> indices;
            int index;
            while (sstream >> index) indices.Append(index);
            if (use_ids && Get_Particles_Id_Array(set))
            {
                Clear_Id_Selection();
                Select_Particles_By_Ids(indices);
            }
            else 
            {
                opengl_points->Clear_Selection();
                opengl_points->Select_Points(indices);
            }
        }
    }
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Command_Prompt()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Command: ", Command_Prompt_Response_CB(),"");
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Next_Set()
{
    set=min(set+1,number_of_sets);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}
template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Previous_Set()
{
    set=max(set-1,1);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Toggle_Draw_Multiple_Particle_Sets()
{
    draw_multiple_particle_sets=!draw_multiple_particle_sets;
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Select_Particle_By_Id(int id, int particle_set)
{
    selected_ids(particle_set).Append(id);
    Apply_Id_Selection();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Select_Particles_By_Ids(const ARRAY<int> &ids)
{
    // TODO: implement for multiple particle sets
    for (int i = 1; i <= ids.m; i++) selected_ids(1).Append(ids(i));
    Apply_Id_Selection();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Clear_Id_Selection()
{
    for(int i=1;i<=selected_ids.m;i++) selected_ids(i).Remove_All();
    Apply_Id_Selection();
}

template<class T,class T_PARTICLES,class RW> void OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Apply_Id_Selection()
{
    for(int current_set=1;current_set<=number_of_sets;current_set++){
        opengl_points_multiple(current_set)->Clear_Selection();
        if (use_ids && valid && opengl_points_multiple(current_set)->points.m == particles_multiple(current_set)->array_collection->Size())
        {
            ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(current_set);
            if(!ids) continue;
            int idx = 1;
            for (int i = 1; i <= particles_multiple(current_set)->array_collection->Size(); i++){
                int dummy;
                if(selected_ids(current_set).Find((*ids)(i),dummy)) opengl_points_multiple(current_set)->Select_Point(idx); 
                idx++;}
        }
    }
}

template<class T,class T_PARTICLES,class RW> ARRAY_VIEW<int>* OPENGL_COMPONENT_PARTICLES_2D<T,T_PARTICLES,RW>::
Get_Particles_Id_Array(int set_number) const
{
    if(!set_number) set_number=set;
    ARRAY_VIEW<int>* ids=particles_multiple(set_number)->array_collection->template Get_Array<int>(ATTRIBUTE_ID_ID);
    if(ids && ids->Size() && (*ids)(1)) return ids; // A hack to ignore ids if the first one equals zero
    return 0;
}

//#####################################################################
// Selection object functions
//#####################################################################
template<class T>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>::
Bounding_Box() const
{
    return object->Selection_Bounding_Box((OPENGL_SELECTION*)this);
}
//#####################################################################
