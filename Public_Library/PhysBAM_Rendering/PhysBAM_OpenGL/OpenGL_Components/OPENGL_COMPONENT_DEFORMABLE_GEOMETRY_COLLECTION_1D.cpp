//#####################################################################
// Copyright 2009, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D(const std::string& prefix,const int start_frame,const bool initialize_geometry)
    :OPENGL_COMPONENT("Deformable Object List"),prefix(prefix),frame_loaded(-1),valid(false),use_active_list(false),display_mode(0),
    incremented_active_object(0),smooth_shading(false),selected_vertex(0),
    own_deformable_geometry(0),deformable_geometry_collection(0),real_selection(0),
    color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map())
{
    if(initialize_geometry){
        own_deformable_geometry=true;
        deformable_geometry_collection=new DEFORMABLE_GEOMETRY_COLLECTION<TV>(*(new GEOMETRY_PARTICLES<TV>()));}

    // check for per frame particles
    if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_structures",prefix.c_str(),start_frame)))
        invalidate_deformable_objects_selection_each_frame=true;
    else invalidate_deformable_objects_selection_each_frame=false;

    is_animation=true;
    Initialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D()
{
    delete color_map;
    if(own_deformable_geometry){
        delete &(deformable_geometry_collection->particles);
        delete deformable_geometry_collection;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Initialize()
{
    draw_velocity_vectors=false;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Reinitialize(bool force,bool read_geometry)
{
    if(!(draw && (force || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)))) return;
    static bool first_time=true;
    std::string frame_string=STRING_UTILITIES::string_sprintf("%s/%d/",prefix.c_str(),frame);
    std::string static_frame_string=frame_string;
    int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:-1;
    bool read_static_variables=static_frame!=-1 || first_time || !(deformable_geometry_collection && deformable_geometry_collection->structures.m);
    if(read_geometry) deformable_geometry_collection->Read(STREAM_TYPE(RW()),prefix,prefix,frame,static_frame,read_static_variables);
        
    if(read_static_variables){
        int m=deformable_geometry_collection->structures.m;active_list.Resize(m,true,true,true);
        point_simplices_1d_objects.Delete_Pointers_And_Clean_Memory();point_simplices_1d_objects.Resize(m);
        for(int i=1;i<=m;i++){
            STRUCTURE<TV>* structure=deformable_geometry_collection->structures(i);
            if(POINT_SIMPLICES_1D<T>* point_simplices_1d=dynamic_cast<POINT_SIMPLICES_1D<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": point simplices 1d\n";
                point_simplices_1d_objects(i)=new OPENGL_POINT_SIMPLICES_1D<T>(*point_simplices_1d,OPENGL_COLOR((T).5,(T).25,0));
                point_simplices_1d_objects(i)->draw_vertices=true;}
            else{if(first_time) LOG::cout<<"object "<<i<<": object unrecognized at geometry level\n";}}}

    Update_Velocity_Field();
    frame_loaded=frame;
    valid=true;
    first_time=false;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_particles",prefix.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    if(draw_input) ARRAYS_COMPUTATIONS::Fill(active_list,true);
    Reinitialize();
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Display(const int in_color) const
{
    if(!draw || !valid) return;

    for(int i=1;i<=point_simplices_1d_objects.m;i++){
        if(!active_list(i)) continue;
        glPushName(i);
        if(point_simplices_1d_objects(i)){
            glPushName(1);
            point_simplices_1d_objects(i)->Display(in_color);
            glPopName();}
        glPopName();}

    if(selected_vertex) OPENGL_SELECTION::Draw_Highlighted_Vertex(deformable_geometry_collection->particles.X(selected_vertex),selected_vertex);

    //if(draw_velocity_vectors) velocity_field.Display(in_color);
}
//#####################################################################
// Function Cycle_Display_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Cycle_Display_Mode()
{
    display_mode=(display_mode+1)%4;
}
//#####################################################################
// Function Show_Only_First
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Show_Only_First()
{
    ARRAYS_COMPUTATIONS::Fill(active_list,false);
    active_list(1)=true;
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Use_Bounding_Box() const
{
    return draw && valid && deformable_geometry_collection->structures.m>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    if(draw && valid && deformable_geometry_collection->structures.m>0){
        for(int i=1;i<=point_simplices_1d_objects.m;i++) if(point_simplices_1d_objects(i))box.Enlarge_To_Include_Box(point_simplices_1d_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Highlight_Particle
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Highlight_Particle()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Enter Particle Number: ",Highlight_Particle_Response_CB());
}
//#####################################################################
// Function Toggle_Active_Value
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Active_Value()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Enter Component Number: ",Toggle_Active_Value_Response_CB());
}
//#####################################################################
// Function Toggle_Use_Active_List
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Use_Active_List()
{
    if(active_list.m<1) return;
    use_active_list=!use_active_list;if(!use_active_list) ARRAYS_COMPUTATIONS::Fill(active_list,true);
    else{ARRAYS_COMPUTATIONS::Fill(active_list,false);incremented_active_object=1;active_list(incremented_active_object)=true;}
}
//#####################################################################
// Function Select_Segment
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Selection_Mode()
{
    if (!real_selection) return;
    assert(real_selection->body_selection);
    if (real_selection->body_selection->type == OPENGL_SELECTION::POINT_SIMPLICES_1D){
        delete real_selection->body_selection;
        real_selection->body_selection=real_selection->saved_selection;}
    else return;
    Highlight_Selection(real_selection);
    glutPostRedisplay();
}
//#####################################################################
// Function Increment_Active_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Increment_Active_Object()
{
    if(active_list.m<1 || !use_active_list) return;
    active_list(incremented_active_object)=false;incremented_active_object++;
    if(incremented_active_object>active_list.m) incremented_active_object=1;
    active_list(incremented_active_object)=true;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>* selection=0;
    if(buffer_size>=1){
        selection=new OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>(this);
        selection->body_index=buffer[0];selection->subobject_type=buffer[1];
        switch(selection->subobject_type){
            case 1:selection->subobject=point_simplices_1d_objects(buffer[0]);break;
            default:PHYSBAM_FATAL_ERROR();}
        selection->body_selection=selection->subobject->Get_Selection(&buffer[2],buffer_size-2);
        if(!selection->body_selection){delete selection;selection=0;}}
    return selection;
}
//#####################################################################
// Function Set_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type!=OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_1D) return;
    real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>*)selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type!=OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_1D) return;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>*)selection;
    if(selection->hide) active_list(real_selection->body_index)=false;
    else real_selection->subobject->Highlight_Selection(real_selection->body_selection);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Clear_Highlight()
{
    for(int i=1;i<=point_simplices_1d_objects.m;i++) if(point_simplices_1d_objects(i) && active_list(i)) point_simplices_1d_objects(i)->Clear_Highlight();
    real_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    if(selection && selection->type==OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_1D && selection->object==this){
        OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>*)selection;
        output_stream<<"Deformable object "<<real_selection->body_index<<std::endl;
        real_selection->subobject->Print_Selection_Info(output_stream,real_selection->body_selection);}
}
//#####################################################################
// Function Toggle_Active_Value_Response
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Active_Value_Response()
{
    bool start_val=true;
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        int index;
        std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);
        sstream>>index;
        if (index!=0) {
            if(active_list.m<index) active_list.Resize(index+1,true,true,start_val);
            active_list(index)=!active_list(index);}}
}
//#####################################################################
// Function Highlight_Particle_Response
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Highlight_Particle_Response()
{
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        int index=0;std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);sstream>>index;
        if(index>0 && index<=deformable_geometry_collection->particles.array_collection->Size()) selected_vertex=index;}
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Vector_Size(double size)
{
    //velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Increase_Vector_Size()
{
    //velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Decrease_Vector_Size()
{
    //velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Update_Velocity_Field()
{
    if(draw_velocity_vectors){
        positions=deformable_geometry_collection->particles.X;
        velocity_vectors=deformable_geometry_collection->particles.V;}
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T,RW>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this && invalidate_deformable_objects_selection_each_frame) delete_selection=true;
    return 0;
}
//#####################################################################
template class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<double,double>;
#endif
