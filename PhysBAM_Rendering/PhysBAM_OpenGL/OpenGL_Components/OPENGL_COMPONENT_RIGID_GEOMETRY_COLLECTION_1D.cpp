//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D(const std::string& basedir_input,const bool initialize_geometry)
    :OPENGL_COMPONENT("Rigid Geometry Collection 1D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_node_velocity_vectors(false),draw_point_simplices(true),
    current_selection(0),need_destroy_rigid_geometry_collection(false)
{
    if(initialize_geometry){
        need_destroy_rigid_geometry_collection=true;
        rigid_geometry_collection=new RIGID_GEOMETRY_COLLECTION<TV>(*new RIGID_GEOMETRY_PARTICLES<TV>(),0);}
    is_animation=true;
    has_init_destroy_information=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
~OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D()
{
    if(need_destroy_rigid_geometry_collection){
        delete &(rigid_geometry_collection->particles);
        delete rigid_geometry_collection;}
}
//#####################################################################
// Function Read_Hints
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Read_Hints(const std::string& filename)
{
    ARRAY<OPENGL_RIGID_BODY_HINTS,int> opengl_hints;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    Read_Binary<RW>(*input,opengl_hints);delete input;

    for(int i(1);i<=opengl_point_simplices.Size();i++) if(opengl_point_simplices(i) && i<=opengl_hints.Size()){
        opengl_point_simplices(i)->color=opengl_hints(i).material.diffuse;
        use_object_bounding_box(i)=opengl_hints(i).include_bounding_box;}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;

        if(read_geometry){
            if(!STREAM_TYPE(RW()).use_doubles)
                Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(STREAM_TYPE(RW()),basedir,frame,*rigid_geometry_collection);
            else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
                Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(STREAM_TYPE(RW()),basedir,frame,*rigid_geometry_collection);
#else
                PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
        }
        if(has_init_destroy_information) for(int i=1;i<=needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        int max_number_of_bodies(max(opengl_point_simplices.Size(),rigid_geometry_collection->particles.array_collection->Size()));
        // only enlarge array as we read in more geometry to memory
        opengl_point_simplices.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(draw_object,true);
        use_object_bounding_box.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(use_object_bounding_box,true);

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=1;i<=needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_geometry_collection->Is_Active(id));
            Create_Geometry(id);}
        else for(int i=1;i<=max_number_of_bodies;i++){if(rigid_geometry_collection->Is_Active(i)) Create_Geometry(i);} // TODO: can we figure out what bodies need_init

        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
            if(rigid_geometry_collection->Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_point_simplices.Size();id++){
            Destroy_Geometry(id);}

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Create_Geometry(const int id)
{
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection->Rigid_Geometry(id);

    if(!opengl_axes(id)) opengl_axes(id)=new OPENGL_AXES<T>();
    if(rigid_geometry.simplicial_object && !opengl_point_simplices(id)){
        opengl_point_simplices(id)=new OPENGL_POINT_SIMPLICES_1D<T>(*rigid_geometry.simplicial_object);
        opengl_point_simplices(id)->draw_vertices=true;
        opengl_point_simplices(id)->Enslave_Transform_To(*opengl_axes(id));}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<float,3> >(Convert_1d_To_3d(rigid_geometry_collection->Rigid_Geometry(id).Frame()));
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Destroy_Geometry(const int id)
{
    if(opengl_point_simplices(id)){delete opengl_point_simplices(id);opengl_point_simplices(id)=0;}
    if(opengl_axes(id)){delete opengl_axes(id);opengl_axes(id)=0;}
    draw_object(id)=false;
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Update_Object_Labels()
{
    for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
        if(draw_object(i)){
            if(opengl_point_simplices(i)){
                if(output_positions){
                    //rigid_geometry_collection->Rigid_Geometry(i).Update_Angular_Velocity();
                    opengl_point_simplices(i)->Set_Name(STRING_UTILITIES::string_sprintf("%s <%.3f>",rigid_geometry_collection->Rigid_Geometry(i).name.c_str(),rigid_geometry_collection->particles.X(i).x));}
                else opengl_point_simplices(i)->Set_Name(rigid_geometry_collection->Rigid_Geometry(i).name);}}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_bodies",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Display(const int in_color) const
{
    if(draw){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);

        if(draw_point_simplices){
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
                if(draw_object(i) && opengl_point_simplices(i)) opengl_point_simplices(i)->Display(in_color);}}

        if(show_object_names){
            glColor3f(1,1,1);
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
                if(draw_object(i) && rigid_geometry_collection->Rigid_Geometry(i).name.length()){
                    OpenGL_String(rigid_geometry_collection->particles.X(i),STRING_UTILITIES::string_sprintf("%s %f",rigid_geometry_collection->Rigid_Geometry(i).name.c_str(),rigid_geometry_collection->particles.twist(i).linear.x));}}}
        glPopAttrib();}
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i(1);i<=opengl_point_simplices.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return (draw && num_drawable_objects>0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box;
    if(draw){
        bool first=true;
        for(int i(1);i<=opengl_point_simplices.Size();i++)
            if(draw_object(i) && use_object_bounding_box(i) && opengl_point_simplices(i)){
                if(first){
                    box=opengl_point_simplices(i)->Bounding_Box();
                    first=false;}
                else{box=RANGE<VECTOR<float,3> >::Combine(box,opengl_point_simplices(i)->Bounding_Box());}}}
    return box;
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Draw_Object(int i,bool draw_it)
{
    if(i>draw_object.Size()) draw_object.Resize(i);
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Set_Object_Color
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Object_Color(int i,const OPENGL_COLOR &color_input)
{
    if(!opengl_point_simplices(i)) return;
    opengl_point_simplices(i)->color=color_input;
    opengl_point_simplices(i)->color_gray=color_input.Grayscale();
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Use_Object_Bounding_Box
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T,RW>::
Toggle_Draw_Mode()
{
    draw_point_simplices=!draw_point_simplices;
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<double,double>;
#endif
