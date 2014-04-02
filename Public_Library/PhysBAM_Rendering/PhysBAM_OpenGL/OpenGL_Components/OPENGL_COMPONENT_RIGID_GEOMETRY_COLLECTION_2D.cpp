//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D(RIGID_GEOMETRY_PARTICLES<TV>* particles_input,const std::string& basedir_input)
    :OPENGL_COMPONENT("Rigid Geometry Collection 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    rigid_geometry_collection(new RIGID_GEOMETRY_COLLECTION<TV>(*particles_input,0,0)),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true),
    current_selection(0),need_destroy_rigid_geometry_collection(true)
{
    is_animation=true;
    has_init_destroy_information=false;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT("Rigid Geometry Collection 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    rigid_geometry_collection(&rigid_geometry_collection),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true),
    current_selection(0),need_destroy_rigid_geometry_collection(false)
{
    is_animation=true;
    has_init_destroy_information=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
~OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D()
{
    if(need_destroy_rigid_geometry_collection) delete &rigid_geometry_collection;
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Draw_Mode(const int mode)
{
    draw_segmented_curve=(mode&0x01)!=0;
    draw_triangulated_area=(mode&0x02)!=0;
    draw_implicit_curve=(mode&0x04)!=0;
}
//#####################################################################
// Function Get_Draw_Mode
//#####################################################################
template<class T,class RW> int OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Get_Draw_Mode() const
{
    return 0x04*(int)(draw_implicit_curve)+0x02*(int)(draw_triangulated_area)+0x01*(int)(draw_segmented_curve);
}
//#####################################################################
// Function Read_Hints
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Read_Hints(const std::string& filename)
{
    ARRAY<OPENGL_RIGID_BODY_HINTS,int> opengl_hints;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    Read_Binary<RW>(*input,opengl_hints);delete input;

    for(int i(1);i<=opengl_segmented_curve.Size();i++) if(opengl_segmented_curve(i) && i<=opengl_hints.Size()){
        opengl_segmented_curve(i)->color=opengl_hints(i).material.diffuse;
        use_object_bounding_box(i)=opengl_hints(i).include_bounding_box;}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;

        if(read_geometry){
            if(!STREAM_TYPE(RW()).use_doubles)
                Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(STREAM_TYPE(RW()),basedir,frame,*rigid_geometry_collection,&needs_init,&needs_destroy);
            else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
                Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(STREAM_TYPE(RW()),basedir,frame,*rigid_geometry_collection,&needs_init,&needs_destroy);
#else
            PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
        }
        if(has_init_destroy_information) for(int i=1;i<=needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        int max_number_of_bodies(max(opengl_segmented_curve.Size(),rigid_geometry_collection->particles.array_collection->Size()));
        // only enlarge array as we read in more geometry to memory
        opengl_segmented_curve.Resize(max_number_of_bodies);
        opengl_triangulated_area.Resize(max_number_of_bodies);
        opengl_levelset.Resize(max_number_of_bodies);
        extra_components.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(draw_object,true);
        use_object_bounding_box.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(use_object_bounding_box,true);

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=1;i<=needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_geometry_collection->Is_Active(id));
            Create_Geometry(id);}

        // Only display real bodies (not ghost bodies)
        if (FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/partition",basedir.c_str(),frame))) {
            ARRAY<int> particles_of_this_partition;
            FILE_UTILITIES::template Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/partition",basedir.c_str(),frame),particles_of_this_partition);
            for (int i(1);i<=max_number_of_bodies;i++)
                draw_object(i)=false;
            for (int i(1);i<=particles_of_this_partition.Size();i++)
                draw_object(particles_of_this_partition(i))=true;}

        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
            if(rigid_geometry_collection->Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_segmented_curve.Size();id++) Destroy_Geometry(id);
        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/colors",basedir.c_str(),frame)))
            FILE_UTILITIES::template Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/colors",basedir.c_str(),frame),colors);
        for(int id(1);id<=colors.m;id++){
            if(colors(id)==1) Set_Object_Color(id,OPENGL_COLOR::Green());
            if(colors(id)==2) Set_Object_Color(id,OPENGL_COLOR::Magenta());}

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Create_Geometry(const int id)
{
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection->Rigid_Geometry(id);

    if(!opengl_axes(id)){opengl_axes(id)=new OPENGL_AXES<T>();}
    if(rigid_geometry.simplicial_object && !opengl_segmented_curve(id)){
        opengl_segmented_curve(id)=new OPENGL_SEGMENTED_CURVE_2D<T>(*rigid_geometry.simplicial_object);
        opengl_segmented_curve(id)->draw_vertices=true;
        opengl_segmented_curve(id)->Enslave_Transform_To(*opengl_axes(id));}

    // add triangulated area
    TRIANGULATED_AREA<T>* triangulated_area=rigid_geometry.template Find_Structure<TRIANGULATED_AREA<T>*>();
    if(triangulated_area && !opengl_triangulated_area(id)){
        triangulated_area->mesh.Initialize_Segment_Mesh();
        opengl_triangulated_area(id)=new OPENGL_TRIANGULATED_AREA<T>(*triangulated_area);
        opengl_triangulated_area(id)->Enslave_Transform_To(*opengl_axes(id));}

    if(rigid_geometry.implicit_object && !opengl_levelset(id)){ // ASSUMES LEVELSET_IMPLICIT_CURVE!
        IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_geometry.implicit_object->object_space_implicit_object;
        if(LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            opengl_levelset(id)=new OPENGL_LEVELSET_2D<T>(levelset_implicit_object->levelset);
            // TODO: fix me
            //opengl_levelset(id)->Set_Color_Mode(OPENGL_LEVELSET_2D<T>::COLOR_GRADIENT);
            opengl_levelset(id)->draw_cells=true;opengl_levelset(id)->draw_curve=false;opengl_levelset(id)->draw_area=false;
            opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));}}

    // add extra components
    if(opengl_triangulated_area(id)){
        std::string filename_pattern=STRING_UTILITIES::string_sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            LOG::cout<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>* component=
                new OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T,RW>(*triangulated_area,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            extra_components(id).Append(component);}}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<float,3> >(Convert_2d_To_3d(rigid_geometry_collection->Rigid_Geometry(id).Frame()));
    if(opengl_triangulated_area(id)){
        std::string color_map_filename=STRING_UTILITIES::string_sprintf("%s/%d/stress_map_of_triangulated_area_%d",basedir.c_str(),frame,id);
        if(FILE_UTILITIES::File_Exists(color_map_filename)){
            if(!opengl_triangulated_area(id)->color_map) opengl_triangulated_area(id)->color_map=new ARRAY<OPENGL_COLOR>;
            FILE_UTILITIES::Read_From_File<RW>(color_map_filename,*opengl_triangulated_area(id)->color_map);}
        else if(opengl_triangulated_area(id)->color_map){delete opengl_triangulated_area(id)->color_map;opengl_triangulated_area(id)->color_map=0;}}
    for(int i=1;i<=extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Destroy_Geometry(const int id)
{
    if(opengl_segmented_curve(id)){delete opengl_segmented_curve(id);opengl_segmented_curve(id)=0;}
    if(opengl_triangulated_area(id)){delete opengl_triangulated_area(id)->color_map;delete opengl_triangulated_area(id);opengl_triangulated_area(id)=0;}
    if(opengl_levelset(id)){delete opengl_levelset(id);opengl_levelset(id)=0;}
    if(opengl_axes(id)){delete opengl_axes(id);opengl_axes(id)=0;}
    extra_components(id).Delete_Pointers_And_Clean_Memory();
    draw_object(id)=false;
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Update_Object_Labels()
{
    int number_of_drawn_bodies=draw_object.Count_Matches(true);
    if(draw_velocity_vectors){
        positions.Resize(number_of_drawn_bodies);
        velocity_vectors.Resize(number_of_drawn_bodies);}
    if(draw_node_velocity_vectors){
        node_positions.Resize(number_of_drawn_bodies*4);
        node_velocity_vectors.Resize(number_of_drawn_bodies*4);}

    int idx=0;
    for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
        if(draw_object(i)){
            if(draw_velocity_vectors || draw_node_velocity_vectors){idx++;
                if(draw_velocity_vectors){
                    positions(idx)=rigid_geometry_collection->particles.X(i);
                    velocity_vectors(idx)=rigid_geometry_collection->particles.twist(i).linear;}
                if(draw_node_velocity_vectors){
                    //only valid for squares . . .
                    RIGID_GEOMETRY<TV>* rigid_geometry=&rigid_geometry_collection->Rigid_Geometry(i);

                    node_positions((idx-1)*4+1)=rigid_geometry->World_Space_Point(VECTOR<T,2>(1,1));
                    node_velocity_vectors((idx-1)*4+1)=rigid_geometry->Pointwise_Object_Velocity(rigid_geometry->World_Space_Vector(VECTOR<T,2>(1,1)));

                    node_positions((idx-1)*4+2)=rigid_geometry->World_Space_Point(VECTOR<T,2>(1,-1));
                    node_velocity_vectors((idx-1)*4+2)=rigid_geometry->Pointwise_Object_Velocity(rigid_geometry->World_Space_Vector(VECTOR<T,2>(1,-1)));

                    node_positions((idx-1)*4+3)=rigid_geometry->World_Space_Point(VECTOR<T,2>(-1,1));
                    node_velocity_vectors((idx-1)*4+3)=rigid_geometry->Pointwise_Object_Velocity(rigid_geometry->World_Space_Vector(VECTOR<T,2>(-1,1)));

                    node_positions((idx-1)*4+4)=rigid_geometry->World_Space_Point(VECTOR<T,2>(-1,-1));
                    node_velocity_vectors((idx-1)*4+4)=rigid_geometry->Pointwise_Object_Velocity(rigid_geometry->World_Space_Vector(VECTOR<T,2>(-1,-1)));}}
            if(opengl_segmented_curve(i)){
                if(output_positions){
                    opengl_segmented_curve(i)->Set_Name(STRING_UTILITIES::string_sprintf("%s <%.3f %.3f> [w=%.3f]",rigid_geometry_collection->Rigid_Geometry(i).name.c_str(),rigid_geometry_collection->particles.X(i).x,rigid_geometry_collection->particles.X(i).y,rigid_geometry_collection->particles.twist(i).angular.x));}
                else opengl_segmented_curve(i)->Set_Name(rigid_geometry_collection->Rigid_Geometry(i).name);}}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Display(const int in_color) const
{
    if(draw){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);

        GLint mode;
        glGetIntegerv(GL_RENDER_MODE,&mode);

        if(draw_segmented_curve){
            glPushName(1);
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_segmented_curve(i)) opengl_segmented_curve(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_triangulated_area){
            glPushName(2);
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_triangulated_area(i)) opengl_triangulated_area(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_implicit_curve){
            glPushName(3);
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++){
                glPushName(Value(i));
                if(draw_object(i) && opengl_levelset(i)) opengl_levelset(i)->Display(in_color);
                glPopName();}
            glPopName();}
        if(draw_individual_axes)
            for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++)
                if(draw_object(i) && opengl_axes(i)) opengl_axes(i)->Display(in_color);

        if(mode!=GL_SELECT){
            if(draw_velocity_vectors) velocity_field.Display(in_color);
            if(draw_node_velocity_vectors) node_velocity_field.Display(in_color);

            for(int i(1);i<=extra_components.Size();i++)
                for(int j=1;j<=extra_components(i).m;j++)
                    extra_components(i)(j)->Display(in_color);

            if(show_object_names){
                glColor3f(1,1,1);
                for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++)
                    if(draw_object(i) && rigid_geometry_collection->Rigid_Geometry(i).name.length())
                        OpenGL_String(rigid_geometry_collection->particles.X(i),rigid_geometry_collection->Rigid_Geometry(i).name);}}
        glPopAttrib();}
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i(1);i<=opengl_segmented_curve.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return (draw && num_drawable_objects>0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box;
    if(draw){
        bool first=true;
        for(int i(1);i<=opengl_segmented_curve.Size();i++)
            if(draw_object(i) && use_object_bounding_box(i) && opengl_segmented_curve(i)){
                if(first){
                    box=opengl_segmented_curve(i)->Bounding_Box();
                    first=false;}
                else{box=RANGE<VECTOR<float,3> >::Combine(box,opengl_segmented_curve(i)->Bounding_Box());}}}
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size>=2){
        if(buffer[0]==1 || buffer[0]==2){
            int body_id(buffer[1]);
            OPENGL_SELECTION* body_selection=0;
            if(buffer[0]==1){ // segmented curve
                PHYSBAM_ASSERT(opengl_segmented_curve(body_id));
                body_selection=opengl_segmented_curve(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
            else{ // triangulated area
                PHYSBAM_ASSERT(opengl_triangulated_area(body_id));
                body_selection=opengl_triangulated_area(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
            if(body_selection) selection=new OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>(this,body_id,body_selection);}}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>*)selection;
        if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D ||
           real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){ // segmented curve
            if(opengl_segmented_curve(real_selection->body_id)) // might have become inactive
                opengl_segmented_curve(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){ // triangulated area
            if(opengl_triangulated_area(real_selection->body_id)) // might have become inactive
                opengl_triangulated_area(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
    for(int i(1);i<=opengl_segmented_curve.Size();i++)
        if(opengl_segmented_curve(i)) opengl_segmented_curve(i)->Clear_Highlight();
    for(int i(1);i<=opengl_triangulated_area.Size();i++)
        if(opengl_triangulated_area(i)) opengl_triangulated_area(i)->Clear_Highlight();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(!selection) return;
    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>*)selection;

        output_stream<<"Rigid body "<<real_selection->body_id<<std::endl;

        // handle case of body no longer being active
        if(!rigid_geometry_collection->Is_Active(real_selection->body_id) || !opengl_segmented_curve(real_selection->body_id)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_GEOMETRY<TV> *body=const_cast<RIGID_GEOMETRY<TV>*>(&rigid_geometry_collection->Rigid_Geometry(real_selection->body_id));
        //body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name = "<<body->name<<std::endl;}
        Read_Write<ARRAY_COLLECTION,RW>::Print(output_stream,*rigid_geometry_collection->particles.array_collection,real_selection->body_id);

        MATRIX<T,3> body_transform=body->Frame().Matrix();

        // Have to do this ourselves here because we need to transform particle positions to world space
        if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D ||
           real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){ // segmented curve
            if(opengl_segmented_curve(real_selection->body_id))
                opengl_segmented_curve(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);

            if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)real_selection->body_selection)->index;
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(index)))<<std::endl;}
            else if(real_selection->body_selection->type==OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D){
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)real_selection->body_selection)->index;
                int node1,node2;body->simplicial_object->mesh.elements(index).Get(node1,node2);
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(node1)))<<std::endl;
                output_stream<<"Pointwise object velocity "<<body->Pointwise_Object_Velocity(body->World_Space_Point(body->simplicial_object->particles.X(node2)))<<std::endl;}}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT ||
                real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE){
            if(opengl_triangulated_area(real_selection->body_id))
                opengl_triangulated_area(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);}}
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Draw_Object(int i,bool draw_it)
{
    if(i>draw_object.Size()) draw_object.Resize(i);
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Set_Object_Color
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Object_Color(int i,const OPENGL_COLOR &color_input)
{
    if(!opengl_segmented_curve(i)) return;
    opengl_segmented_curve(i)->color=color_input;
    opengl_segmented_curve(i)->vertex_color=color_input;
    opengl_segmented_curve(i)->color_gray=color_input.Grayscale();
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Use_Object_Bounding_Box
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
    node_velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Individual_Axes
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Individual_Axes()
{
    draw_individual_axes=!draw_individual_axes;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Toggle_Node_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Node_Velocity_Vectors()
{
    draw_node_velocity_vectors=!draw_node_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
    node_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
    node_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T,RW>::
Toggle_Draw_Mode()
{
    int mode=Get_Draw_Mode();
    mode=(mode+1)%8;
    Set_Draw_Mode(mode);
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<double,double>;
#endif
