//#####################################################################
// Copyright 2004-2009, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D(const std::string& prefix,const int start_frame,const bool initialize_geometry)
    :OPENGL_COMPONENT("Deformable Object List"),prefix(prefix),frame_loaded(-1),valid(false),use_active_list(false),hide_unselected(false),display_mode(0),display_relative_velocity_mode(0),number_of_segmented_curve(1),
    incremented_active_object(0),smooth_shading(false),selected_vertex(0),
    own_deformable_geometry(0),deformable_geometry(0),real_selection(0),
    has_tetrahedralized_volumes(false),has_hexahedralized_volumes(false),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),.25,false,false),
    color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map())
{
    if(initialize_geometry){
        own_deformable_geometry=true;
        deformable_geometry=new DEFORMABLE_GEOMETRY_COLLECTION<TV>(*(new GEOMETRY_PARTICLES<TV>()));}

    // check for per frame particles
    if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_structures",prefix.c_str(),start_frame)))
        invalidate_deformable_objects_selection_each_frame=true;
    else invalidate_deformable_objects_selection_each_frame=false;

    color_map_relative_velocity=new OPENGL_COLOR_RAMP<T>();
    color_map_relative_velocity->Add_Color(-.05,OPENGL_COLOR(1,0,0));
    color_map_relative_velocity->Add_Color(0,OPENGL_COLOR(0.5,0.5,0.5));
    color_map_relative_velocity->Add_Color(.05,OPENGL_COLOR(0,0,1));

    is_animation=true;
    Initialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D()
{
    delete color_map;
    if(own_deformable_geometry){
        delete &(deformable_geometry->particles);
        delete deformable_geometry;}
    delete color_map_relative_velocity;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Initialize()
{
    draw_velocity_vectors=false;
}
//#####################################################################
// Function Set_Material
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Material(const int object,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material)
{
    if(triangulated_surface_objects(object)){
        triangulated_surface_objects(object)->front_material=front_material;
        triangulated_surface_objects(object)->back_material=back_material;
        triangulated_surface_objects(object)->Set_Two_Sided(true);}
    else if(tetrahedralized_volume_objects(object))
        tetrahedralized_volume_objects(object)->material=front_material;
    else if(hexahedralized_volume_objects(object))
        hexahedralized_volume_objects(object)->material=front_material;
}
//#####################################################################
// Function Set_All_Materials
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_All_Materials(const OPENGL_MATERIAL& meshfront,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material)
{
    for(int i=1;i<=deformable_geometry->structures.m;i++) Set_Material(i,front_material,back_material);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Reinitialize(bool force,bool read_geometry)
{
    if(!(draw && (force || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)))) return;
    static bool first_time=true;
    std::string frame_string=STRING_UTILITIES::string_sprintf("%s/%d/",prefix.c_str(),frame);
    std::string static_frame_string=frame_string;
    int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:-1;
    bool read_static_variables=static_frame!=-1 || first_time || !(deformable_geometry && deformable_geometry->structures.m);
    if(read_geometry) deformable_geometry->Read(STREAM_TYPE(RW()),prefix,prefix,frame,static_frame,read_static_variables);
        
    if(read_static_variables){
        int m=deformable_geometry->structures.m;active_list.Resize(m,true,true,true);
        segmented_curve_objects.Delete_Pointers_And_Clean_Memory();segmented_curve_objects.Resize(m);
        triangulated_surface_objects.Delete_Pointers_And_Clean_Memory();triangulated_surface_objects.Resize(m);
        tetrahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();tetrahedralized_volume_objects.Resize(m);
        hexahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();hexahedralized_volume_objects.Resize(m);
        free_particles_objects.Delete_Pointers_And_Clean_Memory();free_particles_objects.Resize(m);
        free_particles_indirect_arrays.Delete_Pointers_And_Clean_Memory();free_particles_indirect_arrays.Resize(m);
        int color_map_index=15;
        for(int i=1;i<=m;i++){
            STRUCTURE<TV>* structure=deformable_geometry->structures(i);
            if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": segmented curve\n";
                segmented_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_3D<T>(*segmented_curve,OPENGL_COLOR((T).5,(T).25,0));
                segmented_curve_objects(i)->use_solid_color=false;}
            else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": triangulated surface, range = "<<ARRAYS_COMPUTATIONS::Min(triangulated_surface->mesh.elements.Flattened())<<" "<<ARRAYS_COMPUTATIONS::Max(triangulated_surface->mesh.elements.Flattened())<<"\n";
                triangulated_surface->mesh.Initialize_Segment_Mesh();
                ARRAY<OPENGL_COLOR> front_colors,back_colors;
                front_colors.Append(OPENGL_COLOR::Yellow());back_colors.Append(OPENGL_COLOR::Magenta());
                front_colors.Append(OPENGL_COLOR::Green());back_colors.Append(OPENGL_COLOR::Cyan());
                triangulated_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(*triangulated_surface,false,
                    OPENGL_MATERIAL::Metal(front_colors((i-1)%front_colors.Size()+1)),OPENGL_MATERIAL::Metal(back_colors((i-1)%back_colors.Size()+1)));}
            else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": tetrahedralized_volume, range = "<<ARRAYS_COMPUTATIONS::Min(tetrahedralized_volume->mesh.elements.Flattened())<<" "<<ARRAYS_COMPUTATIONS::Max(tetrahedralized_volume->mesh.elements.Flattened())<<"\n";
                tetrahedralized_volume_objects(i)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&tetrahedralized_volume->mesh,&(deformable_geometry->particles),
                    OPENGL_MATERIAL::Metal(OPENGL_COLOR::Red(.7f)));}
            else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": hexahedralized_volume\n";
                hexahedralized_volume_objects(i)=new OPENGL_HEXAHEDRALIZED_VOLUME<T>(&hexahedralized_volume->mesh,&(deformable_geometry->particles),
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));}
            else if(FREE_PARTICLES<TV>* fp=dynamic_cast<FREE_PARTICLES<TV>*>(structure)){
                free_particles_indirect_arrays(i)=new INDIRECT_ARRAY<ARRAY_VIEW<TV> >(deformable_geometry->particles.X,fp->nodes);
                free_particles_objects(i)=new OPENGL_FREE_PARTICLES<TV>(*deformable_geometry,*free_particles_indirect_arrays(i),color_map->Lookup(color_map_index--));}
            else{if(first_time) LOG::cout<<"object "<<i<<": object unrecognized at geometry level\n";}}}
    for(int i=1;i<=deformable_geometry->structures.m;i++){
        std::string i_string=STRING_UTILITIES::string_sprintf("%d",i);
        if(tetrahedralized_volume_objects(i)) has_tetrahedralized_volumes=true;
        if(hexahedralized_volume_objects(i)) has_hexahedralized_volumes=true;
        if(tetrahedralized_volume_objects(i)){
            std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/subset_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File<RW>(filename,tetrahedralized_volume_objects(i)->subset);
            filename=STRING_UTILITIES::string_sprintf("%s/%d/colliding_nodes_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File<RW>(filename,tetrahedralized_volume_objects(i)->subset_particles);}
        else if(hexahedralized_volume_objects(i)){
            std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/subset_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File<RW>(filename,hexahedralized_volume_objects(i)->subset);
            filename=STRING_UTILITIES::string_sprintf("%s/%d/colliding_nodes_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File<RW>(filename,hexahedralized_volume_objects(i)->subset_particles);
            filename=STRING_UTILITIES::string_sprintf("%s/%d/directions_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File<RW>(filename,hexahedralized_volume_objects(i)->vectors_at_hex_centers);}}
    if(smooth_shading){
        for(int i=1;i<=triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i)) triangulated_surface_objects(i)->Initialize_Vertex_Normals();
        for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Initialize_Vertex_Normals();}
    for(int i=1;i<=free_particles_indirect_arrays.m;i++) if(free_particles_indirect_arrays(i)){
        ARRAY_VIEW<TV> tmp(deformable_geometry->particles.X);
        free_particles_indirect_arrays(i)->array.Exchange(tmp);}

    Update_Velocity_Field();
    frame_loaded=frame;
    valid=true;
    first_time=false;

    for(number_of_segmented_curve=1;deformable_geometry->template Find_Structure<SEGMENTED_CURVE<TV>*>(number_of_segmented_curve);number_of_segmented_curve++);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_particles",prefix.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    if(draw_input) ARRAYS_COMPUTATIONS::Fill(active_list,true);
    Reinitialize();
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Set_Display_Modes_For_Geometry_Collection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Display_Modes_For_Geometry_Collection(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
        bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects) const
{
    display_triangulated_surface_objects=display_mode!=1;
    display_tetrahedralized_volume_objects=display_mode!=2;
    display_hexahedralized_volume_objects=display_mode!=2;
    display_free_particles_objects=display_mode!=1;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Display(const int in_color) const
{
    if(!draw || !valid) return;
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    bool display_triangulated_surface_objects,display_tetrahedralized_volume_objects,display_hexahedralized_volume_objects,display_free_particles_objects;
    Set_Display_Modes_For_Geometry_Collection(display_triangulated_surface_objects,display_tetrahedralized_volume_objects,display_hexahedralized_volume_objects,
            display_free_particles_objects);
    for(int i=1;i<=segmented_curve_objects.m;i++){
        if(!active_list(i)) continue;
        glPushName(i);
        if(segmented_curve_objects(i)){
            glPushName(1);
            segmented_curve_objects(i)->hide_unselected=hide_unselected;
            if(real_selection) segmented_curve_objects(i)->parent_curve=dynamic_cast<OPENGL_SEGMENTED_CURVE_3D<T>*>(real_selection->subobject);
            else segmented_curve_objects(i)->parent_curve=0;
            segmented_curve_objects(i)->Display(in_color);
            glPopName();}
        if(triangulated_surface_objects(i) && display_triangulated_surface_objects){
            glPushName(2);triangulated_surface_objects(i)->Display(in_color);glPopName();}
        if(tetrahedralized_volume_objects(i) && display_tetrahedralized_volume_objects){
            glPushName(3);tetrahedralized_volume_objects(i)->Display(in_color);glPopName();}
        if(hexahedralized_volume_objects(i) && display_hexahedralized_volume_objects){
            glPushName(3);hexahedralized_volume_objects(i)->Display(in_color);glPopName();}
        if(free_particles_objects(i) && display_free_particles_objects){glPushName(6);free_particles_objects(i)->Display(in_color);glPopName();}
        glPopName();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();

    if(selected_vertex) OPENGL_SELECTION::Draw_Highlighted_Vertex(deformable_geometry->particles.X(selected_vertex),selected_vertex);

    if(draw_velocity_vectors) velocity_field.Display(in_color);

    // Visualize relative velocity on edges
    if(display_relative_velocity_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        SEGMENTED_CURVE<TV>* segmented_curve=deformable_geometry->template Find_Structure<SEGMENTED_CURVE<TV>*>(display_relative_velocity_mode);
        for(int j=1;j<=segmented_curve->mesh.elements.m;j++){int p1=segmented_curve->mesh.elements(j)(1),p2=segmented_curve->mesh.elements(j)(2);
            TV relative_velocity=deformable_geometry->particles.V(p2)-deformable_geometry->particles.V(p1);
            TV edge_vector=deformable_geometry->particles.X(p2)-deformable_geometry->particles.X(p1);
            T edge_length=edge_vector.Magnitude();
            if(edge_length>(T)0){
                OPENGL_COLOR edge_color=color_map_relative_velocity->Lookup(TV::Dot_Product(relative_velocity,edge_vector)/edge_length);
                edge_color.Send_To_GL_Pipeline();
                OpenGL_Line(deformable_geometry->particles.X(p1),deformable_geometry->particles.X(p2));}}
        OpenGL_End();
        glPopAttrib();}
}
//#####################################################################
// Function Cycle_Display_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Cycle_Display_Mode()
{
    display_mode=(display_mode+1)%4;
}
//#####################################################################
// Function Cycle_Cutaway_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Cycle_Cutaway_Mode()
{
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Cycle_Cutaway_Mode();
}
//#####################################################################
// Function Show_Only_First
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Show_Only_First()
{
    ARRAYS_COMPUTATIONS::Fill(active_list,false);
    active_list(1)=true;
}
//#####################################################################
// Function Decrease_Cutaway_Fraction
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Decrease_Cutaway_Fraction()
{
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Decrease_Cutaway_Fraction();
}
//#####################################################################
// Function Increase_Cutaway_Fraction
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Increase_Cutaway_Fraction()
{
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Increase_Cutaway_Fraction();
}
//#####################################################################
// Function Create_One_Big_Triangulated_Surface_And_Write_To_File
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Create_One_Big_Triangulated_Surface_And_Write_To_File()
{
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=deformable_geometry->particles;TRIANGLE_MESH mesh;
    TRIANGULATED_SURFACE<T> triangulated_surface(mesh,particles);
    for(int index=1;index<=deformable_geometry->structures.m;index++){
        LOG::cout<<"Adding "<<index<<"th triangulated surface out of "<<deformable_geometry->structures.m<<" surfaces"<<std::endl;
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_geometry->structures(index));
        if(!tetrahedralized_volume) continue;
        tetrahedralized_volume->Initialize_Triangulated_Surface();
        TRIANGULATED_SURFACE<T>& individual_boundary_triangulated_surface=*tetrahedralized_volume->triangulated_surface;
        int individual_number_of_triangles=individual_boundary_triangulated_surface.mesh.elements.m;
        for(int t=1;t<=individual_number_of_triangles;t++) mesh.elements.Append(individual_boundary_triangulated_surface.mesh.elements(t));}
    mesh.number_nodes=particles.array_collection->Size();
    FILE_UTILITIES::Write_To_File<T>("one_big_tri_surface.tri",triangulated_surface);
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Use_Bounding_Box() const
{
    return draw && valid && deformable_geometry->structures.m>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    if(draw && valid && deformable_geometry->structures.m>0){
        for(int i=1;i<=segmented_curve_objects.m;i++) if(segmented_curve_objects(i))box.Enlarge_To_Include_Box(segmented_curve_objects(i)->Bounding_Box());
        for(int i=1;i<=triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i))box.Enlarge_To_Include_Box(triangulated_surface_objects(i)->Bounding_Box());
        for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i))box.Enlarge_To_Include_Box(tetrahedralized_volume_objects(i)->Bounding_Box());
        for(int i=1;i<=hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i))box.Enlarge_To_Include_Box(hexahedralized_volume_objects(i)->Bounding_Box());        
        for(int i=1;i<=free_particles_objects.m;i++) if(free_particles_objects(i)) box.Enlarge_To_Include_Box(free_particles_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Highlight_Particle
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Highlight_Particle()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Enter Particle Number: ",Highlight_Particle_Response_CB());
}
//#####################################################################
// Function Toggle_Active_Value
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Active_Value()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Enter Component Number: ",Toggle_Active_Value_Response_CB());
}
//#####################################################################
// Function Toggle_Hide_Unselected
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Hide_Unselected()
{
    if(real_selection) hide_unselected=!hide_unselected;
}
//#####################################################################
// Function Toggle_Draw_Interior
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Draw_Interior()
{
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i)) tetrahedralized_volume_objects(i)->Toggle_Boundary_Only();
    for(int i=1;i<=hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i)) hexahedralized_volume_objects(i)->Toggle_Boundary_Only();
}
//#####################################################################
// Function Toggle_Draw_Subsets
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Draw_Subsets()
{
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))
        tetrahedralized_volume_objects(i)->draw_subsets=!tetrahedralized_volume_objects(i)->draw_subsets;
    for(int i=1;i<=hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i))
        hexahedralized_volume_objects(i)->draw_subsets=!hexahedralized_volume_objects(i)->draw_subsets;
}
//#####################################################################
// Function Toggle_Use_Active_List
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Use_Active_List()
{
    if(active_list.m<1) return;
    use_active_list=!use_active_list;if(!use_active_list) ARRAYS_COMPUTATIONS::Fill(active_list,true);
    else{ARRAYS_COMPUTATIONS::Fill(active_list,false);incremented_active_object=1;active_list(incremented_active_object)=true;}
}
//#####################################################################
// Function Select_Segment
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Selection_Mode()
{
    if (!real_selection) return;
    assert(real_selection->body_selection);
    if (real_selection->body_selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_3D){
        delete real_selection->body_selection;
        real_selection->body_selection=real_selection->saved_selection;}
    else if (real_selection->body_selection->type == OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_3D){
        real_selection->saved_selection=real_selection->body_selection;
        real_selection->body_selection=((OPENGL_SEGMENTED_CURVE_3D<T>*)real_selection->subobject)->Get_Curve_Selection(
                ((OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>*)real_selection->body_selection)->index);}
    else return;
    Highlight_Selection(real_selection);
    glutPostRedisplay();
}
//#####################################################################
// Function Increment_Active_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
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
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* selection=0;
    if(buffer_size>=1){
        selection=new OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>(this);
        selection->body_index=buffer[0];selection->subobject_type=buffer[1];
        switch(selection->subobject_type){
            case 1:selection->subobject=segmented_curve_objects(buffer[0]);break;
            case 2:selection->subobject=triangulated_surface_objects(buffer[0]);break;
            case 3:selection->subobject=tetrahedralized_volume_objects(buffer[0]);break;
            case 6:selection->subobject=free_particles_objects(buffer[0]);break;
            default:PHYSBAM_FATAL_ERROR();}
        selection->body_selection=selection->subobject->Get_Selection(&buffer[2],buffer_size-2);
        if(!selection->body_selection){delete selection;selection=0;}}
    return selection;
}
//#####################################################################
// Function Set_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type!=OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D) return;
    real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type!=OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D) return;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
    if(selection->hide) active_list(real_selection->body_index)=false;
    else real_selection->subobject->Highlight_Selection(real_selection->body_selection);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Clear_Highlight()
{
    for(int i=1;i<=segmented_curve_objects.m;i++) if(segmented_curve_objects(i) && active_list(i)) segmented_curve_objects(i)->Clear_Highlight();
    for(int i=1;i<=triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i) && active_list(i)) triangulated_surface_objects(i)->Clear_Highlight();
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i)) tetrahedralized_volume_objects(i)->Clear_Highlight();
    for(int i=1;i<=hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i)) hexahedralized_volume_objects(i)->Clear_Highlight();
    for(int i=1;i<=free_particles_objects.m;i++) if(free_particles_objects(i) && active_list(i))free_particles_objects(i)->Clear_Highlight();
    if(hide_unselected) hide_unselected=false;
    real_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    if(selection && selection->type==OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D && selection->object==this){
        OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
        output_stream<<"Deformable object "<<real_selection->body_index<<std::endl;
        real_selection->subobject->Print_Selection_Info(output_stream,real_selection->body_selection);}
}
//#####################################################################
// Function Toggle_Active_Value_Response
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
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
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Highlight_Particle_Response()
{
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        int index=0;std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);sstream>>index;
        if(index>0 && index<=deformable_geometry->particles.array_collection->Size()) selected_vertex=index;}
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Update_Velocity_Field()
{
    if(draw_velocity_vectors){
        positions=deformable_geometry->particles.X;
        velocity_vectors=deformable_geometry->particles.V;}
}
//#####################################################################
// Function Cycle_Relative_Velocity_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Cycle_Relative_Velocity_Mode()
{
    display_relative_velocity_mode=(display_relative_velocity_mode+1)%number_of_segmented_curve;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T,RW>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this && invalidate_deformable_objects_selection_each_frame) delete_selection=true;
    return 0;
}
//#####################################################################
template class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<double,double>;
#endif
