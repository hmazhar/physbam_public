//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AXES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD.h>
#include <sstream>
#include <string>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D(const std::string& basedir_input,bool use_display_lists,const bool initialize_geometry)
    :OPENGL_COMPONENT("Rigid Geometry Collection"),basedir(basedir_input),use_display_lists(use_display_lists),frame_loaded(-1),valid(false),
    velocity_field(velocity_vectors,positions,OPENGL_COLOR::Cyan(),.25,true,true),
    angular_velocity_field(angular_velocity_vectors,positions,OPENGL_COLOR::Magenta(),.25,true,true),need_destroy_rigid_geometry_collection(false),one_sided(false),
    front_color_map(0),back_color_map(0)
{
    if(initialize_geometry){
        need_destroy_rigid_geometry_collection=true;
        rigid_geometry_collection=new RIGID_GEOMETRY_COLLECTION<TV>(*new RIGID_GEOMETRY_PARTICLES<TV>(),0);}
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
~OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D()
{
    if(need_destroy_rigid_geometry_collection){
        delete &(rigid_geometry_collection->particles);
        delete rigid_geometry_collection;}
    delete front_color_map;delete back_color_map;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Initialize()
{
    is_animation=true;
    show_object_names=false;
    output_positions=true;
    draw_velocity_vectors=false;
    draw_angular_velocity_vectors=false;
    draw_individual_axes=false;
    draw_triangulated_surface=true;
    draw_tetrahedralized_volume=true;
    draw_implicit_surface=false;
    read_triangulated_surface=true;
    read_implicit_surface=false;
    read_tetrahedralized_volume=false;
    draw_simplicial_object_particles=false;
    has_init_destroy_information=false;

    front_color_map=OPENGL_INDEXED_COLOR_MAP::Rigid_Body_Color_Map();
    back_color_map=OPENGL_INDEXED_COLOR_MAP::Rigid_Body_Color_Map();
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Draw_Mode(const int mode)
{
    draw_triangulated_surface=(mode&0x01)!=0;
    draw_tetrahedralized_volume=(mode&0x02)!=0;
    draw_implicit_surface=(mode&0x04)!=0;
}
//#####################################################################
// Function Get_Draw_Mode
//#####################################################################
template<class T,class RW> int OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Get_Draw_Mode() const
{
    return 0x04*(int)(draw_implicit_surface)+0x02*(int)(draw_tetrahedralized_volume)+0x01*(int)(draw_triangulated_surface);
}
//#####################################################################
// Function Read_Hints
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Read_Hints(const std::string& filename)
{
    ARRAY<OPENGL_RIGID_BODY_HINTS,int> opengl_hints;
    FILE_UTILITIES::Read_From_File<RW>(filename,opengl_hints);
    for(int i(1);i<=opengl_triangulated_surface.Size();i++) if(opengl_triangulated_surface(i) && i<=opengl_hints.Size()){
        opengl_triangulated_surface(i)->Set_Front_Material(opengl_hints(i).material);
        use_object_bounding_box(i)=opengl_hints(i).include_bounding_box;}
}
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Resize_Structures(const int size)
{
    opengl_triangulated_surface.Resize(size);
    opengl_tetrahedralized_volume.Resize(size);
    opengl_levelset.Resize(size);
    opengl_octree_levelset_surface.Resize(size);
    opengl_axes.Resize(size);
    draw_object.Resize(size);
    use_object_bounding_box.Resize(size);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame))) return;

        // TODO: currently reads in all structures, should only read in certain kinds based on read_triangulated_surface,read_implicit_surface,read_tetrahedralized_volume
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

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(opengl_triangulated_surface.Size(),rigid_geometry_collection->particles.array_collection->Size());

        std::string rigid_body_colors_file=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_colors",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(rigid_body_colors_file)) FILE_UTILITIES::Read_From_File<T>(rigid_body_colors_file,opengl_colors);
        else{opengl_colors.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(opengl_colors,OPENGL_COLOR::Cyan());}
        Resize_Structures(max_number_of_bodies);

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=1;i<=needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_geometry_collection->Is_Active(id));
            Create_Geometry(id);}
        else for(int i=1;i<=max_number_of_bodies;i++){if(rigid_geometry_collection->Is_Active(i)) Create_Geometry(i);} // TODO: can we figure out what bodies need_init

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
            if(rigid_geometry_collection->Is_Active(id)){
                Update_Geometry(id);
                RIGID_GEOMETRY<TV>& body=rigid_geometry_collection->Rigid_Geometry(id);
                IMPLICIT_OBJECT<TV>* object_space_implicit_object=body.implicit_object?body.implicit_object->object_space_implicit_object:0;
                if(body.name=="ground" || (object_space_implicit_object && typeid(*object_space_implicit_object)==typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >)))
                    Set_Object_Material(id,OPENGL_COLOR::Ground_Tan((T)1));
                else{
                    if(one_sided) Set_Object_Material(id,front_color_map->Lookup(Value(id)-1));
                    else Set_Object_Material(id,front_color_map->Lookup(Value(id)-1),back_color_map->Lookup(Value(id)-1));}}
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Reinitialize_Without_Files
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Reinitialize_Without_Files(const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::cout<<"Reinit being called\n";
        valid=false;

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(opengl_triangulated_surface.Size(),rigid_geometry_collection->particles.array_collection->Size());
        Resize_Structures(max_number_of_bodies);
        
        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
            if(rigid_geometry_collection->Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Initialize_One_Body
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Initialize_One_Body(const int body_id,const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::cout<<"Init for one body being called\n";
        valid=false;

        int max_number_of_bodies=max(opengl_triangulated_surface.Size(),rigid_geometry_collection->particles.array_collection->Size());
        // only enlarge array as we read in more geometry to memory
        opengl_colors.Resize(max_number_of_bodies);ARRAYS_COMPUTATIONS::Fill(opengl_colors,OPENGL_COLOR::Cyan());
        opengl_triangulated_surface.Resize(max_number_of_bodies);
        opengl_tetrahedralized_volume.Resize(max_number_of_bodies);
        opengl_levelset.Resize(max_number_of_bodies);
        opengl_octree_levelset_surface.Resize(max_number_of_bodies);
        //extra_components.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);
        use_object_bounding_box.Resize(max_number_of_bodies);

        // Initialize bodies which have become active
        PHYSBAM_ASSERT(rigid_geometry_collection->Is_Active(body_id));
        Create_Geometry(body_id);

        // Update active bodies / remove inactive bodies
        for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
            if(rigid_geometry_collection->Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Update_Bodies
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Update_Bodies(const bool update_arb_points)
{
    // Update active bodies / remove inactive bodies
    for(int id(1);id<=rigid_geometry_collection->particles.array_collection->Size();id++){
        if(rigid_geometry_collection->Is_Active(id)) Update_Geometry(id);
        else Destroy_Geometry(id);}
    for(int id=rigid_geometry_collection->particles.array_collection->Size()+1;id<=opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);
    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Create_Geometry(const int id)
{
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection->Rigid_Geometry(id);
    draw_object(id)=true;use_object_bounding_box(id)=true;
    if(rigid_geometry.name=="ground") use_object_bounding_box(id)=false; // don't use the ground bounding box
    if(!opengl_axes(id)){opengl_axes(id)=new OPENGL_AXES<T>();}

    // add tetrahedralized volume
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=rigid_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>();
    if(tetrahedralized_volume && !opengl_tetrahedralized_volume(id)){
        tetrahedralized_volume->mesh.Initialize_Triangle_Mesh();
        opengl_tetrahedralized_volume(id)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&tetrahedralized_volume->mesh,&tetrahedralized_volume->particles,
            OPENGL_MATERIAL::Metal(OPENGL_COLOR::Magenta(1,1)),true);
        opengl_tetrahedralized_volume(id)->Enslave_Transform_To(*opengl_axes(id));}

    // add implicit object
    if(rigid_geometry.implicit_object){
        IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_geometry.implicit_object->object_space_implicit_object;
        if(LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            if(!opengl_levelset(id)){
                opengl_levelset(id)=new OPENGL_LEVELSET_MULTIVIEW<T,RW>();
                opengl_levelset(id)->Set_Levelset(levelset_implicit_object->levelset);
                opengl_levelset(id)->Set_Slice(slice);
                opengl_levelset(id)->Generate_Triangulated_Surface();
                opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));}}
        if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object_space_implicit_object)){
            if(!opengl_levelset(id)){
                opengl_levelset(id)=new OPENGL_LEVELSET_MULTIVIEW<T,RW>();
                opengl_levelset(id)->Set_Levelset(dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object_transformed->object_space_implicit_object)->levelset);
                opengl_levelset(id)->Set_Slice(slice);
                opengl_levelset(id)->Generate_Triangulated_Surface();
                opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));
                opengl_levelset(id)->implicit_object_transform=implicit_object_transformed->transform;}}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        else if(DYADIC_IMPLICIT_OBJECT<TV>* dyadic_implicit_object=dynamic_cast<DYADIC_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            if(!opengl_octree_levelset_surface(id)){
                DUALCONTOUR_OCTREE<T> contour(&dyadic_implicit_object->levelset);
                opengl_octree_levelset_surface(id)=new OPENGL_TRIANGULATED_SURFACE<T>(*contour.Get_Triangulated_Surface());
                opengl_octree_levelset_surface(id)->Enslave_Transform_To(*opengl_axes(id));}}
#endif
        else if(!rigid_geometry.simplicial_object && !opengl_triangulated_surface(id)){
            if(typeid(*object_space_implicit_object)==typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >))
                use_object_bounding_box(id)=false; // don't use the ground bounding box
            TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(*object_space_implicit_object);
            if(surface) rigid_geometry.Add_Structure(*surface);}}

    // add triangulated surface
    if(rigid_geometry.simplicial_object && !opengl_triangulated_surface(id)){
        LOG::cout<<"name = "<<rigid_geometry.name<<std::endl;
        OPENGL_COLOR color=opengl_colors(id);
        opengl_triangulated_surface(id)=new OPENGL_TRIANGULATED_SURFACE<T>(*rigid_geometry.simplicial_object,false,OPENGL_MATERIAL::Plastic(color),
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Green()));
        opengl_triangulated_surface(id)->Enslave_Transform_To(*opengl_axes(id));
        opengl_triangulated_surface(id)->draw_particles=draw_simplicial_object_particles;}

    // add extra components
    if(opengl_tetrahedralized_volume(id)){
        std::string filename_pattern=STRING_UTILITIES::string_sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            LOG::cout<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T,RW>* component=
                new OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T,RW>(*tetrahedralized_volume,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            //extra_components(id).Append(component);
        }}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<float,3> >(rigid_geometry_collection->Rigid_Geometry(id).Frame());
    if(opengl_tetrahedralized_volume(id)){
        std::string color_map_filename=STRING_UTILITIES::string_sprintf("%s/%d/stress_map_of_tetrahedralized_volume_%d",basedir.c_str(),frame,id);
        if(FILE_UTILITIES::File_Exists(color_map_filename)){
            if(!opengl_tetrahedralized_volume(id)->color_map) opengl_tetrahedralized_volume(id)->color_map=new ARRAY<OPENGL_COLOR>;
            FILE_UTILITIES::Read_From_File<RW>(color_map_filename,*opengl_tetrahedralized_volume(id)->color_map);}
        else if(opengl_tetrahedralized_volume(id)->color_map){delete opengl_tetrahedralized_volume(id)->color_map;opengl_tetrahedralized_volume(id)->color_map=0;}}
    RIGID_GEOMETRY<TV> &rigid_geometry=rigid_geometry_collection->Rigid_Geometry(id);
    if(rigid_geometry.implicit_object && rigid_geometry.implicit_object->object_space_implicit_object->update_every_frame)
    {
        if(opengl_levelset(id))
            opengl_levelset(id)->Reset_Surface();
        if(opengl_triangulated_surface(id))
        {
            delete opengl_triangulated_surface(id);
            TRIANGULATED_SURFACE<T>* surface=rigid_geometry.simplicial_object;
            rigid_geometry.Add_Structure(*TESSELLATION::Generate_Triangles(*rigid_geometry.implicit_object->object_space_implicit_object));
            delete surface;

            OPENGL_COLOR color=opengl_colors(id);
            opengl_triangulated_surface(id)=new OPENGL_TRIANGULATED_SURFACE<T>(*rigid_geometry.simplicial_object,false,OPENGL_MATERIAL::Plastic(color),
                OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Green()));
            opengl_triangulated_surface(id)->Enslave_Transform_To(*opengl_axes(id));
            opengl_triangulated_surface(id)->draw_particles=draw_simplicial_object_particles;
        }
    }

    //for(int i=1;i<=extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Destroy_Geometry(const int id)
{
    if(draw_object.Size()<id)return; // it's possible that we try to delete geometry that we've never added because of substepping in the simulator
    delete opengl_triangulated_surface(id);opengl_triangulated_surface(id)=0;delete opengl_axes(id);opengl_axes(id)=0;
    if(opengl_tetrahedralized_volume(id)){delete opengl_tetrahedralized_volume(id)->color_map;delete opengl_tetrahedralized_volume(id);opengl_tetrahedralized_volume(id)=0;}
    delete opengl_levelset(id);opengl_levelset(id)=0;delete opengl_axes(id);opengl_axes(id)=0;
    //extra_components(id).Delete_Pointers_And_Clean_Memory();
    draw_object(id)=false;
}
//#####################################################################
// Function Initialize_Display_Lists
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Initialize_Display_Lists()
{
#if 0
    cout<<"Initializing display lists..."<<endl;
    int i;
    // First create display lists for original triangulated surface
    for(i=1;i<=rigid_body_list.triangulated_surface_list.preallocated.m;i++)
        if(!rigid_body_list.triangulated_surface_list.preallocated(i)){
            opengl_triangulated_surface(i)->Create_Display_List();}
    for(i=1;i<=rigid_body_list.triangulated_surface_list.preallocated.m;i++){
        int first_instance=rigid_body_list.triangulated_surface_list.preallocated(i);
        if(first_instance>0){
            int first_instance_id=opengl_triangulated_surface(first_instance)->Get_Display_List_Id();
            opengl_triangulated_surface(i)->Use_Display_List(first_instance_id);}}
#endif
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Update_Object_Labels()
{
    int number_of_drawn_bodies=draw_object.Count_Matches(true);
    if(draw_velocity_vectors || draw_angular_velocity_vectors) positions.Resize(number_of_drawn_bodies);
    if(draw_velocity_vectors) velocity_vectors.Resize(number_of_drawn_bodies);
    if(draw_angular_velocity_vectors) angular_velocity_vectors.Resize(number_of_drawn_bodies);

    int idx=0;
    for(int i(1);i<=rigid_geometry_collection->particles.array_collection->Size();i++) if(draw_object(i)){
        if(draw_velocity_vectors || draw_angular_velocity_vectors){idx++;
            positions(idx)=rigid_geometry_collection->particles.X(i);
            if(draw_velocity_vectors)
                velocity_vectors(idx)=rigid_geometry_collection->particles.twist(i).linear;
            if(draw_angular_velocity_vectors){
                //rigid_geometry_collection->Rigid_Geometry(i).Update_Angular_Velocity();
                angular_velocity_vectors(idx)=rigid_geometry_collection->particles.twist(i).angular;}}
        if(opengl_triangulated_surface(i)){
            if(output_positions)
                opengl_triangulated_surface(i)->Set_Name(STRING_UTILITIES::string_sprintf("%s <%.3f %.3f %.3f>",rigid_geometry_collection->Rigid_Geometry(i).name.c_str(),rigid_geometry_collection->particles.X(i).x,rigid_geometry_collection->particles.X(i).y,rigid_geometry_collection->particles.X(i).z));
            else opengl_triangulated_surface(i)->Set_Name(rigid_geometry_collection->Rigid_Geometry(i).name);}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    if(draw_input) ARRAYS_COMPUTATIONS::Fill(draw_object,true);
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Display(const int in_color) const
{
    if(!draw) return;
    GLint mode;glGetIntegerv(GL_RENDER_MODE,&mode);
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}
    if(draw_triangulated_surface){
        glPushName(1);
        for(int i(1);i<=opengl_triangulated_surface.Size();i++) if(draw_object(i) && opengl_triangulated_surface(i)){
            glPushName(Value(i));opengl_triangulated_surface(i)->Display(in_color);glPopName();}
        glPopName();}
    if(draw_tetrahedralized_volume){
        glPushName(2);
        for(int i(1);i<=opengl_tetrahedralized_volume.Size();i++) if(draw_object(i) && opengl_tetrahedralized_volume(i)){
            glPushName(Value(i));opengl_tetrahedralized_volume(i)->Display(in_color);glPopName();}
        glPopName();}
    if(draw_implicit_surface){
        glPushName(3);
        int levelset_count=0;
        for(int i(1);i<=opengl_levelset.Size();i++) if(draw_object(i)){
            glPushName(Value(i));
            if(opengl_levelset(i)){
                if(++levelset_count>50){
                    PHYSBAM_WARNING("Refusing to draw more than 10 levelsets to save memory.");
                    OPENGL_WORLD::Singleton()->Add_String("WARNING: Refusing to draw more than 10 levelsets to save memory.");
                    break;}
                opengl_levelset(i)->Update();
                opengl_levelset(i)->Display(in_color);}
            if(opengl_octree_levelset_surface(i)) opengl_octree_levelset_surface(i)->Display(in_color);
            glPopName();}
        glPopName();}
    if(draw_individual_axes)
        for(int i(1);i<=opengl_axes.Size();i++){
            if(draw_object(i) && opengl_axes(i)){
                opengl_axes(i)->box.max_corner.x=opengl_axes(i)->box.max_corner.y=opengl_axes(i)->box.max_corner.z=2*rigid_geometry_collection->Rigid_Geometry(i).Object_Space_Bounding_Box().Edge_Lengths().Min();
                opengl_axes(i)->Display(in_color);}}
    if(draw_velocity_vectors) velocity_field.Display(in_color);
    if(draw_angular_velocity_vectors) angular_velocity_field.Display(in_color);

    if(show_object_names){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glColor3f(1,1,1);
        for(int i(1);i<=opengl_triangulated_surface.Size();i++)
            if(draw_object(i) && rigid_geometry_collection->Rigid_Geometry(i).name.length())
                OpenGL_String(rigid_geometry_collection->particles.X(i),rigid_geometry_collection->Rigid_Geometry(i).name);
        glPopAttrib();}
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i(1);i<=opengl_triangulated_surface.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return draw && num_drawable_objects>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    if(draw)
        for(int i(1);i<=opengl_triangulated_surface.Size();i++)
            if(rigid_geometry_collection->Is_Active(i) && rigid_geometry_collection->Rigid_Geometry(i).name!="ground")
                if(draw_object(i) && use_object_bounding_box(i) && opengl_triangulated_surface(i))
                    box=RANGE<VECTOR<float,3> >::Combine(box,opengl_triangulated_surface(i)->Bounding_Box());
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size>=2){
        int body_id(buffer[1]);
        OPENGL_SELECTION* body_selection=0;
        if(buffer[0]==1){ // segmented curve
            PHYSBAM_ASSERT(opengl_triangulated_surface(body_id));
            body_selection=opengl_triangulated_surface(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
        else if(buffer[0]==2){ // tetrahedralized_volume
            PHYSBAM_ASSERT(opengl_tetrahedralized_volume(body_id));
            body_selection=opengl_tetrahedralized_volume(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
        if(body_selection) selection=new OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>(this,body_id,body_selection);}

    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*)selection;
        if(selection->hide) draw_object(real_selection->body_id)=false;
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX ||
           real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT ||
           real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE){ // triangulated surface
            if(opengl_triangulated_surface(real_selection->body_id)) // might have become inactive
                opengl_triangulated_surface(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){ // tetrahedralized volume
            if(opengl_tetrahedralized_volume(real_selection->body_id)) // might have become inactive
                opengl_tetrahedralized_volume(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Clear_Highlight()
{
    for(int i(1);i<=opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Clear_Highlight();
    for(int i(1);i<=opengl_tetrahedralized_volume.Size();i++)
        if(opengl_tetrahedralized_volume(i)) opengl_tetrahedralized_volume(i)->Clear_Highlight();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(!selection || selection->object!=this) return;

    if(selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*)selection;

        output_stream<<"Rigid Geometry "<<real_selection->body_id<<std::endl;

        // handle case of body no longer being active
        if(!rigid_geometry_collection->Is_Active(real_selection->body_id) || !opengl_triangulated_surface(real_selection->body_id)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_GEOMETRY<TV> *body=&rigid_geometry_collection->Rigid_Geometry(real_selection->body_id);

        if(!body->name.empty()){output_stream<<"Name ="<<body->name<<std::endl;}
        Read_Write<ARRAY_COLLECTION,RW>::Print(output_stream,*rigid_geometry_collection->particles.array_collection,real_selection->body_id);

        MATRIX<T,4> body_transform=body->Frame().Matrix();

        if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX ||
            real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT ||
            real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE){ // triangulated surface
            if(opengl_triangulated_surface(real_selection->body_id))
                opengl_triangulated_surface(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);
            if(real_selection->body_selection->type==OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX){
                VECTOR<T,3> position=body->simplicial_object->particles.X(((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>*)real_selection->body_selection)->index);
                output_stream<<"Pointwise velocity = "<<body->Pointwise_Object_Velocity(body->World_Space_Point(position))<<std::endl;}}
        else if(real_selection->body_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
            if(opengl_tetrahedralized_volume(real_selection->body_id))
                opengl_tetrahedralized_volume(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);}}
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Turn_Smooth_Shading_On()
{
    for(int i(1);i<=opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Turn_Smooth_Shading_Off()
{
    for(int i(1);i<=opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Turn_Smooth_Shading_Off();
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Slice_Has_Changed()
{
    for(int i(1);i<=opengl_levelset.Size();i++){
        if(opengl_levelset(i)) opengl_levelset(i)->Set_Slice(slice);
        if(opengl_levelset(i)) opengl_levelset(i)->Slice_Has_Changed();}
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Draw_Object(int i,bool draw_it)
{
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Object_Material
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Object_Material(int i,const OPENGL_MATERIAL &front_material_input)
{
    if(!opengl_triangulated_surface(i)) return;
    opengl_triangulated_surface(i)->Set_Front_Material(front_material_input);
    opengl_triangulated_surface(i)->Set_Two_Sided(false);
}
//#####################################################################
// Function Set_Object_Material
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Object_Material(int i,const OPENGL_MATERIAL &front_material_input,const OPENGL_MATERIAL &back_material_input)
{
    if(!opengl_triangulated_surface(i)) return;
    opengl_triangulated_surface(i)->Set_Front_Material(front_material_input);
    opengl_triangulated_surface(i)->Set_Back_Material(back_material_input);
    opengl_triangulated_surface(i)->Set_Two_Sided(true);
}
//#####################################################################
// Function Set_Use_Bounding_Box
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
    angular_velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Angular_Velocity_Vectors
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Angular_Velocity_Vectors()
{
    draw_angular_velocity_vectors=!draw_angular_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Individual_Axes
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Individual_Axes()
{
    draw_individual_axes=!draw_individual_axes;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
    angular_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
    angular_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Draw_Values
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Draw_Values()
{
    velocity_field.draw_value=!velocity_field.draw_value;
    angular_velocity_field.draw_value=!angular_velocity_field.draw_value;
}
//#####################################################################
// Function Toggle_Draw_Values
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_One_Sided()
{
    one_sided=!one_sided;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Particles
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Draw_Particles()
{
    draw_simplicial_object_particles=!draw_simplicial_object_particles;
    for(int i=1;i<=opengl_triangulated_surface.m;i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->draw_particles=draw_simplicial_object_particles;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Toggle_Draw_Mode()
{
    int mode=Get_Draw_Mode();
    mode=(mode+1)%8;
    Set_Draw_Mode(mode);
}
//#####################################################################
// Function Turn_Off_Individual_Smooth_Shading_Prompt
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Turn_Off_Individual_Smooth_Shading_Prompt()
{
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        int object_id;
        STRING_UTILITIES::String_To_Value(OPENGL_WORLD::Singleton()->prompt_response,object_id);
        if(1<=object_id && object_id<=rigid_geometry_collection->particles.array_collection->Size() && opengl_triangulated_surface(object_id))
            opengl_triangulated_surface(object_id)->Turn_Smooth_Shading_Off();}
}
//#####################################################################
// Function Turn_Off_Individual_Smooth_Shading
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Turn_Off_Individual_Smooth_Shading()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Turn off smooth shading for object: ",Turn_Off_Individual_Smooth_Shading_Prompt_CB(),"");
}
//#####################################################################
// Function Manipulate_Individual_Body_Prompt
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Manipulate_Individual_Body_Prompt()
{
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        int object_id;
        std::string command;
        std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);
        sstream>>command;
        sstream>>object_id;
        if(1<=object_id && object_id<=rigid_geometry_collection->particles.array_collection->Size() && opengl_triangulated_surface(object_id)){
            if(command=="s"){
                VECTOR<T,3> scale;
                sstream>>scale;
                LOG::cout<<"Scaling object "<<object_id<<" by "<<scale<<std::endl;
                opengl_triangulated_surface(object_id)->Rescale(scale.x,scale.y,scale.z);}}}
}
//#####################################################################
// Function Manipulate_Individual_Body
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T,RW>::
Manipulate_Individual_Body()
{
    OPENGL_WORLD::Singleton()->Prompt_User("Manipulate: ",Manipulate_Individual_Body_Prompt_CB(),"");
}
//#####################################################################
// Selection Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<double,double>;
#endif
