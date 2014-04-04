//#####################################################################
// Copyright 2004-2008, Zhaosheng Bao, Eilene Hao, Jeong-Mo Hong, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Jerry Talton, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SEGMENTED_CURVE.h>
#include "GENERIC_RENDER_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Function Apply_General_Parameters
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Apply_General_Parameters(RENDERING_OBJECT<T>& rendering_object,PARAMETER_LIST& parameters)
{
    const T minimum_surface_roughness=(T)1e-5,surface_roughness_factor=4;
    rendering_object.small_number=surface_roughness_factor*minimum_surface_roughness;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    std::string bssrdf_shader=parameters.Get_Parameter("BSSRDF_Shader",std::string("<unknown>"));
    if(bssrdf_shader!="<unknown>") rendering_object.bssrdf_shader=dynamic_cast<RENDERING_BSSRDF_SHADER<T>*>(shaders.Get(bssrdf_shader));
#endif
}
//#####################################################################
// Function Apply_Triangulated_Surface_Parameters
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Apply_Triangulated_Surface_Parameters(RENDERING_TRIANGULATED_SURFACE<T>& rendering_triangulated_surface,PARAMETER_LIST& parameters)
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=rendering_triangulated_surface.triangulated_surface;
    rendering_triangulated_surface.add_to_spatial_partition=parameters.Get_Parameter("Spatial_Partition",true);
    if(parameters.Get_Parameter("Smooth_Normals",true)) triangulated_surface.Use_Vertex_Normals();
    rendering_triangulated_surface.add_triangles_to_acceleration_structure=parameters.Get_Parameter("Accel_By_Triangle",true);
    triangulated_surface.avoid_normal_interpolation_across_sharp_edges=parameters.Get_Parameter("Preserve_Creases",true);
    T scale=parameters.Get_Parameter("Scale",(T)1);if(scale!=1) triangulated_surface.Rescale(scale);
    if(parameters.Get_Parameter("Subdivide_Geometry",false)){triangulated_surface.Loop_Subdivide();triangulated_surface.Refresh_Auxiliary_Structures();}
    std::string texture_coordinate_filename=parameters.Get_Parameter("Texture_Coordinate_File",std::string("unknown"));
    if(texture_coordinate_filename!="unknown"){
        rendering_triangulated_surface.template Read_Texture_Coordinates<float>(texture_coordinate_filename);
        LOG::cout<<"Reading text file "<<texture_coordinate_filename<<"\n";}
    if(triangulated_surface.use_vertex_normals) triangulated_surface.Update_Vertex_Normals();
}
//#####################################################################
// Function Get_Rigid_Objects
//#####################################################################
template<class T,class RW> RIGID_GEOMETRY_COLLECTION<VECTOR<T,3> >& GENERIC_RENDER_EXAMPLE<T,RW>::
Get_Rigid_Objects(PARAMETER_LIST& parameters,int frame,std::string& rigid_body_collection_name)
{
    RIGID_GEOMETRY_COLLECTION<TV>* rigid_body_collection;
    std::string unknown="<unknown>",prefix=parameters.Get_Parameter("Prefix",unknown),particles_name=parameters.Get_Parameter("List_Name",unknown);
    if(particles_name!=unknown && rigid_body_collection_list.Get(particles_name,rigid_body_collection)){rigid_body_collection_name=particles_name;return *rigid_body_collection;}
    PHYSBAM_ASSERT(prefix!=unknown);
    std::string name=parameters.Get_Parameter("Name",unknown);
    PHYSBAM_ASSERT(name!=unknown);
    if(!rigid_body_collection_cached.Get(prefix,rigid_body_collection)){
        rigid_body_collection=new RIGID_GEOMETRY_COLLECTION<TV>(0,0);
        // bool load_implicit_surfaces;load_implicit_surfaces=parameters.Get_Parameter("Load_Implicit_Surfaces",false); // TODO: this is currently ignored
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(STREAM_TYPE(RW()),prefix,frame,*rigid_body_collection);
        rigid_body_collection_cached.Set(prefix,rigid_body_collection);
        rigid_body_collection_list.Set(name,rigid_body_collection);}
    rigid_body_collection_name=name;
    return *rigid_body_collection;
}
//#####################################################################
// Function List_Object
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
List_Object(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    std::string type=parameters.Get_Parameter("Type",std::string("Rigid_Body_List"));
    std::string name=parameters.Get_Parameter("Name",std::string("<unknown>"));
    std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
    std::string range=parameters.Get_Parameter("Range",std::string("<unknown>"));
    int triangles_per_hierarchy_group=parameters.Get_Parameter("Triangles_Per_Hierarchy_Group",0);
    int priority=parameters.Get_Parameter("Priority",0);
    bool support_transparent_overlapping_objects=parameters.Get_Parameter("Support_Transparent_Overlapping_Objects",false);
    std::string rigid_body_collection_name;
    LOG::cout<<"Object has priority "<<priority<<" and support for overlapping transparent objects is: "<<support_transparent_overlapping_objects<<std::endl;
    ARRAY<int> integer_list;
    if(name=="<unknown>"){LOG::cerr<<"No list object name specified"<<std::endl;PHYSBAM_FATAL_ERROR();}
    if(shader_name=="<unknown>"){LOG::cerr<<"No shader (volume or surface) specified for object "<<name<<std::endl;PHYSBAM_FATAL_ERROR();}
    if(!shaders.Contains(shader_name)){LOG::cerr<<"Invalid shader '"<<shader_name<<"' specified for object "<<name<<std::endl;PHYSBAM_FATAL_ERROR();}

    RENDERING_OBJECT<T> *object=0;
    if(type=="Rigid_Body_List" || type=="Rigid_Body_Instance"){
        std::string sample_locations_filename=parameters.Get_Parameter("Sample_Locations_File",std::string("unknown"));
        bool two_sided=parameters.Get_Parameter("Two_Sided",true);
        ARRAY<int> id_list;
        ARRAY<int> parents;// parsing parent data
        ARRAY<int> rigid_body_parents;

        if(type=="Rigid_Body_List"){
            std::string rigid_body_particles_file_name=parameters.Get_Parameter("Prefix",std::string("unknown"));
            std::string parents_input=parameters.Get_Parameter("Parents",std::string(""));
            STRING_UTILITIES::Parse_Integer_List(parents_input,parents);
            if(parents.m) FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_parents",rigid_body_particles_file_name.c_str(),frame),rigid_body_parents);}
        RIGID_GEOMETRY_COLLECTION<TV> &rigid_body_collection=Get_Rigid_Objects(parameters,frame,rigid_body_collection_name);

        if(range!="<unknown>"){
            STRING_UTILITIES::Parse_Integer_List(range,integer_list);
            id_list.Resize(integer_list.Size());
            for(int i=1;i<=integer_list.m;i++) id_list(i)=int(integer_list(i));}
        else if(type=="Rigid_Body_List") for(int i(1);i<=rigid_body_collection.particles.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)) id_list.Append(i);
        else PHYSBAM_FATAL_ERROR("A Range is Required for a Rigid_Body_Instance.");

        if(type=="Rigid_Body_List" && parents.m){
            for(int index=1;index<=id_list.m;index++){int id=id_list(index);
                if(rigid_body_parents(index)){
                    int parent_number=0;
                    for(int k=1;k<=parents.m;k++) if(rigid_body_parents(index)==parents(k)) parent_number++;
                    if(parent_number==0){rigid_body_collection.particles.Remove_Geometry(id);}}}}

        for(int index=1;index<=id_list.m;index++){int id=id_list(index);
            if(parents.m && rigid_body_parents(index)){
                int parent_number=0;
                for(int k=1;k<=parents.m;k++)
                    if(rigid_body_parents(index)==parents(k)) parent_number++;
                if(parent_number==0) continue;}// rigid_body_parents(index)==0 for rigid bodies non 0 for deformable for breakable bodies.
            std::string render_object=parameters.Get_Parameter("Render_Object",std::string("Null")); // object to use instead of an extra object
            if(render_object=="Null"){
                RENDERING_TRIANGULATED_SURFACE<T>* rendering_triangulated_surface=0;
                std::string file_name=rigid_body_collection_name+"_"+STRING_UTILITIES::string_sprintf("%d",id);
                if(type=="Rigid_Body_Instance" && objects.Get(file_name,object)){rendering_triangulated_surface=dynamic_cast<RENDERING_TRIANGULATED_SURFACE<T>*>(object);object->name=name;}
                if(!rendering_triangulated_surface){
                    if(!rigid_body_collection.Rigid_Geometry(id).simplicial_object->mesh.elements.m){
                        IMPLICIT_OBJECT<TV>* implicit_object=rigid_body_collection.Rigid_Geometry(id).implicit_object->object_space_implicit_object;
                        FRAME<TV> levelset_frame=FRAME<TV>();
                        while(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* transform=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(implicit_object)){
                            levelset_frame=levelset_frame**transform->transform;
                            implicit_object=transform->object_space_implicit_object;}
                        LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object);
                        if(!lio) PHYSBAM_FATAL_ERROR("Should have levelset implicit object");
                        TRIANGULATED_SURFACE<T>* surface_from_levelset=DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(lio->levelset);
                        object=rendering_triangulated_surface=new RENDERING_TRIANGULATED_SURFACE<T>(*surface_from_levelset->Create_Compact_Copy(),triangles_per_hierarchy_group);
                        rigid_body_collection.Rigid_Geometry(id).Set_Frame(rigid_body_collection.Rigid_Geometry(id).Frame()*levelset_frame);
                        delete surface_from_levelset;}
                    else object=rendering_triangulated_surface=new RENDERING_TRIANGULATED_SURFACE<T>(*rigid_body_collection.Rigid_Geometry(id).simplicial_object->Create_Compact_Copy(),triangles_per_hierarchy_group);
                    rendering_triangulated_surface->triangulated_surface.particles.Store_Velocity();
                    object->name=file_name;
                    objects.Set(object->name,object);}

                Apply_General_Parameters(*rendering_triangulated_surface,parameters);
                Apply_Triangulated_Surface_Parameters(*rendering_triangulated_surface,parameters);

                shaders.Get(shader_name,object->material_shader);
                GENERIC_RENDER_EXAMPLE<T,RW>::Add_Solid_Texture(object,parameters);
                object->Update_Transform(current_transform*rigid_body_collection.Rigid_Geometry(id).Frame().Matrix());
                object->two_sided=two_sided;
                if(sample_locations_filename!="unknown") rendering_triangulated_surface->sample_locations_file=sample_locations_filename;
                object->priority=priority;object->support_transparent_overlapping_objects=support_transparent_overlapping_objects;}
            else{
                RENDERING_OBJECT<T>* object=objects.Get(render_object);
                object->Update_Transform(current_transform*rigid_body_collection.Rigid_Geometry(id).Frame().Matrix());}
            LOG::cout<<"Processed Rigid Body "<<id<<std::endl;}}
    else if(type=="Deformable_Object"){
        DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_body_collection=*new DEFORMABLE_GEOMETRY_COLLECTION<TV>(*new GEOMETRY_PARTICLES<TV>);
        int local_frame=parameters.Get_Parameter("Frame",frame);
        bool split_object=parameters.Get_Parameter("Split_Object",false);
        std::string sample_locations_filename=parameters.Get_Parameter("Sample_Locations_File",std::string("unknown"));
        bool two_sided=parameters.Get_Parameter("Two_Sided",true);
        std::string prefix=parameters.Get_Parameter("Prefix",std::string("<unknown>"));
        std::string static_frame_prefix=parameters.Get_Parameter("Static_Frame_Prefix",prefix);
        prefix+="/";static_frame_prefix+="/";
        std::string frame_string=STRING_UTILITIES::string_sprintf("%d/",local_frame);
        int static_frame=FILE_UTILITIES::File_Exists(static_frame_prefix+frame_string+"deformable_object_structures")?frame:-1;
        std::string free_particles_geometry=parameters.Get_Parameter("Free_Particles_Geometry",std::string("Null")); // object to use instead of an extra object
        std::string free_particles_range=parameters.Get_Parameter("Free_Particles_Range",std::string("<unknown>"));
        deformable_body_collection.Read(STREAM_TYPE(RW()),prefix,static_frame_prefix,local_frame,static_frame,true);
        if(range!="<unknown>") STRING_UTILITIES::Parse_Integer_List(range,integer_list);
        else integer_list=IDENTITY_ARRAY<>(deformable_body_collection.structures.m);
        if(split_object){
            if(deformable_body_collection.structures.m!=1) PHYSBAM_NOT_IMPLEMENTED("Split_Object for more than one object");
            if(TETRAHEDRALIZED_VOLUME<T>* volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>()){
                deformable_body_collection.structures.Remove_All();
                ARRAY<int> label;volume->mesh.Identify_Edge_Connected_Components(label);
                HASHTABLE<int> component_done;
                for(int p=1;p<=label.m;p++) if(component_done.Set(label(p))){
                    TETRAHEDRALIZED_VOLUME<T>* component_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
                    for(int t=1;t<=volume->mesh.elements.m;t++) if(label.Subset(volume->mesh.elements(t)).Contains(label(p)))
                        component_volume->mesh.elements.Append(volume->mesh.elements(t));
                    component_volume->Update_Number_Nodes();
                    deformable_body_collection.Add_Structure(component_volume);}
                LOG::cout<<"Split object into "<<component_done.Size()<<" components"<<std::endl;
                delete volume;}
            else PHYSBAM_NOT_IMPLEMENTED("Split_Object for non-tetrahedralized volumes");}
        if(ARRAYS_COMPUTATIONS::Min(integer_list)<1 || ARRAYS_COMPUTATIONS::Max(integer_list)>deformable_body_collection.structures.m){
            LOG::cerr<<"Range out of bound for object "<<name<<std::endl;PHYSBAM_FATAL_ERROR();}
        for(int index=1;index<=integer_list.m;index++){
            int i=integer_list(index);
            object=0;
            RENDERING_TRIANGULATED_SURFACE<T>* surface=0;
            STRUCTURE<TV>* structure=deformable_body_collection.structures(i);
            if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(structure)){
                LOG::cout<<"object "<<i<<": segmented curve"<<std::endl;
                T thickness=parameters.Get_Parameter("Thickness",(T).01);
                RENDERING_SEGMENTED_CURVE<T> *rendering_curve=new RENDERING_SEGMENTED_CURVE<T>(*segmented_curve,thickness);
                object=rendering_curve;}
            else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
                LOG::cout<<"object "<<i<<": triangulated surface"<<std::endl;
                triangulated_surface->mesh.Initialize_Segment_Mesh();
                surface=new RENDERING_TRIANGULATED_SURFACE<T>(*triangulated_surface->Create_Compact_Copy(),triangles_per_hierarchy_group);}
            else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
                LOG::cout<<"object "<<i<<": tetrahedralized_volume"<<std::endl;
                tetrahedralized_volume->Update_Number_Nodes();
                tetrahedralized_volume->Initialize_Triangulated_Surface();
                surface=new RENDERING_TRIANGULATED_SURFACE<T>(*tetrahedralized_volume->triangulated_surface->Create_Compact_Copy(),triangles_per_hierarchy_group);}
            else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(structure)){
                LOG::cout<<"object "<<i<<": hexahedralized_volume"<<std::endl;
                hexahedralized_volume->Initialize_Triangulated_Surface();
                surface=new RENDERING_TRIANGULATED_SURFACE<T>(*hexahedralized_volume->triangulated_surface->Create_Compact_Copy(),triangles_per_hierarchy_group);}
            else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(structure)){
                ARRAY<int> free_particles_list;
                if(free_particles_range!="<unknown>") STRING_UTILITIES::Parse_Integer_List(free_particles_range,free_particles_list);
                else free_particles_list=IDENTITY_ARRAY<>(free_particles->nodes.m);
                if(free_particles_list.m!=1) PHYSBAM_FATAL_ERROR("FREE_PARTICLES only supports one particle at a time");
                RENDERING_OBJECT<T>* free_particles_object=objects.Get(free_particles_geometry);
                if(!free_particles_object){PHYSBAM_FATAL_ERROR("Rendering FREE_PARTICLES requires a valid Free_Particles_Geometry rendering object.");}
                LOG::cout<<"object "<<i<<": free particles "<<free_particles_list;
                int free_particle=free_particles->nodes(free_particles_list(1));
                TV point=deformable_body_collection.particles.X(free_particle);
                free_particles_object->transform.Add_To_Submatrix(1,4,point);
                free_particles_object->inverse_transform=free_particles_object->transform.Inverse();}
            else{LOG::cout<<"Weird object "<<i<<std::endl;PHYSBAM_FATAL_ERROR();}
            if(!object) object=surface;
            if(object){
                Apply_General_Parameters(*object,parameters);
                if(RENDERING_TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<RENDERING_TRIANGULATED_SURFACE<T>*>(object)){
                    triangulated_surface->triangulated_surface.particles.Store_Velocity();
                    Apply_Triangulated_Surface_Parameters(*triangulated_surface,parameters);}

                object->Update_Transform(current_transform);
                std::string texture_coordinate_file=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Texture_Coordinate_File%d",i),std::string("unknown"));
                if(texture_coordinate_file!="unknown"){
                    std::string filename;
                    filename=STRING_UTILITIES::string_sprintf(texture_coordinate_file.c_str(),local_frame);
                    surface->template Read_Texture_Coordinates<RW>(filename);}
                object->name=name+"_"+STRING_UTILITIES::string_sprintf("%d",i);
                object->material_shader=shaders.Get(parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Shader%d",i),shader_name));
                object->priority=priority;object->support_transparent_overlapping_objects=support_transparent_overlapping_objects;
                object->two_sided=two_sided;
                if(sample_locations_filename!="unknown") dynamic_cast<RENDERING_TRIANGULATED_SURFACE<T>&>(*object).sample_locations_file=sample_locations_filename;
                GENERIC_RENDER_EXAMPLE<T,RW>::Add_Solid_Texture(object,parameters);
                objects.Set(object->name,object);object->add_to_spatial_partition=true;}}
        LOG::cout<<"Deformable_Object '"<<name<<"' Number of Structures="<<deformable_body_collection.structures.Size()<<std::endl;}
    else if(type=="Deformable_Object_Instance"){
        std::string object_name=parameters.Get_Parameter("Object_Name",std::string("<unknown>"));
        if(object_name==std::string("<unknown>")){LOG::cout<<"Unknown Deformable Object"<<std::endl;PHYSBAM_FATAL_ERROR();}
        if(range!="<unknown>") STRING_UTILITIES::Parse_Integer_List(range,integer_list);
        else{LOG::cout<<"A Range is Required for a Deformable_Object_Instance."<<std::endl;PHYSBAM_FATAL_ERROR();}
        for(int index=1;index<=integer_list.m;index++){
            int i=integer_list(index);
            std::string structure_name=object_name+"_"+STRING_UTILITIES::string_sprintf("%d",i);
            if(!objects.Get(structure_name,object)){LOG::cerr<<"Structure "<<structure_name<<" not found.  Check original range for object "<<name<<std::endl;PHYSBAM_FATAL_ERROR();}
            object->name=name;object->material_shader=shaders.Get(shader_name);
            /*if(add_to_collisions){
                PHYSBAM_NOT_IMPLEMENTED(); // TODO: fix once triangle collisions are cleaned up
                //DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_body_collection=*deformable_body_collections[object_name];
                //deformable_body_collection.Initialize_Collision_Geometry();
                //body_list->Append(&deformable_body_collection_list.deformable_body_collections(i)->collisions);
                // deformable_body_collection_list.deformable_body_collections(i)->collisions.roughness=minimum_surface_roughness;
                object->small_number=surface_roughness_factor*minimum_surface_roughness;}*/}}
}
//#####################################################################
// Function List_Object_Compute_Acceleration_Structures
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
List_Object_Compute_Acceleration_Structures()
{
    //if(collision_body_list) collision_body_list->Update_Intersection_Acceleration_Structures(false);
}
//#####################################################################
template void GENERIC_RENDER_EXAMPLE<float,float>::List_Object(RENDER_WORLD<float>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<float,float>::List_Object_Compute_Acceleration_Structures();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void GENERIC_RENDER_EXAMPLE<double,float>::List_Object(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,double>::List_Object(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,float>::List_Object_Compute_Acceleration_Structures();
template void GENERIC_RENDER_EXAMPLE<double,double>::List_Object_Compute_Acceleration_Structures();
#endif
