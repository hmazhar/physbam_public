//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERIC_RENDER_EXAMPLE
//#####################################################################
#ifndef __GENERIC_RENDER_EXAMPLE__
#define __GENERIC_RENDER_EXAMPLE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/GENERIC_PARSER.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_GRID_ACCELERATOR.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
#include <list>
#include "../RAY_TRACING_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW>
class GENERIC_RENDER_EXAMPLE:public RAY_TRACING_EXAMPLE<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RAY_TRACING_EXAMPLE<T>::output_filename;using RAY_TRACING_EXAMPLE<T>::keep_old_files;
    using RAY_TRACING_EXAMPLE<T>::alpha_filename;
    using RAY_TRACING_EXAMPLE<T>::gamma_correction;using RAY_TRACING_EXAMPLE<T>::clipping_region;using RAY_TRACING_EXAMPLE<T>::frame;

    MATRIX<T,4> current_transform;                    // list of saved transformations
    STACK<MATRIX<T,4> > transform_stack;              // current transformation for objects/lights
    std::list<RENDERING_LIGHT<T>*> lights;              // all lights in scene
    HASHTABLE<std::string,MATERIAL_SHADER<T>*> shaders;  // all shaders in scene
    HASHTABLE<std::string,VOLUMETRIC_SHADER<T>*> volume_shaders;  // all volumetric shaders in scene
    HASHTABLE<std::string,RENDERING_OBJECT<T>*> objects; // all render objects in scene
    HASHTABLE<std::string,DEFORMABLE_GEOMETRY_COLLECTION<TV>*> deformable_body_collections;
    HASHTABLE<std::string,RIGID_GEOMETRY_COLLECTION<TV>*> rigid_body_collection_list;
    HASHTABLE<std::string,RIGID_GEOMETRY_COLLECTION<TV>*> rigid_body_collection_cached;
    RENDERING_UNIFORM_GRID_ACCELERATOR<T> accelerator;
    bool use_spatial_partition;

    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* body_list;
    COLLISION_GEOMETRY_COLLECTION<TV> collision_body_list_dummy;
    GENERIC_PARSER<T>* parser;

    GENERIC_RENDER_EXAMPLE(const std::string& filename_input,const int frame_input)
        :RAY_TRACING_EXAMPLE<T>(),use_spatial_partition(false),body_list(0)
    {
        Initialize_Geometry_Particle();
        Initialize_Read_Write_Structures();
        frame=frame_input;parser=new GENERIC_PARSER<T>(filename_input,frame);
    }

    std::string Animated_Filename(const std::string& filename_input,const int substitution_frame)
    {return STRING_UTILITIES::string_sprintf(filename_input.c_str(),substitution_frame);}

//#####################################################################
    void Initialize_Scene(RENDER_WORLD<T>& world,const int frame) PHYSBAM_OVERRIDE;
    void Camera(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Get_Camera_Frame(const int camera_frame,const std::string& motion_filename,TV& current_location,TV& current_look_at,TV& current_pseudo_up);
    void Options(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Light(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Material(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Volume_Material(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Transform(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void Object(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void List_Object(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters);
    void List_Object_Compute_Acceleration_Structures();
    static RENDERING_TRIANGULATED_SURFACE<T>* Get_New_Render_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool smooth_normals,const bool preserve_creases);
    static void Add_Solid_Texture(RENDERING_OBJECT<T>* object,PARAMETER_LIST& parameters);
    void Apply_General_Parameters(RENDERING_OBJECT<T>& rendering_object,PARAMETER_LIST& parameters);
    void Apply_Triangulated_Surface_Parameters(RENDERING_TRIANGULATED_SURFACE<T>& rendering_triangulated_surface,PARAMETER_LIST& parameters);
    RIGID_GEOMETRY_COLLECTION<TV>& Get_Rigid_Objects(PARAMETER_LIST& parameters,int frame,std::string& rigid_body_collection_name);
//#####################################################################
};
}
#endif
