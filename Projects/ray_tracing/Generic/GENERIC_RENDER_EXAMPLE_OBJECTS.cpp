//#####################################################################
// Copyright 2004-2008, Zhaosheng Bao, Eilene Hao, Jeong-Mo Hong, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/GRADIENT_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/EXTRAPOLATION_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_RLE_Level_Sets/READ_WRITE_LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_CYLINDER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_IMPLICIT_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OCTREE_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SEGMENTED_CURVE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SHOCKS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SPHERE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_WALL.h>
#include "GENERIC_RENDER_EXAMPLE.h"
#include "RLE_IMPLICIT_SURFACE_TWO_LEVEL.h"
using namespace PhysBAM;
//#####################################################################
// Object
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Object(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    std::string type=parameters.Get_Parameter("Type",std::string("Null"));
    std::string name=parameters.Get_Parameter("Name",std::string("Name"));
    int priority=parameters.Get_Parameter("Priority",0);
    bool support_transparent_overlapping_objects=parameters.Get_Parameter("Support_Transparent_Overlapping_Objects",false);
    LOG::cout<<"Object has priority "<<priority<<" and support for overlapping transparent objects is: "<< support_transparent_overlapping_objects<<std::endl;
    RENDERING_OBJECT<T> *object=0;
    ARRAY<RENDERING_OBJECT<T>*> object_group;
    LOG::cout<<"---=-=-=-=Parsing object..."<<std::endl;
    if(type=="Sphere"){
        TV position=parameters.Get_Parameter("Position",TV(0,0,0));
        T radius=parameters.Get_Parameter("Radius",(T)1);
        object=new RENDERING_SPHERE<T>(position,radius);object->add_to_spatial_partition=true;
        LOG::cout<<"Object '"<<name<<"'Sphere radius="<<radius<<" position="<<position<<std::endl;}
    else if(type=="Plane"){
        LOG::cout<<"Plane"<<std::endl;
        TV position=parameters.Get_Parameter("Position",TV(0,0,0));
        TV normal=parameters.Get_Parameter("Normal",TV(0,1,0));
        RENDERING_PLANE<T>* plane=new RENDERING_PLANE<T>(normal,position);
        object=plane;
        plane->texture_vector1=parameters.Get_Parameter("Texture_Vector1",TV(1,0,0));
        plane->texture_vector2=parameters.Get_Parameter("Texture_Vector2",TV(0,0,1));}
    else if(type=="Box"){
        LOG::cout<<"Box"<<std::endl;
        T xmin=parameters.Get_Parameter("Xmin",(T)-1);T xmax=parameters.Get_Parameter("Xmax",(T)1);
        T ymin=parameters.Get_Parameter("Ymin",(T)-1);T ymax=parameters.Get_Parameter("Ymax",(T)1);
        T zmin=parameters.Get_Parameter("Zmin",(T)-1);T zmax=parameters.Get_Parameter("Zmax",(T)1);
        object=new RENDERING_BOX<T>(xmin,xmax,ymin,ymax,zmin,zmax);}
     else if(type=="Cylinder"){
        TV X1=parameters.Get_Parameter("X1",TV(0,-1,0));
        TV X2=parameters.Get_Parameter("X2",TV(0,1,0));
        T radius=parameters.Get_Parameter("Radius",(T)1);
        LOG::cout<<"Cylinder principal points X1="<<X1<<" X2="<<X2<<" radius="<<radius<<std::endl;
        object=new RENDERING_CYLINDER<T>(X1,X2,radius);}
    else if(type=="Triangle"){
        LOG::cout<<"Triangle"<<std::endl;
        TV x1=parameters.Get_Parameter("X1",TV(1,0,0));
        TV x2=parameters.Get_Parameter("X2",TV(0,1,0));
        TV x3=parameters.Get_Parameter("X3",TV(0,0,1));
        object=new RENDERING_TRIANGLE<T>(x1,x2,x3);}
    else if(type=="Wall"){
        LOG::cout<<"Wall"<<std::endl;
        std::string control=parameters.Get_Parameter("Control",std::string("Null"));
        bool xmin_wall=parameters.Get_Parameter("Show_Xmin", true); bool xmax_wall=parameters.Get_Parameter("Show_Xmax",true);
        bool ymin_wall=parameters.Get_Parameter("Show_Ymin", true); bool ymax_wall=parameters.Get_Parameter("Show_Ymax",true);
        bool zmin_wall=parameters.Get_Parameter("Show_Zmin", true); bool zmax_wall=parameters.Get_Parameter("Show_Zmax",true);
        bool just_frame=parameters.Get_Parameter("Just_Frame", false);
        BOX<TV>* box=0;
        if (control=="Null"||control=="Free"){
            T xmin=parameters.Get_Parameter("Xmin",(T)-1);T xmax=parameters.Get_Parameter("Xmax",(T)1);
            T ymin=parameters.Get_Parameter("Ymin",(T)-1);T ymax=parameters.Get_Parameter("Ymax",(T)1);
            T zmin=parameters.Get_Parameter("Zmin",(T)-1);T zmax=parameters.Get_Parameter("Zmax",(T)1);
            box=new BOX<TV>(xmin,xmax,ymin,ymax,zmin,zmax);}
        else if(control=="Object"){ // any object has a bounding box so just use it...
            std::string object_name=parameters.Get_Parameter("Object",std::string("Null"));
            RENDERING_OBJECT<T>* object=0;
            if(!objects.Get(object_name,object)){LOG::cout<<"Unknown rendering object '"<<object_name<<"'specified for Wall"<<std::endl;exit(1);}
            box=new BOX<TV>(object->World_Space_Bounding_Box());}
        T shrink=parameters.Get_Parameter("Shrink", (T)0);
        box->Change_Size(-shrink);
        T thickness=parameters.Get_Parameter("Thickness",(T)0);
        if(thickness){
            // thick walls are just a group of rendering boxes
            BOX<TV> outer_box=*box;
            T perturb_value=(T)0.0005;
            if(xmin_wall){
                outer_box.min_corner.x=box->min_corner.x-thickness;outer_box.max_corner.x=box->min_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z-thickness;outer_box.max_corner.z=box->min_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->min_corner.x-thickness;outer_box.max_corner.x=box->min_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->max_corner.z-thickness;outer_box.max_corner.z=box->max_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->min_corner.x-thickness;outer_box.max_corner.x=box->min_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->min_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->min_corner.x-thickness;outer_box.max_corner.x=box->min_corner.x+thickness;
                outer_box.min_corner.y=box->max_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                if(!just_frame){
                  outer_box.min_corner.x=box->min_corner.x-thickness;outer_box.max_corner.x=box->min_corner.x+thickness;
                  outer_box.min_corner.y=box->min_corner.y+thickness-perturb_value;outer_box.max_corner.y=box->max_corner.y-thickness+perturb_value;
                  outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
            if(xmax_wall){
                outer_box.min_corner.x=box->max_corner.x-thickness;outer_box.max_corner.x=box->max_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z-thickness;outer_box.max_corner.z=box->min_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->max_corner.x-thickness;outer_box.max_corner.x=box->max_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->max_corner.z-thickness;outer_box.max_corner.z=box->max_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->max_corner.x-thickness;outer_box.max_corner.x=box->max_corner.x+thickness;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->min_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->max_corner.x-thickness;outer_box.max_corner.x=box->max_corner.x+thickness;
                outer_box.min_corner.y=box->max_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                if(!just_frame){
                  outer_box.min_corner.x=box->max_corner.x-thickness;outer_box.max_corner.x=box->max_corner.x+thickness;
                  outer_box.min_corner.y=box->min_corner.y+thickness-perturb_value;outer_box.max_corner.y=box->max_corner.y-thickness+perturb_value;
                  outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
            if(ymin_wall){
                if(!just_frame){
                  outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                  outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->min_corner.y+thickness;
                  outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
            if(ymax_wall){
                if(!just_frame){
                  outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                  outer_box.min_corner.y=box->max_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                  outer_box.min_corner.z=box->min_corner.z+thickness-perturb_value;outer_box.max_corner.z=box->max_corner.z-thickness+perturb_value;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
            if(zmin_wall){
                outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                outer_box.min_corner.y=box->max_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z-thickness;outer_box.max_corner.z=box->min_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->min_corner.y+thickness;
                outer_box.min_corner.z=box->min_corner.z-thickness;outer_box.max_corner.z=box->min_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                if(!just_frame){
                  outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                  outer_box.min_corner.y=box->min_corner.y+thickness-perturb_value;outer_box.max_corner.y=box->max_corner.y-thickness+perturb_value;
                  outer_box.min_corner.z=box->min_corner.z-thickness;outer_box.max_corner.z=box->min_corner.z+thickness;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
            if(zmax_wall){
                outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                outer_box.min_corner.y=box->max_corner.y-thickness;outer_box.max_corner.y=box->max_corner.y+thickness;
                outer_box.min_corner.z=box->max_corner.z-thickness;outer_box.max_corner.z=box->max_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                outer_box.min_corner.y=box->min_corner.y-thickness;outer_box.max_corner.y=box->min_corner.y+thickness;
                outer_box.min_corner.z=box->max_corner.z-thickness;outer_box.max_corner.z=box->max_corner.z+thickness;
                object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));

                if(!just_frame){
                  outer_box.min_corner.x=box->min_corner.x+thickness-perturb_value;outer_box.max_corner.x=box->max_corner.x-thickness+perturb_value;
                  outer_box.min_corner.y=box->min_corner.y+thickness-perturb_value;outer_box.max_corner.y=box->max_corner.y-thickness+perturb_value;
                  outer_box.min_corner.z=box->max_corner.z-thickness;outer_box.max_corner.z=box->max_corner.z+thickness;
                  object_group.Append(new RENDERING_BOX<T>(outer_box.min_corner.x,outer_box.max_corner.x,outer_box.min_corner.y,outer_box.max_corner.y,outer_box.min_corner.z,outer_box.max_corner.z));
                }
            }
        }
        else{
            T texture_scale=parameters.Get_Parameter("Texture_Scale", (T)1);
            object=new RENDERING_WALL<T>(*box, xmin_wall, xmax_wall, ymin_wall, ymax_wall, zmin_wall, zmax_wall,texture_scale);}}
    else if(type=="Segmented_Curve"){
        LOG::cout<<"Segmented Curve"<<std::endl;
        //// Get parameters
        std::string filename=parameters.Get_Parameter("Filename",std::string("unknown"));
        T thickness=parameters.Get_Parameter("Thickness",T(0.05));
        SEGMENTED_CURVE<TV>* curve;
        std::string substituted_filename=Animated_Filename(filename,frame);
        LOG::cout<<"Triangulated Surface Filename="<<filename<<" trying to read "<<substituted_filename<<std::endl;
        FILE_UTILITIES::Create_From_File<RW>(substituted_filename, curve);
        RENDERING_SEGMENTED_CURVE<T>* rendering_curve=new RENDERING_SEGMENTED_CURVE<T>(*curve, thickness);
        object=rendering_curve;
        object->add_to_spatial_partition=true;}
    else if(type=="Triangulated_Surface"){
        LOG::cout<<"Triangulated Surface"<<std::endl;
        //// Get parameters
        std::string filename=parameters.Get_Parameter("Filename",std::string("unknown"));
        T scale=parameters.Get_Parameter("Scale",T(1.0));
        bool two_sided=parameters.Get_Parameter("Two_Sided",true);
        bool smooth_normals=parameters.Get_Parameter("Smooth_Normals",true);
        bool preserve_creases=parameters.Get_Parameter("Preserve_Creases",true);
        T smooth_normals_threshold=parameters.Get_Parameter("Preserve_Creases_Threshold",T(0.1));
        std::string bssrdf_shader=parameters.Get_Parameter("BSSRDF_Shader",std::string("<unknown>"));
        bool flip_normal=parameters.Get_Parameter("Flip_Normal",false);
        std::string texture_coordinate_filename=parameters.Get_Parameter("Texture_Coordinate_File",std::string("unknown"));
        int num_subdivisions=parameters.Get_Parameter("Times_To_Subdivide",0);
        std::string bump_map_filename=parameters.Get_Parameter("Bump_Map_File",std::string("unknown"));
        T perturb_factor=parameters.Get_Parameter("Perturb_Factor",T(0.1));
        T perturb_power=parameters.Get_Parameter("Perturb_Power",T(1.0));
        int triangles_per_hierarchy_group=parameters.Get_Parameter("Triangles_Per_Hierarchy_Group",0);
        std::string sample_locations_filename=parameters.Get_Parameter("Sample_Locations_File",std::string("unknown"));
        //// create triangulated surface
        // substitute the frame in
        std::string substituted_filename=Animated_Filename(filename,frame);
        LOG::cout<<"Triangulated Surface Filename="<<filename<<" trying to read "<<substituted_filename<<std::endl;
        TRIANGULATED_SURFACE<T>* surface;FILE_UTILITIES::Create_From_File<RW>(substituted_filename,surface);
        std::cout<<"    Contains "<<surface->mesh.number_nodes<<" nodes "<<surface->mesh.elements.m<<" triangles"<<std::endl;
        surface->particles.array_collection->Resize(surface->mesh.number_nodes);
        surface->Rescale(scale);
        surface->normal_variance_threshold=smooth_normals_threshold;
        surface->avoid_normal_interpolation_across_sharp_edges=preserve_creases;
        if(smooth_normals)surface->Use_Vertex_Normals();
        for (int i=0;i<num_subdivisions;++i)surface->Linearly_Subdivide();
        RENDERING_TRIANGULATED_SURFACE<T>* triangulated_object=new RENDERING_TRIANGULATED_SURFACE<T>(*surface,triangles_per_hierarchy_group);
        std::cout<<"     two_sided="<<two_sided<<std::endl;
        triangulated_object->two_sided=two_sided;
        if(texture_coordinate_filename!="unknown"){
            triangulated_object->template Read_Texture_Coordinates<float>(texture_coordinate_filename.c_str());
            T texture_scaling_factor=parameters.Get_Parameter("Texture_Scaling_Factor",(T)0.0);
            if(texture_scaling_factor>0)triangulated_object->Rescale_Texture_Coordinates(texture_scaling_factor);}
        if(bump_map_filename!="unknown"){
            triangulated_object->Initialize_Bump_Map(bump_map_filename);
            triangulated_object->Do_Displacement_Map_Per_Vertex(perturb_factor,perturb_power);}
        if(sample_locations_filename!="unknown") triangulated_object->sample_locations_file=sample_locations_filename;
        object=triangulated_object;object->add_to_spatial_partition=true;
        if(bssrdf_shader!="<unknown>") object->bssrdf_shader=(RENDERING_BSSRDF_SHADER<T>*)shaders.Get(bssrdf_shader);
        object->flip_normal=flip_normal;}
    else if(type=="Voxel_Data"){
        LOG::cout<<"Voxel"<<std::endl;
        bool precompute=parameters.Get_Parameter("Precompute_Single_Scattering",false);
        std::string grid_filename=parameters.Get_Parameter("Grid_Filename",std::string("<unknown>"));
        std::string coarse_grid_filename=parameters.Get_Parameter("Coarse_Grid_Filename",std::string("<unknown>"));
        std::string density_filename=parameters.Get_Parameter("Density_Filename",std::string("<unknown>"));
        std::string density_fraction_filename=parameters.Get_Parameter("Density_Fraction_Filename",std::string("<unknown>"));
        std::string temperature_filename=parameters.Get_Parameter("Temperature_Filename",std::string("<unknown>"));

        bool use_density_gradient=parameters.Get_Parameter("Use_Density_Gradient",false);
        bool use_cubic_for_density=parameters.Get_Parameter("Use_Density_Cubic",false);
        T density_scale=parameters.Get_Parameter("Density_Scale",(T)1);
        T density_offset=parameters.Get_Parameter("Density_Offset",(T)0);
        bool clamp_low_density=parameters.Get_Parameter("Clamp_Low_Density",false);
        T density_lowest=parameters.Get_Parameter("Density_Lowest",(T)0);
        bool clamp_high_density=parameters.Get_Parameter("Clamp_High_Density",false);
        T density_highest=parameters.Get_Parameter("Density_Highest",(T)0);

        T temperature_scale=parameters.Get_Parameter("Temperature_Scale",(T)1);
        T temperature_offset=parameters.Get_Parameter("Temperature_Offset",(T)0);
        bool clamp_low_temperature=parameters.Get_Parameter("Clamp_Low_Temperature",false);
        T temperature_lowest=parameters.Get_Parameter("Temperature_Lowest",(T)0);
        bool clamp_high_temperature=parameters.Get_Parameter("Clamp_High_Temperature",false);
        T temperature_highest=parameters.Get_Parameter("Temperature_Highest",(T)0);

        bool collision_aware=parameters.Get_Parameter("Use_Collision_Aware_Interpolation",false);
        T volume_step=parameters.Get_Parameter("Volume_Step",(T)0.1);
        // always read densities but sometimes read temperatures
        GRID<TV>* grid=new GRID<TV>;ARRAY<T,VECTOR<int,3> >* density_data=new ARRAY<T,VECTOR<int,3> >;
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(grid_filename,frame),*grid);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(density_filename,frame),*density_data);
        if(use_density_gradient){
            LOG::cout<<"Using density gradient"<<std::endl;
            ARRAY<T,VECTOR<int,3> >* density_gradient_data=new ARRAY<T,VECTOR<int,3> >(*density_data);
            GRADIENT::Compute_Magnitude(*grid,0,*density_data,*density_gradient_data);
            delete density_data;
            density_data=density_gradient_data;}
        if(density_fraction_filename!="<unknown>"){
            LOG::cout<<"Using density fraction from "<<density_fraction_filename<<std::endl;
            ARRAY<T,VECTOR<int,3> > density_fraction_data;
            FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(density_fraction_filename,frame),density_fraction_data);
            *density_data*=density_fraction_data;}
        if(use_cubic_for_density){
            LOG::cout<<"Using density cubic"<<std::endl;
            for(typename GRID<TV>::CELL_ITERATOR iterator(*grid);iterator.Valid();iterator.Next()){VECTOR<int,TV::dimension> cell=iterator.Cell_Index();(*density_data)(cell)=(*density_data)(cell)*(*density_data)(cell)*(*density_data)(cell);}}
        LOG::cout<<"  Read grid " <<*grid;
        RENDERING_UNIFORM_VOXELS<T> *voxels;
        if(coarse_grid_filename!="<unknown>"){
            GRID<TV>* coarse_grid=new GRID<TV>;FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(coarse_grid_filename,frame),*coarse_grid);
            voxels=new RENDERING_UNIFORM_VOXELS<T>(*grid,*coarse_grid,*density_data,(T)volume_step);}
        else voxels=new RENDERING_UNIFORM_VOXELS<T>(*grid,*density_data,(T)volume_step);
        if(temperature_filename!="<unknown>"){
            ARRAY<T,VECTOR<int,3> >* temperature_data=new ARRAY<T,VECTOR<int,3> >;voxels->data.Append(temperature_data);
            FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(temperature_filename,frame),*temperature_data);}

        voxels->data_scale.Append(density_scale);
        voxels->data_scale.Append(temperature_scale);

        voxels->data_offset.Append(density_offset);
        voxels->data_offset.Append(temperature_offset);

        voxels->data_clamp_low_value.Append(clamp_low_density);
        voxels->data_clamp_low_value.Append(clamp_low_temperature);

        voxels->data_lowest_value.Append(density_lowest);
        voxels->data_lowest_value.Append(temperature_lowest);

        voxels->data_clamp_high_value.Append(clamp_high_density);
        voxels->data_clamp_high_value.Append(clamp_high_temperature);

        voxels->data_highest_value.Append(density_highest);
        voxels->data_highest_value.Append(temperature_highest);

        voxels->Print_Parameters();
        
        if(collision_aware){
            if(!body_list){LOG::cout<<"Error: Voxel_Data Use_Collision_Aware_Interpolation set to true, but collision list is null"<<std::endl;exit(1);}
            LOG::cout<<"Using collidable thin shell interpolation..."<<std::endl;
            GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >* fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >(*grid);
            for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);
            ARRAY<bool,VECTOR<int,3> >* cell_valid_mask=new ARRAY<bool,VECTOR<int,3> >(grid->Domain_Indices(3),false);cell_valid_mask->Fill(true);
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>* linear=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*fluid_collision_body_list,cell_valid_mask,0);
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,TV>* linear_vector=
                new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,TV>(*fluid_collision_body_list,cell_valid_mask,TV());
            GRID<TV> occupied_grid=grid->Get_MAC_Grid();
            fluid_collision_body_list->Rasterize_Objects();
            fluid_collision_body_list->Compute_Occupied_Blocks(false,(T)2*grid->Minimum_Edge_Length(),5);
            voxels->Set_Custom_Source_Interpolation(linear);
            voxels->Set_Custom_Light_Interpolation(linear_vector);}
        voxels->number_of_smoothing_steps=parameters.Get_Parameter("Number_Of_Smoothing_Steps",(int)3);
        voxels->precompute_single_scattering=precompute;
        object=voxels;}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    else if(type=="Octree_Voxel_Data"){
        LOG::cout<<"Octree Voxel Data"<<std::endl;
        bool precompute=parameters.Get_Parameter("Precompute_Single_Scattering",false);
        std::string octree_grid_filename=parameters.Get_Parameter("Octree_Grid_Filename",std::string("<unknown>"));
        std::string density_filename=parameters.Get_Parameter("Density_Filename",std::string("<unknown>"));
        std::string temperature_filename=parameters.Get_Parameter("Temperature_Filename",std::string("<unknown>"));
        bool collision_aware=parameters.Get_Parameter("Use_Collision_Aware_Interpolation",false);
        T volume_step=parameters.Get_Parameter("Volume_Step",(T)0.05);
        // always read densities but sometimes read temperatures
        OCTREE_GRID<T>* octree_grid=new OCTREE_GRID<T>;ARRAY<T>* density_data=new ARRAY<T>;
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(octree_grid_filename,frame),*octree_grid);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(density_filename,frame),*density_data);
        RENDERING_OCTREE_VOXELS<T> *voxels=new RENDERING_OCTREE_VOXELS<T>(*octree_grid,*density_data,volume_step);
        if(temperature_filename!="<unknown>"){
            ARRAY<T>* temperature_data=new ARRAY<T>;voxels->data.Append(temperature_data);
            FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(temperature_filename,frame),*temperature_data);}
        if(collision_aware){
            if(!body_list){LOG::cout<<"Error: Octree_Voxel_Data Use_Collision_Aware_Interpolation set to true, but collision list is null"<<std::endl;exit(1);}
            LOG::cout<<"Using collidable thin shell interpolation..."<<std::endl;
            GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >(*octree_grid);
            for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);

            ARRAY<bool>* cell_valid_mask=new ARRAY<bool>(octree_grid->number_of_cells,false);cell_valid_mask->Fill(false);
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,T>* linear=
                new LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,T>(*fluid_collision_body_list,cell_valid_mask,T());
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,TV>* linear_vector=
                new LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,TV>(*fluid_collision_body_list,cell_valid_mask,TV());
            fluid_collision_body_list->Rasterize_Objects();
            voxels->Set_Custom_Source_Interpolation(linear);
            voxels->Set_Custom_Light_Interpolation(linear_vector);}
        voxels->precompute_single_scattering=precompute;
        object=voxels;}
#endif
    else if(type=="Shock"){
        LOG::cout<<"Shock"<<std::endl;
        std::string grid_filename=parameters.Get_Parameter("Grid_Filename",std::string("<unknown>"));
        std::string density_filename=parameters.Get_Parameter("Density_Filename",std::string("<unknown>"));
        std::string pressure_filename=parameters.Get_Parameter("Pressure_Filename",std::string("<unknown>"));
        T gradient_threshold=parameters.Get_Parameter("Gradient_Threshold",(T)100);
        T refraction_multiplier=parameters.Get_Parameter("Refraction_Multiplier",(T)1);
        T volume_step=parameters.Get_Parameter("Volume_Step",(T)0.1);
        T fine_volumetric_step=parameters.Get_Parameter("Fine_Volumetric_Step",volume_step);
        T skip_next_intersection_factor=parameters.Get_Parameter("Skip_Next_Intersection_Factor",(T)10);
        bool use_pressure_for_intersection=parameters.Get_Parameter("Use_Pressure_For_Intersection",false);
        bool use_pressure_for_rarefaction=parameters.Get_Parameter("Use_Pressure_For_Rarefaction",false);

        GRID<TV>* grid=new GRID<TV>;
        ARRAY<T,VECTOR<int,3> >* density=new ARRAY<T,VECTOR<int,3> >;
        ARRAY<T,VECTOR<int,3> >* pressure=new ARRAY<T,VECTOR<int,3> >;
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(grid_filename,frame),*grid);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(density_filename,frame),*density);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(pressure_filename,frame),*pressure);
        LOG::cout<<"  Read grid " <<*grid <<std::endl;
        RENDERING_SHOCKS<T>* rendering_shocks=new RENDERING_SHOCKS<T>(*grid,*density,*pressure,gradient_threshold,refraction_multiplier,volume_step,fine_volumetric_step,skip_next_intersection_factor,use_pressure_for_intersection,use_pressure_for_rarefaction);
        object=rendering_shocks;}
    else if(type=="Levelset"){
        LOG::cout<<"Levelset"<<std::endl;
        ARRAY<T,VECTOR<int,3> >* phi=new ARRAY<T,VECTOR<int,3> >;GRID<TV>* grid=new GRID<TV>;
        RENDERING_IMPLICIT_SURFACE<T>* rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(*grid,*phi);
        LEVELSET_IMPLICIT_OBJECT<TV>* implicit_surface=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(rendering_implicit_surface->implicit_surface);
        std::string raw_filename=parameters.Get_Parameter("Filename",std::string("unknown"));
        bool negate=parameters.Get_Parameter("Negate",(bool)false);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(raw_filename,frame),*implicit_surface);
        LOG::cout<<"Grid "<<(*grid)<<std::endl;
        if(negate) for(int i=1;i<=grid->counts.x;i++)for(int j=1;j<=grid->counts.y;j++)for(int ij=1;ij<=grid->counts.z;ij++)(*phi)(i,j,ij)=-(*phi)(i,j,ij);
        T contour=parameters.Get_Parameter("Contour",(T)0);
        if(contour) *phi-=contour;
        int reinitialization_band=parameters.Get_Parameter("Reinitialization_Band",(int)5);
        std::string particle_processing_mode=STRING_UTILITIES::toupper(parameters.Get_Parameter("Particle_Processing_Mode",std::string("NONE")));

        bool collision_aware=parameters.Get_Parameter("Use_Collision_Aware_Interpolation",false);
        bool revalidate_for_subdivision=parameters.Get_Parameter("Revalidate_For_Subdivision",false);
        if(collision_aware){
            // remove ghost values (will get crash otherwise)
            if(phi->counts != grid->counts){
                LOG::cout<<"Removing phi ghost values"<<std::endl;
                phi->Resize(1,grid->counts.x,1,grid->counts.y,1,grid->counts.z);}
            LOG::cout<<"Using collidable thin shell interpolation..."<<std::endl;
            GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >* fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >(*grid);
            for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);

            if(revalidate_for_subdivision){ // TODO: clean this up
                ARRAY<bool,VECTOR<int,3> > current_arrays_valid(grid->Domain_Indices(3)),next_arrays_valid(grid->Domain_Indices(3));
                for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) (*body_list)(i)->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
                fluid_collision_body_list->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,(T)0);
                for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) (*body_list)(i)->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
                fluid_collision_body_list->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,(T)1);
                ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T> advection(*fluid_collision_body_list,current_arrays_valid,next_arrays_valid,(T)1e-5,true);

                fluid_collision_body_list->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);

                fluid_collision_body_list->Initialize_Grids();
                // get swept structures
                fluid_collision_body_list->Compute_Occupied_Blocks(false,(T)2*(*grid).DX().Max(),10);
                fluid_collision_body_list->Compute_Occupied_Blocks(true,(T).5*(*grid).DX().Max(),10);

                fluid_collision_body_list->Rasterize_Objects();
                fluid_collision_body_list->Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
                    COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);

                const T dt=(T)1;
                int count=0;
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(*grid);iterator.Valid();iterator.Next()){
                    VECTOR<int,3> cell=iterator.Cell_Index();
                    if(fluid_collision_body_list->Swept_Occupied_Cell_Center(cell) && fluid_collision_body_list->Latest_Cell_Crossover(cell,dt)){
                        count++;
                        current_arrays_valid(cell)=false;}
                    else current_arrays_valid(cell)=true;}

                LOG::cout<<count<<" invalidated cells"<<std::endl;
                FILE_UTILITIES::Write_To_File<T>("phi_before.phi",implicit_surface->levelset);
                fluid_collision_body_list->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
                fluid_collision_body_list->Compute_Grid_Visibility();
                advection.Average_To_Invalidated_Cells(*grid,(T)1e-5,implicit_surface->levelset.phi);
                FILE_UTILITIES::Write_To_File<T>("phi_after.phi",implicit_surface->levelset);
                fluid_collision_body_list->Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
                    COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);}

            ARRAY<bool,VECTOR<int,3> >* cell_valid_mask=new ARRAY<bool,VECTOR<int,3> >(grid->Domain_Indices(3),false);cell_valid_mask->Fill(true);
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV> ,T>* linear=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV> ,T>(*fluid_collision_body_list,cell_valid_mask,
                implicit_surface->levelset.collidable_phi_replacement_value);
            GRID<TV> occupied_grid=grid->Get_MAC_Grid();
            fluid_collision_body_list->Rasterize_Objects();
            fluid_collision_body_list->Compute_Occupied_Blocks(false,(T)2*grid->Minimum_Edge_Length(),5);
            implicit_surface->Set_Custom_Secondary_Interpolation(*linear);
            implicit_surface->Use_Secondary_Interpolation(true);}

        bool reinitialize=parameters.Get_Parameter("Reinitialize",false);
        if(reinitialize){
            LOG::cout<<"Reinitializing levelset"<<std::endl;
            implicit_surface->levelset.Fast_Marching_Method(0,reinitialization_band*grid->min_dX);}

        object=rendering_implicit_surface;}
    else if(type=="Levelset_Multiple"){
        LOG::cout<<"Multiple Levelset"<<std::endl;
        int number_of_regions=parameters.Get_Parameter("Number_Regions", 0);
        GRID<TV>* grid=new GRID<TV>;ARRAY<ARRAY<T,VECTOR<int,3> > >* phis=new ARRAY<ARRAY<T,VECTOR<int,3> > >(number_of_regions);
        RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<TV> > >* rendering_levelset_multiple_object=new RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<TV> > >(*grid,*phis);
        // read in the data
        for(int i=1;i<=number_of_regions;i++){
            LOG::cout<<"Reading Region("<<i<<")"<<std::endl;
            std::string raw_filename=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Filename%d",i), std::string("unknown"));
            FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(raw_filename,frame),*rendering_levelset_multiple_object->levelset_multiple.levelsets(i));
            *grid=rendering_levelset_multiple_object->levelset_multiple.levelsets(i)->grid;}
        for(int i=1;i<=number_of_regions;i++){
            // setup shaders before the rendering_multiple_implicit_surface is casted to rendering_object;
            std::string shader_name=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Shader%d",i),std::string("<unknown>"));
            std::string volume_shader_name=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Volume_Shader%d",i),std::string("<unknown>"));
            std::string name=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Name%d",i),std::string("<unknown>"));
            // check to make sure at least one shader specified
            if(shader_name=="<unknown>"&&volume_shader_name=="<unknown>") PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No shader (volume or surface) specified for object %s\n",name.c_str()));
            RENDERING_OBJECT<T> *object_local=rendering_levelset_multiple_object->rendering_levelset_multiple_region_objects(i);
            object_local->index_of_refraction=parameters.Get_Parameter(STRING_UTILITIES::string_sprintf("Index_Of_Refraction%d",i),(T)1);
            object_local->name=name;
            object=rendering_levelset_multiple_object;object->Update_Transform(current_transform);object->priority=priority;
            object->support_transparent_overlapping_objects=support_transparent_overlapping_objects;
            // volume shader (volume shader should go first not to add object_local to standard object list of world
            if(volume_shader_name!="<unknown>"){
                if(volume_shaders.Get(volume_shader_name,object_local->volumetric_shader)){
                    object_local->add_to_spatial_partition=false;
                    world.Add_Object(object_local);}
                else PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid volume shader '%s' specified for object %s\n",volume_shader_name.c_str(),name.c_str()));}
            // surface shader
            if(shader_name!="<unknown>"){
                if(shaders.Get(shader_name,object->material_shader)) object_local->material_shader=object->material_shader;
                else PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid shader '%s' specified for object %s\n",shader_name.c_str(),name.c_str()));}

            object_local->name=name;object_local->Update_Transform(current_transform);object_local->priority=priority;
            object_local->support_transparent_overlapping_objects=support_transparent_overlapping_objects;
        }
        object->name=parameters.Get_Parameter("Name",std::string("<unknown>"));
        objects.Set(object->name,object);
    }
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    else if(type=="Octree_Levelset"){
        LOG::cout<<"Octree_Levelset"<<std::endl;
        ARRAY<T>* phi=new ARRAY<T>;OCTREE_GRID<T>* grid=new OCTREE_GRID<T>;
        std::string raw_grid_filename=parameters.Get_Parameter("Grid_Filename",std::string("unknown"));
        std::string raw_phi_filename=parameters.Get_Parameter("Phi_Filename",std::string("unknown"));
        bool negate=parameters.Get_Parameter("Negate",(bool)false);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(raw_grid_filename,frame),*grid);
        FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(raw_phi_filename,frame),*phi);
        if(negate){for(int i=1;i<=grid->number_of_nodes;i++)(*phi)(i)=-(*phi)(i);}

        std::string particle_processing_mode=STRING_UTILITIES::toupper(parameters.Get_Parameter("Particle_Processing_Mode",std::string("NONE")));
        if(particle_processing_mode!="NONE"){
            std::string raw_removed_particles_filename=parameters.Get_Parameter("Removed_Negative_Particles_Filename",std::string("unknown"));
            if(raw_removed_particles_filename!="unknown"){ // merge with removed negative particles
                LOG::cout<<"Reading removed negative particles"<<std::endl;
                ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*> particle_array;
                FILE_UTILITIES::template Read_From_File<RW>(Animated_Filename(raw_removed_particles_filename,frame),particle_array);
                OCTREE_REMOVED_PARTICLES_PROCESSING<T> particle_processing(*grid,*phi,particle_array);
                particle_processing.blending_parameter=parameters.Get_Parameter("Particle_Blending_Parameter",(T).8);
                particle_processing.scale=parameters.Get_Parameter("Particle_Scale",(T).25);
                particle_processing.octree_maximum_depth=parameters.Get_Parameter("Particle_Octree_Maximum_Depth",(int)2);
                particle_processing.particle_power=parameters.Get_Parameter("Particle_Power",(T)1);
                particle_processing.use_velocity_scaling=parameters.Get_Parameter("Particle_Use_Velocity_Scaling",(bool)true);
                particle_processing.dt=parameters.Get_Parameter("Particle_Dt",(T)1/24);
                particle_processing.preserve_volume=parameters.Get_Parameter("Particle_Preserve_Volume",(bool)true);
                bool collision_aware_particle_processing=parameters.Get_Parameter("Collision_Aware_Particle_Processing",(bool)true);
                if(collision_aware_particle_processing && body_list){
                    LOG::cout<<"Particle processing will be collision aware"<<std::endl;
                    GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >(*grid);
                    for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);
                    particle_processing.Set_Collision_Aware(fluid_collision_body_list);}
                particle_processing.Refine_And_Create_Particle_Phi();
                particle_array.Delete_Pointers_And_Clean_Memory();

                ARRAY<T>* combined_phi=new ARRAY<T>;
                if(particle_processing_mode=="MERGE"){
                    LOG::cout<<"Merging particles with level set"<<std::endl;
                    particle_processing.Merge_Phi(*combined_phi);}
                else if(particle_processing_mode=="UNION"){
                    LOG::cout<<"Taking union of particles with level set"<<std::endl;
                    particle_processing.Union_Phi(*combined_phi);}
                else if(particle_processing_mode=="BLEND"){
                    T blend_cells=parameters.Get_Parameter("Particle_Blend_Cells",(T)5);
                    LOG::cout<<"Taking blend of particles with level set (blend_cells="<<blend_cells<<")"<<std::endl;
                    particle_processing.Blend_Phi(*combined_phi,blend_cells);}
                else{
                    LOG::cout<<"Unrecognized particle processing mode: "<<particle_processing_mode<<std::endl;
                    exit(1);}

                // switch over to new phi
                delete phi;phi=combined_phi;

                std::string raw_combined_grid_output_filename=parameters.Get_Parameter("Combined_Grid_Output_Filename",std::string("unknown"));
                if(raw_combined_grid_output_filename!="unknown"){
                    std::string filename=Animated_Filename(raw_combined_grid_output_filename,frame);
                    LOG::cout<<"Writing combined grid to "<<filename<<std::endl;
                    FILE_UTILITIES::template Write_To_File<RW>(filename,*grid);}
                std::string raw_combined_phi_output_filename=parameters.Get_Parameter("Combined_Phi_Output_Filename",std::string("unknown"));
                if(raw_combined_phi_output_filename!="unknown"){
                    std::string filename=Animated_Filename(raw_combined_phi_output_filename,frame);
                    LOG::cout<<"Writing combined phi to "<<filename<<std::endl;
                    FILE_UTILITIES::template Write_To_File<RW>(filename,*phi);}
            }
        }
        RENDERING_IMPLICIT_SURFACE<T>* rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(*grid,*phi);
        DYADIC_IMPLICIT_OBJECT<TV>& implicit_surface=dynamic_cast<DYADIC_IMPLICIT_OBJECT<TV>&>(*rendering_implicit_surface->implicit_surface);
        bool collision_aware=parameters.Get_Parameter("Use_Collision_Aware_Interpolation",false);
        bool reinitialize=parameters.Get_Parameter("Reinitialize",false);
        int reinitialization_band=parameters.Get_Parameter("Reinitialization_Band",(int)5);
        if(reinitialize){
            LOG::cout<<"Reinitializing levelset"<<std::endl;
            implicit_surface.levelset.Fast_Marching_Method(0,reinitialization_band*grid->Minimum_Edge_Length());}
        if(collision_aware){
            LOG::cout<<"Using collidable thin shell interpolation..."<<std::endl;
            GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >(*grid);
            for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);
            ARRAY<bool>* cell_valid_mask=new ARRAY<bool>(grid->number_of_cells,false);cell_valid_mask->Fill(true);
            LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,T>* linear=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<OCTREE_GRID<T>,T>(*fluid_collision_body_list,
                cell_valid_mask,implicit_surface.levelset.collidable_phi_replacement_value);
            fluid_collision_body_list->Rasterize_Objects();
            implicit_surface.Set_Custom_Secondary_Interpolation(*linear);
            implicit_surface.Use_Secondary_Interpolation(true);
            assert(rendering_implicit_surface->implicit_surface->use_secondary_interpolation);}

        object=rendering_implicit_surface;}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    else if(type=="RLE_Levelset"){
        LOG::cout<<"RLE_Levelset"<<std::endl;
        std::string raw_grid_filename=parameters.Get_Parameter("Grid_Filename",std::string("<unknown_grid_file>"));
        std::string raw_phi_filename=parameters.Get_Parameter("Phi_Filename",std::string("<unknown_phi_file>"));
        RLE_GRID_3D<T>* grid=new RLE_GRID_3D<T>;ARRAY<T>* phi=new ARRAY<T>;
        FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_grid_filename,frame),*grid);
        FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_phi_filename,frame),*phi);
        LOG::cout<<"uniform grid "<<grid->uniform_grid<<", number_of_cells = "<<grid->number_of_cells<<std::endl;
        RENDERING_IMPLICIT_SURFACE<T>* rendering_implicit_surface;
        LEVELSET_RLE<RLE_GRID_3D<T> >* levelset;
        bool use_fine_levelset=parameters.Get_Parameter("Use_Fine_Levelset",false);
        if(!use_fine_levelset){
            rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(*grid,*phi);
            levelset=&((RLE_IMPLICIT_OBJECT<TV>*)rendering_implicit_surface->implicit_surface)->levelset;}
        else{
            std::string raw_fine_grid_filename=parameters.Get_Parameter("Fine_Grid_Filename",std::string("<unknown_fine_grid_file>"));
            std::string raw_fine_phi_filename=parameters.Get_Parameter("Fine_Phi_Filename",std::string("<unknown_fine_phi_file>"));
            RLE_GRID_3D<T>* fine_grid=new RLE_GRID_3D<T>;ARRAY<T>* fine_phi=new ARRAY<T>;
            FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_fine_grid_filename,frame),*fine_grid);
            FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_fine_phi_filename,frame),*fine_phi);
            LOG::cout<<"fine uniform grid "<<fine_grid->uniform_grid<<", number_of_cells = "<<fine_grid->number_of_cells<<std::endl;
            RLE_IMPLICIT_SURFACE_TWO_LEVEL<T>* rle_implicit_surface_two_level=new RLE_IMPLICIT_SURFACE_TWO_LEVEL<T>(*grid,*phi,*fine_grid,*fine_phi);
            rendering_implicit_surface=new RENDERING_IMPLICIT_SURFACE<T>(rle_implicit_surface_two_level);
            levelset=&rle_implicit_surface_two_level->levelset;}
        bool negate=parameters.Get_Parameter("Negate",(bool)false);
        if(negate) *phi*=-1;
        T contour=parameters.Get_Parameter("Contour",(T)0);
        if(contour) *phi-=contour;
        bool reinitialize=parameters.Get_Parameter("Reinitialize",false);
        if(reinitialize){
            LOG::cout<<"Reinitializing levelset"<<std::endl;
            levelset->Fast_Marching_Method(0);}
        bool extrapolate_into_objects=parameters.Get_Parameter("Extrapolate_Into_Objects",false);
        bool use_collision_body_list_for_extrapolation=parameters.Get_Parameter("Use_Collision_Body_List_For_Extrapolation",false);
        bool subtract_objects=parameters.Get_Parameter("Subtract_Objects",false);
        bool collision_aware=parameters.Get_Parameter("Use_Collision_Aware_Interpolation",false);
        GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_3D<T> >* fluid_collision_body_list=0;
        if(use_collision_body_list_for_extrapolation || collision_aware){
            fluid_collision_body_list=new GRID_BASED_COLLISION_GEOMETRY_RLE<RLE_GRID_3D<T> >(*grid);
            if(!body_list) PHYSBAM_FATAL_ERROR();
            for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++) if((*body_list)(i)) fluid_collision_body_list->collision_geometry_collection.Add_Body((*body_list)(i),0,false);}
        if(extrapolate_into_objects || subtract_objects){
            ARRAY<T> phi_object;
            if(!use_collision_body_list_for_extrapolation){
                std::string raw_phi_object_filename=parameters.Get_Parameter("Phi_Object_Filename",std::string(""));
                FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_phi_object_filename,frame),phi_object);}
            else{
                phi_object.Resize(grid->number_of_cells,false,false);phi_object.Fill(1);
                if(!body_list) PHYSBAM_FATAL_ERROR();
                T contour=3*grid->Minimum_Edge_Length();
                for(COLLISION_GEOMETRY_ID i(1);i<=body_list->m;i++)
                    for(RLE_GRID_ITERATOR_CELL_3D<T> cell(*grid,grid->number_of_ghost_cells);cell;cell++){int c=cell.Cell();TV X=cell.X();
                        T phi;if(fluid_collision_body_list->collision_geometry_collection.bodies(i)->Implicit_Geometry_Lazy_Inside_And_Value(X,phi,contour)) phi_object(c)=min(phi_object(c),phi);}}
            phi_object*=-1; // extrapolate into negative phi_object
            if(extrapolate_into_objects){
                LOG::cout<<"extrapolating levelset into object"<<std::endl;
                EXTRAPOLATION_RLE<RLE_GRID_3D<T>,T> extrapolate(*grid,phi_object);extrapolate.Extrapolate_Cells(*phi,0);}
            if(subtract_objects) for(int c=1;c<=grid->number_of_cells;c++)(*phi)(c)=max((*phi)(c),phi_object(c));}
        if(collision_aware){
            LOG::cout<<"Using collidable thin shell interpolation..."<<std::endl;
            levelset->valid_mask_current.Resize(grid->number_of_cells,false);levelset->valid_mask_current.Fill(true);
            levelset->Set_Collision_Body_List(*fluid_collision_body_list,true);
            fluid_collision_body_list->Rasterize_Objects();
            rendering_implicit_surface->implicit_surface->Use_Secondary_Interpolation(true);}
        object=rendering_implicit_surface;}
#endif
    else if(type=="Surface_Of_Revolution"){
        LOG::cout<<"Surface_Of_Revolution"<<std::endl;
        std::string filename=parameters.Get_Parameter("Filename",std::string("unknown"));
        T width=parameters.Get_Parameter("Width",(T).0412);
        T height=parameters.Get_Parameter("Height",(T).1612);
        T tolerance=parameters.Get_Parameter("Tolerance",(T)1e-6);
        SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T>* surface_of_revolution=new SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T>(filename,width,height,tolerance);
        object=new RENDERING_IMPLICIT_SURFACE<T>(surface_of_revolution);
        object->priority=priority;object->support_transparent_overlapping_objects=support_transparent_overlapping_objects;}
    else if(type=="Particle"){
        LOG::cout<<"Particle"<<std::endl;
        T scale(parameters.Get_Parameter("Scale",T(1.0)));
        std::string particle_levelset(parameters.Get_Parameter("Levelset",std::string("unknown")));
        if(particle_levelset=="unknown"){LOG::cerr<<"Error: must specify levelset for particle"<<std::endl;exit(1);}
        ARRAY<POINT_CLOUD<TV>*,VECTOR<int,3> > particles_array;
        std::string raw_filename(parameters.Get_Parameter("Filename",std::string("unknown")));
        FILE_UTILITIES::Read_From_File<RW>(Animated_Filename(raw_filename,frame),particles_array);
        RENDERING_OBJECT<T>* particle_levelset_object;
        if(!objects.Get(particle_levelset,particle_levelset_object)){LOG::cerr<<"Error: levelset "<<particle_levelset<<" not defined"<<std::endl;exit(1);}
        RENDERING_IMPLICIT_SURFACE<T>* rendering_levelset=dynamic_cast<RENDERING_IMPLICIT_SURFACE<T>*>(particle_levelset_object);
        if(rendering_levelset==NULL){LOG::cerr<<"Error: particle "<<name<<", "<<particle_levelset<<" is not a valid levelset"<<std::endl;exit(1);}
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_surface=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(rendering_levelset->implicit_surface);
        if(levelset_implicit_surface==NULL){LOG::cerr<<"Error: particle "<<name<<", "<<particle_levelset<<" is not a valid levelset"<<std::endl;exit(1);}
        RENDERING_PARTICLES<T>* rendering_particle=new RENDERING_PARTICLES<T>(particles_array,levelset_implicit_surface->levelset.grid,scale);
        object=rendering_particle; object->add_to_spatial_partition=true;}
    else PHYSBAM_FATAL_ERROR("Unknown object type "+type);
    // Register the object if it was created
    if(object && type!="Levelset_Multiple"){    // shaders of levelset-multiple are already assigned to region objects.
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        std::string volume_shader_name=parameters.Get_Parameter("Volume_Shader",std::string("<unknown>"));
        object->index_of_refraction=parameters.Get_Parameter("Index_Of_Refraction",(T)1);
        // check to make sure at least one shader specified
        if(shader_name=="<unknown>"&&volume_shader_name=="<unknown>") PHYSBAM_FATAL_ERROR("No shader (volume or surface) specified for object "+name);
        // surface shader
        if(shader_name!="<unknown>"){
            if(!shaders.Get(shader_name,object->material_shader)) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid shader '%s' specified for object %s\n",shader_name.c_str(),name.c_str()));}
        // volume shader
        if(volume_shader_name!="<unknown>"){
            if(!volume_shaders.Get(volume_shader_name,object->volumetric_shader))
                PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid volume shader '%s' specified for object %s\n",volume_shader_name.c_str(),name.c_str()));
            object->add_to_spatial_partition=false;}
        object->name=name;object->Update_Transform(current_transform);objects.Set(name,object);
        object->priority=priority;object->support_transparent_overlapping_objects=support_transparent_overlapping_objects;}
    else if(object_group.m){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        std::string volume_shader_name=parameters.Get_Parameter("Volume_Shader",std::string("<unknown>"));
        T index_of_refraction=parameters.Get_Parameter("Index_Of_Refraction",(T)1);
        for(int i=1;i<=object_group.m;i++){
            object_group(i)->index_of_refraction=index_of_refraction;
            object_group(i)->priority=priority;object_group(i)->support_transparent_overlapping_objects=support_transparent_overlapping_objects;}
        // check to make sure at least one shader specified
        if(shader_name=="<unknown>"&&volume_shader_name=="<unknown>") PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No shader (volume or surface) specified for object %s\n",name.c_str()));
        // surface shader
        if(shader_name!="<unknown>"){
            MATERIAL_SHADER<T>* material_shader=0;
            if(!shaders.Get(shader_name,material_shader)) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid shader '%s' specified for object %s\n",shader_name.c_str(),name.c_str()));
            for(int i=1;i<=object_group.m;i++) object_group(i)->material_shader=material_shader;}
        // volume shader
        if(volume_shader_name!="<unknown>"){
            VOLUMETRIC_SHADER<T>* volumetric_shader=0;
            if(!volume_shaders.Get(volume_shader_name,volumetric_shader)) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid volume shader '%s' specified for object %s\n",volume_shader_name.c_str(),name.c_str()));
            for(int i=1;i<=object_group.m;i++){object_group(i)->volumetric_shader=volumetric_shader;object_group(i)->add_to_spatial_partition=false;}}
        for(int i=1;i<=object_group.m;i++){
            std::string instance_name=STRING_UTILITIES::string_sprintf("%s_%d",name.c_str(),i);
            object_group(i)->name=instance_name;object_group(i)->Update_Transform(current_transform);objects.Set(instance_name,object_group(i));}}
    if(object) GENERIC_RENDER_EXAMPLE<T,RW>::Add_Solid_Texture(object,parameters);
}
//#####################################################################
template void GENERIC_RENDER_EXAMPLE<float,float>::Object(RENDER_WORLD<float>& world,const int frame,PARAMETER_LIST& parameters);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void GENERIC_RENDER_EXAMPLE<double,float>::Object(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,double>::Object(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
#endif
