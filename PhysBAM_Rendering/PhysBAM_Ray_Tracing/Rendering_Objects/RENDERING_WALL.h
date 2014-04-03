//#####################################################################
// Copyright 2005, Andrew Selle, Fen Zhao.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_3D
//#####################################################################
#ifndef __RENDERING_WALL__
#define __RENDERING_WALL__

#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class RENDERING_WALL:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    PLANE<T> xmin,xmax,ymin,ymax,zmin,zmax;
    bool show_xmin,show_xmax,show_ymin,show_ymax,show_zmin,show_zmax;
    RANGE<TV> box;
    TV texture_vector_x_1,texture_vector_x_2,texture_vector_y_1,texture_vector_y_2,texture_vector_z_1,texture_vector_z_2;//Texture vectors

    RENDERING_WALL()
    {show_xmax=show_xmin=show_ymax=show_ymin=show_zmax=show_zmin=true;}
    
    RENDERING_WALL(const RANGE<TV>& box_input,const bool show_xmin_input,const bool show_xmax_input,const bool show_ymin_input,const bool show_ymax_input,const bool show_zmin_input,const bool show_zmax_input,const T texture_scale)
        :xmin(TV(-1,0,0),TV(box_input.min_corner.x,box_input.min_corner.y,box_input.min_corner.z)),xmax(TV(1,0,0),TV(box_input.max_corner.x,box_input.min_corner.y,box_input.min_corner.z)),
         ymin(TV(0,-1,0),TV(box_input.min_corner.x,box_input.min_corner.y,box_input.min_corner.z)),ymax(TV(0,1,0),TV(box_input.min_corner.x,box_input.max_corner.y,box_input.min_corner.z)),
         zmin(TV(0,0,-1),TV(box_input.min_corner.x,box_input.min_corner.y,box_input.min_corner.z)),zmax(TV(0,0,1),TV(box_input.min_corner.x,box_input.min_corner.y,box_input.max_corner.z)),box(box_input),
         texture_vector_x_1(0,0,texture_scale),texture_vector_x_2(0,texture_scale,0),texture_vector_y_1(texture_scale,0,0),texture_vector_y_2(0,0,texture_scale),texture_vector_z_1(texture_scale,0,0),texture_vector_z_2(0,texture_scale,0)
    {show_xmax=show_xmax_input;show_xmin=show_xmin_input;show_ymax=show_ymax_input;show_ymin=show_ymin_input;show_zmax=show_zmax_input;show_zmin=show_zmin_input;}

    virtual ~RENDERING_WALL()
    {}
    
    bool Intersection(RAY<TV>& ray) const  PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(Intersection_With_Shown_Planes(object_space_ray)){ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    bool Intersection_With_Shown_Planes(RAY<TV>& ray) const
    {T distance;bool intersected=false;
    T xmin_minus=box.min_corner.x-small_number,xmax_plus=box.max_corner.x+small_number,ymin_minus=box.min_corner.y-small_number,ymax_plus=box.max_corner.y+small_number,zmin_minus=box.min_corner.z-small_number,zmax_plus=box.max_corner.z+small_number;
        
    if(show_xmin || show_xmax){
        T rate_of_approach=-ray.direction.x;
        if(show_xmin){
            distance=ray.endpoint.x-box.min_corner.x;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,xmin,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.y >= ymin_minus && point.y <= ymax_plus && point.z >= zmin_minus && point.z <= zmax_plus){intersected=true;temp_ray.aggregate_id=1;temp_ray.Save_Intersection_Information(ray);}}}
        if(show_xmax){
            distance=ray.endpoint.x-box.max_corner.x;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,xmax,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.y >= ymin_minus && point.y <= ymax_plus && point.z >= zmin_minus && point.z <= zmax_plus){intersected=true;temp_ray.aggregate_id=2;temp_ray.Save_Intersection_Information(ray);}}}}

    if(show_ymin || show_ymax){
        T rate_of_approach=-ray.direction.y;
        if(show_ymin){
            distance=ray.endpoint.y-box.min_corner.y;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,ymin,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.x >= xmin_minus && point.x <= xmax_plus && point.z >= zmin_minus && point.z <= zmax_plus){intersected=true;temp_ray.aggregate_id=3;temp_ray.Save_Intersection_Information(ray);}}}
        if(show_ymax){
            distance=ray.endpoint.y-box.max_corner.y;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,ymax,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.x >= xmin_minus && point.x <= xmax_plus && point.z >= zmin_minus && point.z <= zmax_plus){intersected=true;temp_ray.aggregate_id=4;temp_ray.Save_Intersection_Information(ray);}}}}

    if(show_zmin || show_zmax){
        T rate_of_approach=-ray.direction.z;
        if(show_zmin){
            distance=ray.endpoint.z-box.min_corner.z;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,zmin,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.x >= xmin_minus && point.x <= xmax_plus &&  point.y >= ymin_minus && point.y <= ymax_plus){intersected=true;temp_ray.aggregate_id=5;temp_ray.Save_Intersection_Information(ray);}}}
        if(show_zmax){
            distance=ray.endpoint.z-box.max_corner.z;RAY<TV> temp_ray=ray;
            if(INTERSECTION::Intersects(temp_ray,zmax,small_number,distance,rate_of_approach)){
                TV point=temp_ray.Point(temp_ray.t_max);
                if(point.x >= xmin_minus && point.x <= xmax_plus &&  point.y >= ymin_minus && point.y <= ymax_plus){intersected=true;temp_ray.aggregate_id=6;temp_ray.Save_Intersection_Information(ray);}}}}
    return intersected;}

    TV Normal(const TV& location,const int aggregate) const PHYSBAM_OVERRIDE
    {assert(aggregate >= 1 && aggregate <= 6);
    if(aggregate == 1) return World_Space_Vector(xmin.Normal());
    else if(aggregate == 2) return World_Space_Vector(xmax.Normal());
    else if(aggregate == 3) return World_Space_Vector(ymin.Normal());
    else if(aggregate == 4) return World_Space_Vector(ymax.Normal());
    else if(aggregate == 5) return World_Space_Vector(zmin.Normal());
    else return World_Space_Vector(zmax.Normal());}
    
    bool Inside(const TV& location) const  PHYSBAM_OVERRIDE
    {return box.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return box.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return box.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& world_location) const PHYSBAM_OVERRIDE
    {TV location=Object_Space_Point(world_location);
    if(box.Lazy_Inside(location)){
        int side=0;T distance=FLT_MAX;
        if(show_xmin && location.x-box.min_corner.x < distance){side=1;distance=location.x-box.min_corner.x;}
        if(show_xmax && box.max_corner.x-location.x < distance){side=2;distance=box.max_corner.x-location.x;}
        if(show_ymin && location.y-box.min_corner.y < distance){side=3;distance=location.y-box.min_corner.y;}
        if(show_ymax && box.max_corner.y-location.y < distance){side=4;distance=box.max_corner.y-location.y;}
        if(show_zmin && location.z-box.min_corner.z < distance){side=5;distance=location.z-box.min_corner.z;}
        if(show_zmax && box.max_corner.z-location.z < distance){side=6;distance=box.max_corner.z-location.z;}
        switch(side){
          case 1:return TV(box.min_corner.x,location.y,location.z);
          case 2:return TV(box.max_corner.x,location.y,location.z);
          case 3:return TV(location.x,box.min_corner.y,location.z);
          case 4:return TV(location.x,box.max_corner.y,location.z);
          case 5:return TV(location.x,location.y,box.min_corner.z);
          case 6:return TV(location.x,location.y,box.max_corner.z);
          default: return TV(clamp(location.x,box.min_corner.x,box.max_corner.x),clamp(location.y,box.min_corner.y,box.max_corner.y),clamp(location.z,box.min_corner.z,box.max_corner.z));}}
    else return TV(clamp(location.x,box.min_corner.x,box.max_corner.x),clamp(location.y,box.min_corner.y,box.max_corner.y),clamp(location.z,box.min_corner.z,box.max_corner.z));}
    
    TRIANGULATED_SURFACE<T>* Generate_Triangles()const PHYSBAM_OVERRIDE
    {TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
   
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    int vertex_1=particles.array_collection->Add_Element(),vertex_2=particles.array_collection->Add_Element(),vertex_3=particles.array_collection->Add_Element(),vertex_4=particles.array_collection->Add_Element();
    int vertex_5=particles.array_collection->Add_Element(),vertex_6=particles.array_collection->Add_Element(),vertex_7=particles.array_collection->Add_Element(),vertex_8=particles.array_collection->Add_Element();
    particles.X(vertex_1)=(TV(box.min_corner.x,box.min_corner.y,box.max_corner.z));particles.X(vertex_2)=(TV(box.max_corner.x,box.min_corner.y,box.max_corner.z));particles.X(vertex_3)=(TV(box.max_corner.x,box.max_corner.y,box.max_corner.z));particles.X(vertex_4)=(TV(box.min_corner.x,box.max_corner.y,box.max_corner.z));
    particles.X(vertex_5)=(TV(box.min_corner.x,box.min_corner.y,box.min_corner.z));particles.X(vertex_6)=(TV(box.max_corner.x,box.min_corner.y,box.min_corner.z));particles.X(vertex_7)=(TV(box.max_corner.x,box.max_corner.y,box.min_corner.z));particles.X(vertex_8)=(TV(box.min_corner.x,box.max_corner.y,box.min_corner.z));
    
    int face_count=0;
    if(show_xmin){face_count++;};if(show_xmax){face_count++;};if(show_ymin){face_count++;};if(show_ymax){face_count++;};if(show_zmin){face_count++;};if(show_zmax){face_count++;};
    int triangle_count=1;ARRAY<VECTOR<int,3> > triangles(2*face_count);
    if(show_xmin){triangles(triangle_count).Set(1,4,8);triangle_count++;triangles(triangle_count).Set(8,5,1);triangle_count++;};
    if(show_xmax){triangles(triangle_count).Set(2,6,7);triangle_count++;triangles(triangle_count).Set(7,3,2);triangle_count++;};
    if(show_ymin){triangles(triangle_count).Set(5,1,2);triangle_count++;triangles(triangle_count).Set(2,6,5);triangle_count++;};
    if(show_ymax){triangles(triangle_count).Set(8,7,3);triangle_count++;triangles(triangle_count).Set(3,4,8);triangle_count++;};
    if(show_zmin){triangles(triangle_count).Set(5,8,7);triangle_count++;triangles(triangle_count).Set(7,6,5);triangle_count++;};
    if(show_zmax){triangles(triangle_count).Set(1,2,3);triangle_count++;triangles(triangle_count).Set(3,4,1);triangle_count++;};

    surface->mesh.Initialize_Mesh(particles.array_collection->Size(),triangles);
    surface->Update_Triangle_List();surface->Update_Vertex_Normals();return surface;}
    
    virtual void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const
    {assert(aggregate >= 1 && aggregate <= 6);
    TV x1;TV texture_vector1,texture_vector2;
    if(aggregate==1){x1=xmin.x1; texture_vector1=texture_vector_x_1;texture_vector2=texture_vector_x_2;}
    else if(aggregate==2){x1=xmax.x1;texture_vector1=texture_vector_x_1;texture_vector2=texture_vector_x_2;}
    else if(aggregate==3){x1=ymin.x1;texture_vector1=texture_vector_y_1;texture_vector2=texture_vector_y_2;}
    else if(aggregate==4){x1=ymax.x1;texture_vector1=texture_vector_y_1;texture_vector2=texture_vector_y_2;}
    else if(aggregate==5){x1=zmin.x1;texture_vector1=texture_vector_z_1;texture_vector2=texture_vector_z_2;}
    else if(aggregate==6){x1=zmax.x1;texture_vector1=texture_vector_z_1;texture_vector2=texture_vector_z_2;}
    TV p=object_space_point-x1;s=TV::Dot_Product(p,texture_vector1);t=TV::Dot_Product(p,texture_vector2);}
    

//#####################################################################
};
}
#endif

