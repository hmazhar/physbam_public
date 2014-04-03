//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_PARTICLES
//##################################################################### 
#ifndef __RENDERING_PARTICLES__
#define __RENDERING_PARTICLES__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM {

template<class T>
class RENDERING_PARTICLES:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::transform;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
private:
    ARRAY<POINT_CLOUD<TV>*,VECTOR<int,3> > particles_array;
    const GRID<TV>& grid;
    T scale;
    ARRAY<ARRAY<int> ,VECTOR<int,3> > particle_to_aggregate_id;
    ARRAY<PAIR<VECTOR<int,3>,int> > aggregate_id_to_particle;
public:

    RENDERING_PARTICLES(const ARRAY<POINT_CLOUD<TV>*,VECTOR<int,3> >& p,const GRID<TV>& g,T s)
        :particles_array(p),grid(g),scale(s),
         particle_to_aggregate_id(particles_array.domain.min_corner.x,particles_array.domain.max_corner.x,particles_array.domain.min_corner.y,particles_array.domain.max_corner.y,particles_array.domain.min_corner.z,particles_array.domain.max_corner.z,true),
         aggregate_id_to_particle()
    {Create_Aggregate_Ids();}


    virtual ~RENDERING_PARTICLES() 
    {}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    bool intersection=false;
    for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++) for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++) for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++) {POINT_CLOUD<TV>* particles=particles_array(i,j,ij);
        if(particles){
            ARRAY_VIEW<T> radius=*particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
            for(int p=1;p<=particles->array_collection->Size();++p) if(INTERSECTION::Intersects(object_space_ray,SPHERE<TV>(particles->X(p),scale*radius(p)),small_number)) {
                object_space_ray.aggregate_id=particle_to_aggregate_id(i,j,ij)(p);intersection=true;}}}
    if(intersection){ray.t_max=object_space_ray.t_max;ray.semi_infinite=object_space_ray.semi_infinite;ray.aggregate_id=object_space_ray.aggregate_id;}
    return intersection;}

    bool Intersection(RAY<TV>& ray,const int aggregate) const  PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);bool intersection=false;
    POINT_CLOUD<TV>* particles=particles_array(aggregate_id_to_particle(aggregate).x);int p=aggregate_id_to_particle(aggregate).y;
    ARRAY_VIEW<T> radius=*particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
    if(INTERSECTION::Intersects(object_space_ray,(SPHERE<TV>(particles->X(p),scale*radius(p))),small_number)){
            ray.t_max=object_space_ray.t_max;ray.semi_infinite=object_space_ray.semi_infinite;ray.aggregate_id=aggregate;intersection=true;}
    return intersection;}

    bool Inside(const TV& location) const  PHYSBAM_OVERRIDE
    {for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++) for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++) for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++) {
        POINT_CLOUD<TV>* particles=particles_array(i,j,ij);
        if(particles){
            ARRAY_VIEW<T> radius=*particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
            for(int p=1;p<=particles->array_collection->Size();++p) if(SPHERE<TV>(particles->X(p),scale*radius(p)).Inside(Object_Space_Point(location),small_number)) return true;}}
        return false;} 

    TV Normal(const TV& location,const int aggregate=1) const  PHYSBAM_OVERRIDE
    {POINT_CLOUD<TV>* particles=particles_array(aggregate_id_to_particle(aggregate).x);int p=aggregate_id_to_particle(aggregate).y;
    ARRAY_VIEW<T> radius=*particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
    return (SPHERE<TV>(particles->X(p),scale*radius(p))).Normal(location);}

    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives)const PHYSBAM_OVERRIDE
    {for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++) for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++) for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++) {
        POINT_CLOUD<TV>* particles=particles_array(i,j,ij);
        if(particles){
            ARRAY_VIEW<T> radius=*particles->array_collection->template Get_Array<T>(ATTRIBUTE_ID(15)); // radius is attribute 15
            for(int p=1;p<=particles->array_collection->Size();++p){
                T radius_scalar=scale*radius(p);
                TV radius_vector(radius_scalar,radius_scalar,radius_scalar),world_center=transform.Homogeneous_Times(particles->X(p));
                RANGE<TV> world_space_bounding_box(world_center-radius_vector,world_center+radius_vector);
                int aggregate_id=particle_to_aggregate_id(i,j,ij)(p);
                primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(world_space_bounding_box,this,aggregate_id));}}}}
private:
    void Create_Aggregate_Ids()
    {int index=1;
    for(int i=particles_array.domain.min_corner.x;i<=particles_array.domain.max_corner.x;i++) for(int j=particles_array.domain.min_corner.y;j<=particles_array.domain.max_corner.y;j++) for(int ij=particles_array.domain.min_corner.z;ij<=particles_array.domain.max_corner.z;ij++){
        POINT_CLOUD<TV>* particles=particles_array(i,j,ij);
        if(particles)
            for(int p=1;p<=particles->array_collection->Size();p++){
                particle_to_aggregate_id(i,j,ij).Append(index);aggregate_id_to_particle.Append(PAIR<VECTOR<int,3>,int>(VECTOR<int,3>(i,j,ij),p));index++;}}}

//##################################################################### 
};

}
#endif
