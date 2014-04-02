//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM  
//##################################################################### 
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM:public ADVECTION<T_GRID,typename T_GRID::SCALAR,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;
public:
    T_GRID_BASED_COLLISION_GEOMETRY& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<T_GRID,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_UNIFORM<T_GRID,T_FACE_LOOKUP> averaging_collidable;
    T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM(T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
        :body_list(body_list_input),averaging_collidable(body_list,0),
         face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM()
    {}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {T_FACE_ARRAYS_BOOL face_velocities_valid_mask_next(grid,3,false);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        T_VECTOR_INT face=iterator.Face_Index();int axis=iterator.Axis();
        if(!body_list.Swept_Occupied_Face_Center(iterator)){
            TV interpolation_point=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities.Nested());
            Z.Component(axis)(face)=linear_interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Nested().Starting_Point_Face(axis,face),interpolation_point);
            face_velocities_valid_mask_next.Component(axis)(face)=true;
            if(Z_min && Z_max){
                VECTOR<T,2> extrema=linear_interpolation.Extrema_Clamped_To_Array_Face_Component(axis,grid,Z_min_ghost->Nested().Starting_Point_Face(axis,face),
                    Z_max_ghost->Nested().Starting_Point_Face(axis,face),interpolation_point);
                Z_min->Component(axis)(face)=extrema.x;Z_max->Component(axis)(face)=extrema.y;}}
        else{
            T face_velocity=0; // avoid uninitialization warning
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*averaging_collidable.Face_To_Face_Vector(grid,axis,face,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            if(body_list.Latest_Velocity_Crossover(axis,face,dt,face_velocity)){
                face_velocities_valid_mask_next.Component(axis)(face)=false;Z.Component(axis)(face)=face_velocity;
                if(Z_min && Z_max){Z_min->Component(axis)(face)=Z_max->Component(axis)(face)=face_velocity;}}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id)) 
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(axis,face);
                Z.Component(axis)(face)=linear_interpolation_collidable.Clamped_To_Array_Face_Component(axis,grid,lookup,interpolation_point);
                face_velocities_valid_mask_next.Component(axis)(face)=lookup.found_valid_point;
                if(Z_min && Z_max){
                    const typename T_FACE_LOOKUP::LOOKUP &Z_min_lookup=Z_min_ghost->Starting_Point_Face(axis,face),&Z_max_lookup=Z_max_ghost->Starting_Point_Face(axis,face);
                    VECTOR<T,2> extrema=linear_interpolation_collidable.Extrema_Clamped_To_Array_Face_Component(axis,grid,Z_min_lookup,Z_max_lookup,interpolation_point);
                    Z_min->Component(axis)(face)=extrema.x;Z_max->Component(axis)(face)=extrema.y;}}}}
    T_FACE_ARRAYS_BOOL::Exchange_Arrays(face_velocities_valid_mask,face_velocities_valid_mask_next);
    // ghost values should always be valid
    for(int axis=1;axis<=T_GRID::dimension;axis++) grid.Get_Face_Grid(axis).Put_Ghost(true,face_velocities_valid_mask.Component(axis),3);}

    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point)) return body_list.Object_Velocity(body_id,aggregate_id,point)[axis];
    else return default_value;}

    void Average_To_Invalidated_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_values)
    {// average values collision aware in Gauss-Jacobi fashion
    typename TV::template REBIND<ARRAY<PAIR<T_VECTOR_INT,bool> > >::TYPE face_invalid_indices; // index and bool true if entry has been validated on iteration
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(!face_velocities_valid_mask.Component(iterator.Axis())(iterator.Face_Index())) 
        face_invalid_indices[iterator.Axis()].Append(PAIR<T_VECTOR_INT,bool>(iterator.Face_Index(),false));
    
    for(int arrays_axis=1;arrays_axis<=T_GRID::dimension;arrays_axis++){
        ARRAY<PAIR<T_VECTOR_INT,bool> >& invalid_indices=face_invalid_indices[arrays_axis];
        T_ARRAYS_BOOL_DIMENSION& neighbors_visible=body_list.face_neighbors_visible.Component(arrays_axis);
        T_ARRAYS_BOOL& valid_points=face_velocities_valid_mask.Component(arrays_axis);T_ARRAYS_BASE& values=face_values.Component(arrays_axis);
        
        bool done=false;
        grid.Put_Ghost(false,valid_points,3); // don't average from boundaries

        while(!done){
            done=true;
            for(int k=1;k<=invalid_indices.m;k++){ 
                T sum=0;int count=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    T_VECTOR_INT min_face=invalid_indices(k).x-T_VECTOR_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+T_VECTOR_INT::Axis_Vector(axis);
                    if(neighbors_visible(min_face)(axis) && valid_points(min_face)){sum+=values(min_face);count++;}
                    if(neighbors_visible(invalid_indices(k).x)(axis) && valid_points(max_face)){sum+=values(max_face);count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}

        // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
        done=false;
        while(!done){
            done=true;
            for(int k=1;k<=invalid_indices.m;k++){ 
                T sum=0;int count=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    T_VECTOR_INT min_face=invalid_indices(k).x-T_VECTOR_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+T_VECTOR_INT::Axis_Vector(axis);
                    if(neighbors_visible(min_face)(axis)){if(valid_points(min_face)){sum+=values(min_face);count++;}}
                    else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(min_face),values(invalid_indices(k).x),T());count++;}
                    if(neighbors_visible(invalid_indices(k).x)(axis)){if(valid_points(max_face)){sum+=values(max_face);count++;}}
                    else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(max_face),values(invalid_indices(k).x),T());count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}
                else values(invalid_indices(k).x)=T();}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
        grid.Put_Ghost(true,valid_points,3);}} // set valid for future advection

//#####################################################################
};
}
#endif
