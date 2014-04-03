//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE_HELPER.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/AVERAGING_COLLIDABLE_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE.h>
namespace PhysBAM{

template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_RLE<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE:public ADVECTION<T_GRID,typename T_GRID::SCALAR,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;
public:
    T_GRID_BASED_COLLISION_GEOMETRY& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_RLE<T_GRID,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_RLE<T_GRID,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_RLE<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_RLE<T_GRID,T_FACE_LOOKUP> averaging_collidable;
    ARRAY<bool>& face_velocities_valid_mask;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE(T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,ARRAY<bool>& face_velocities_valid_mask_input)
        :body_list(body_list_input),averaging_collidable(body_list,0),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE()
    {}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,ARRAY<T>& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T>* Z_min,ARRAY<T>* Z_max)
    {assert(!Z_min && !Z_max);
    ARRAY<bool> face_velocities_valid_mask_next(grid.number_of_faces,false);face_velocities_valid_mask_next.Fill(true);
    T_GRID::template Face_Loop<Update_Advection_Equation_Face_Helper>(*this,grid,Z,Z_ghost,face_velocities,dt,face_velocities_valid_mask_next);
    ARRAY<bool>::Exchange_Arrays(face_velocities_valid_mask,face_velocities_valid_mask_next);
    grid.Put_Ghost_Faces(true,face_velocities_valid_mask);}

private:
    struct Update_Advection_Equation_Face_Helper{template<class T_FACE> static void
    Apply(ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE<T_GRID,T_FACE_LOOKUP>& advection,const T_GRID& grid,ARRAY<T>& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& face_velocities,
        const T dt,ARRAY<bool>& face_velocities_valid_mask_next)
    {T_GRID_BASED_COLLISION_GEOMETRY& body_list=advection.body_list;
    int axis=T_FACE::Axis();
    for(T_FACE iterator(grid,0);iterator;iterator++)if(iterator.Both_Cells_Short()){int face=iterator.Face();
        T_BLOCK face_block(grid,iterator.cell1.I());
        if(!face_block || !body_list.Swept_Occupied_Block(face_block)){
            TV sample_location=iterator.X()-dt*advection.averaging.Face_To_Face_Vector(iterator,face_velocities.Nested());
            T_BLOCK block(grid,sample_location);if(!block) continue;
            Z(face)=advection.linear_interpolation.From_Block_Face_Component(axis,block,Z_ghost.Nested().Starting_Point_Face(iterator),sample_location);}
        else{
            TV grid_point_location=iterator.X(),length_and_direction=-dt*advection.averaging_collidable.Face_To_Face_Vector(face_block,iterator,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            T face_velocity;
            if(body_list.Latest_Velocity_Crossover(iterator,dt,face_velocity)){
                face_velocities_valid_mask_next(face)=false;Z(face)=face_velocity;}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(iterator);
                T_BLOCK block(grid,interpolation_point);if(!block) continue;
                Z(face)=advection.linear_interpolation_collidable.From_Block_Face_Component(axis,block,lookup,interpolation_point);
                face_velocities_valid_mask_next(face)=lookup.found_valid_point;}}}}};
public:

    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point))
        return body_list.Object_Velocity(T_BLOCK(body_list.grid,point),body_id,aggregate_id,point)[axis];
    else return default_value;}

    void Average_To_Invalidated_Face(const T_GRID& grid,ARRAY<T>& values)
    {ARRAY<PAIR<int,bool> > invalid_indices;T_GRID::template Face_Loop<Find_Invalid_Faces>(grid,face_velocities_valid_mask,invalid_indices);
    const ARRAY<VECTOR<bool,T_GRID::dimension> >& face_neighbors_visible=body_list.face_neighbors_visible;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& face_neighbors=grid.Short_Face_Neighbors();
    const ARRAY<TV>& face_locations=grid.Short_Face_Locations();
    bool done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int face=invalid_indices(k).x;T sum=0;int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                int neighbor1=face_neighbors(face)(2*axis-1),neighbor2=face_neighbors(face)(2*axis);
                if(neighbor1 && face_neighbors_visible(neighbor1)(axis) && face_velocities_valid_mask(neighbor1)){sum+=values(neighbor1);count++;}
                if(neighbor2 && face_neighbors_visible(face)(axis) && face_velocities_valid_mask(neighbor2)){sum+=values(neighbor2);count++;}}
            if(count){values(face)=sum/count;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){face_velocities_valid_mask(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
    done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int face=invalid_indices(k).x;T sum=0;int count=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                int neighbor1=face_neighbors(face)(2*axis-1),neighbor2=face_neighbors(face)(2*axis);
                if(neighbor1){
                    if(face_neighbors_visible(neighbor1)(axis)){if(face_velocities_valid_mask(neighbor1)){sum+=values(neighbor1);count++;}}
                    else{sum+=Compute_Revalidation_Value(grid.Face_Axis(face),face_locations(face),face_locations(neighbor1),values(face),T());count++;}}
                if(neighbor2){
                    if(face_neighbors_visible(face)(axis)){if(face_velocities_valid_mask(neighbor2)){sum+=values(neighbor2);count++;}}
                    else{sum+=Compute_Revalidation_Value(grid.Face_Axis(face),face_locations(face),face_locations(neighbor2),values(face),T());count++;}}}
            if(count){values(face)=sum/(T)count;invalid_indices(k).y=true;done=false;}
            /* TODO: cell revalidation has else values(cell)=default_value here*/}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){face_velocities_valid_mask(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}}

private:
    struct Find_Invalid_Faces{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<bool>& face_valid,ARRAY<PAIR<int,bool> >& invalid_faces)
    {for(T_FACE face(grid,0);face;face++)if(face.Both_Cells_Short() && !face_valid(face.Face())) invalid_faces.Append(PAIR<int,bool>(face.Face(),false));}};

//#####################################################################
};
}
#endif
#endif
