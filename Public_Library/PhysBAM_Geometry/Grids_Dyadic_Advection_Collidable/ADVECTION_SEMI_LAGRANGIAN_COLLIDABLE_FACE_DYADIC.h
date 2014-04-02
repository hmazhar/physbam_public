//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC  
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/AVERAGING_COLLIDABLE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC.h>
namespace PhysBAM{
    
template<class T_GRID,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC:public ADVECTION<T_GRID,typename T_GRID::SCALAR,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;typedef typename T_GRID::CELL T_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;
public:    
    T_GRID_BASED_COLLISION_GEOMETRY& body_list;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_DYADIC<T_GRID,T,T_FACE_LOOKUP> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_DYADIC<T_GRID,T,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_DYADIC<T_GRID,typename T_FACE_LOOKUP::NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_DYADIC<T_GRID,T_FACE_LOOKUP> averaging_collidable;
    T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC(T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
        :body_list(body_list_input),averaging_collidable(body_list,0),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC()
    {}
 
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {assert(!Z_min && !Z_max);
    ARRAY<TV> Z_node_ghost(grid.number_of_nodes,false);LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Faces_To_Nodes(grid,Z_ghost.Raw_Data(),Z_node_ghost);
    ARRAY<TV> V_node(grid.number_of_nodes,false);LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Faces_To_Nodes(grid,face_velocities.Raw_Data(),V_node);
    ARRAY<bool> face_velocities_valid_mask_next(grid.number_of_faces,false);
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){int face=iterator.Face_Index(),axis=iterator.Axis();
        if(!body_list.Swept_Occupied_Face_Center(iterator)){
            TV sample_location=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,iterator.Face_Index(),face_velocities.Nested(),V_node);
            T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(iterator.Deepest_Cell(),sample_location);
            if(base_cell) Z(face)=linear_interpolation.From_Block_Face_Component(axis,grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),Z_ghost.Nested().Starting_Point_Face(face),sample_location);
            else{
                T_CELL* leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(iterator.Deepest_Cell(),sample_location);
                Z(face)=linear_interpolation.From_Leaf_Cell_Face_Component(axis,grid,leaf_cell,Z_ghost.Raw_Data(),Z_node_ghost,sample_location);}
            face_velocities_valid_mask_next(iterator.Face_Index())=true;}
        else{
            TV grid_point_location=iterator.Location(),length_and_direction=-dt*averaging_collidable.Face_To_Face_Vector(grid,face,face_velocities),
                interpolation_point=grid_point_location+length_and_direction;
            T face_velocity;
            if(body_list.Latest_Velocity_Crossover(face,dt,face_velocity)){
                face_velocities_valid_mask_next(iterator.Face_Index())=false;Z(iterator.Face_Index())=face_velocity;}
            else{
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(grid_point_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                    interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(iterator.Face_Index());
                T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(iterator.Deepest_Cell(),interpolation_point);assert(base_cell);
                BLOCK_DYADIC<T_GRID> block(grid,base_cell);
                Z(face)=linear_interpolation_collidable.From_Block_Face_Component(axis,grid,block,lookup,interpolation_point);
                face_velocities_valid_mask_next(face)=lookup.found_valid_point;}}}
    ARRAY<bool>::Exchange_Arrays(face_velocities_valid_mask,face_velocities_valid_mask_next);
    // ghost values should always be valid
    Put_Ghost(true,face_velocities_valid_mask);}
    
    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point)) return body_list.Object_Velocity(body_id,aggregate_id,point)[axis];
    else return default_value;}

    void Average_To_Invalidated_Face(const T_GRID& grid,ARRAY<T>& values)
    {ARRAY<PAIR<int,bool> > invalid_indices;ARRAY<VECTOR<bool,T_GRID::dimension> >& face_neighbors_visible=body_list.face_neighbors_visible;
    Put_Ghost(false,face_velocities_valid_mask); // don't revalidate from ghost regions
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()) if(!face_velocities_valid_mask(iterator.Face_Index())) 
        invalid_indices.Append(PAIR<int,bool>(iterator.Face_Index(),false));
    bool done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int face=invalid_indices(k).x;T sum=0;int count=0;
            FACE_ITERATOR iterator(grid,face);
            for(int i=0;i<T_GRID::dimension;i++){
                int min_other_face=iterator.Minimal_Neighbor(2*i),max_other_face=iterator.Minimal_Neighbor(2*i+1);
                if(min_other_face && face_neighbors_visible(min_other_face)(i+1) && face_velocities_valid_mask(min_other_face)){sum+=values(min_other_face);count++;}
                if(max_other_face && face_neighbors_visible(face)(i+1) && face_velocities_valid_mask(max_other_face)){sum+=values(max_other_face);count++;}}
            if(count){values(iterator.Face_Index())=sum/count;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){face_velocities_valid_mask(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
    done=false;
    while(!done){
        done=true;
        for(int k=1;k<=invalid_indices.m;k++){
            int face=invalid_indices(k).x;T sum=0;int count=0;
            FACE_ITERATOR iterator(grid,face);
            for(int i=0;i<T_GRID::dimension;i++){
                int min_other_face=iterator.Minimal_Neighbor(2*i),max_other_face=iterator.Minimal_Neighbor(2*i+1);
                if(min_other_face){
                    if(face_neighbors_visible(min_other_face)(i+1)&&face_velocities_valid_mask(min_other_face)){sum+=values(min_other_face);count++;}
                    else{sum+=Compute_Revalidation_Value(iterator.Axis(),grid.Face_Location(invalid_indices(k).x),grid.Face_Location(min_other_face),values(face),T());count++;}}
                else{sum+=Compute_Revalidation_Value(iterator.Axis(),grid.Face_Location(invalid_indices(k).x),grid.Face_Location(invalid_indices(k).x),values(invalid_indices(k).x),T());count++;}
                if(max_other_face){
                    if(face_neighbors_visible(face)(i+1) && face_velocities_valid_mask(max_other_face)){sum+=values(max_other_face);count++;}
                    else{sum+=Compute_Revalidation_Value(iterator.Axis(),grid.Face_Location(invalid_indices(k).x),grid.Face_Location(max_other_face),values(face),T());count++;}}}
            if(count){values(iterator.Face_Index())=sum/count;invalid_indices(k).y=true;done=false;}}
        if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){face_velocities_valid_mask(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
    Put_Ghost(true,face_velocities_valid_mask);} // restore ghost values to valid

private:
    void Put_Ghost(const bool constant,ARRAY<bool>& array)
    {for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(body_list.grid,body_list.grid.Map_Ghost_Faces());iterator.Valid();iterator.Next()) array(iterator.Face_Index())=constant;}

//##################################################################### 
};
}
#endif
#endif
