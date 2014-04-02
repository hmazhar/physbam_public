#if 0
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM  
//##################################################################### 
#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T_FACE_LOOKUP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<T_GRID,T_FACE_LOOKUP>::
ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM(T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
    :body_list(body_list_input),averaging_collidable(body_list,0),face_velocities_valid_mask(dynamic_cast<T_FACE_ARRAYS_SLIP_BOOL&>(face_velocities_valid_mask_input))
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T_FACE_LOOKUP> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<T_GRID,T_FACE_LOOKUP>::
~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM()
{
}
//#####################################################################
// Function Update_Advection_Equation_Face_Lookup
//#####################################################################
template<class T_GRID,class T_FACE_LOOKUP> void ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<T_GRID,T_FACE_LOOKUP>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    T_FACE_ARRAYS_SLIP_BOOL face_velocities_valid_mask_next(grid,3,false);
    ARRAY<T,SIDED_FACE_INDEX<TV::dimension> >& Z_cast=dynamic_cast<ARRAY<T,SIDED_FACE_INDEX<TV::dimension> >&>(Z);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
    {
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        for(int side=1;side<=2;side++)
        {
            //body_list.side=side; // TODO: make this nicer/avoid
            if(!body_list.Occupied_Face_Center(iterator))
            {
                TV grid_point_location=iterator.Location(),length_and_direction=-dt*averaging.Face_To_Face_Vector(grid,axis,face,face_velocities.Nested()),
                    interpolation_point=grid_point_location+length_and_direction;
                Z_cast(side,axis,face)=linear_interpolation.Clamped_To_Array_Face_Component(axis,grid,Z_ghost.Nested().Starting_Point_Face(axis,face),interpolation_point);
                face_velocities_valid_mask_next(side,axis,face)=true;
            }
            else
            {
                TV_INT adjacent_cell_center=side&1?iterator.First_Cell_Index():iterator.Second_Cell_Index();
                TV_INT other_cell_center=side&1?iterator.Second_Cell_Index():iterator.First_Cell_Index();
                TV cell_center_location=grid.Center(adjacent_cell_center);
                TV grid_point_location=iterator.Location();
                TV offset=grid_point_location-cell_center_location;
                //if(side==2 && axis==1 && face(1)==1 && face(2)==4)
                //    DEBUG_UTILITIES::Debug_Breakpoint();
                TV length_and_direction=-dt*averaging_collidable.Face_To_Face_Vector(grid,side,axis,face,face_velocities);

                // properly, what this needs to check is whether the segment between the face center and the cell center is crossed over.  It should be sufficient to check whether
                // it's presently occluded as well as crossover for the two points.
                //if(side==1 && axis==2 && face(1)==3 && face(2)==5)
                //    DEBUG_UTILITIES::Debug_Breakpoint();
                TV velocity;
                if(body_list.Latest_Cell_Crossover_And_Velocity(adjacent_cell_center,dt,velocity))
                    length_and_direction=-dt*velocity;
                TV interpolation_point=grid_point_location+length_and_direction;
                length_and_direction=interpolation_point-cell_center_location;
                RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                if(RAY<TV>::Create_Non_Degenerate_Ray(cell_center_location,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id)) 
                {
                    int aggregate_id=0;
                    body_list.collision_geometry_collection.Intersection_Between_Points(cell_center_location,cell_center_location+length_and_direction,body_id,aggregate_id,interpolation_point);
                    Z_cast(side,axis,face)=body_list.Object_Velocity(body_id,aggregate_id,interpolation_point)[axis];
                }
                else{
                    const typename T_FACE_LOOKUP::LOOKUP& lookup=Z_ghost.Starting_Point_Face(side,axis,face);
                    lookup.Set_Reference_Point(interpolation_point);
                    Z_cast(side,axis,face)=linear_interpolation_collidable.interpolation.From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<T_GRID>(grid,interpolation_point,lookup.Number_Of_Ghost_Cells()),lookup,interpolation_point);
                    lookup.Clear_Reference_Point();
                    face_velocities_valid_mask_next(side,axis,face)=lookup.found_valid_point;
                }
            }
        }
    }
    T_FACE_ARRAYS_SLIP_BOOL::Exchange_Arrays(face_velocities_valid_mask,face_velocities_valid_mask_next);
    // ghost values should always be valid
    for(int side=1;side<=2;side++) for(int axis=1;axis<=TV::dimension;axis++) grid.Put_Ghost(true,face_velocities_valid_mask.Component(side,axis),3);
}
//#####################################################################
// Function Average_To_Invalidated_Face
//#####################################################################
template<class T_GRID,class T_FACE_LOOKUP> void ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<T_GRID,T_FACE_LOOKUP>::
Average_To_Invalidated_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_values,T_FACE_ARRAYS_BOOL* faces_not_to_revalidate)
{
    return;
    PHYSBAM_FATAL_ERROR("Not doing this");
    // average values collision aware in Gauss-Jacobi fashion
    ARRAY<T,SIDED_FACE_INDEX<TV::dimension> >& face_values_slip=dynamic_cast<ARRAY<T,SIDED_FACE_INDEX<TV::dimension> >&>(face_values);
    typename TV::template REBIND<ARRAY<TRIPLE<int,TV_INT,bool> > >::TYPE face_invalid_indices; // index and bool true if entry has been validated on iteration
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        if(!face_velocities_valid_mask(1,iterator.Axis(),iterator.Face_Index()))
            face_invalid_indices[iterator.Axis()].Append(TRIPLE<int,TV_INT,bool>(1,iterator.Face_Index(),false));
        if(!face_velocities_valid_mask(2,iterator.Axis(),iterator.Face_Index()))
            face_invalid_indices[iterator.Axis()].Append(TRIPLE<int,TV_INT,bool>(2,iterator.Face_Index(),false));}

    for(int side=1;side<=2;side++) for(int axis=1;axis<=TV::dimension;axis++) grid.Put_Ghost(false,face_velocities_valid_mask.Component(side,axis),3);

    T_ARRAYS_BOOL_DIMENSION& cell_neighbors_visible=body_list.cell_neighbors_visible;
    int expensive_count=0;

    for(int arrays_axis=1;arrays_axis<=T_GRID::dimension;arrays_axis++){
        ARRAY<TRIPLE<int,TV_INT,bool> >& invalid_indices=face_invalid_indices[arrays_axis];
        // T_ARRAYS_BOOL_DIMENSION& neighbors_visible=body_list.face_neighbors_visible.Component(arrays_axis);

        bool done=false;
        //T_ARRAYS_BOOL::Put_Ghost(false,valid_points,grid,3); // don't average from boundaries

        while(!done){
            done=true;
            for(int k=1;k<=invalid_indices.m;k++){ 
                if(invalid_indices(k).z || (faces_not_to_revalidate && (*faces_not_to_revalidate)(arrays_axis,invalid_indices(k).y))) continue; // if we're already valid just keep going
                T sum=0;int count=0;
                int side=invalid_indices(k).x; // across cell always visible, neighboring cell visible iff cells can see each other
                TV_INT face_index=invalid_indices(k).y;
                int other_side=side%2+1;
                if(face_velocities_valid_mask(other_side,arrays_axis,face_index) && cell_neighbors_visible(face_index-TV_INT::Axis_Vector(arrays_axis))(arrays_axis)){
                    face_values_slip(side,arrays_axis,face_index)=face_values_slip(other_side,arrays_axis,face_index);
                }else{
                    if(side==1 && face_velocities_valid_mask(2,arrays_axis,face_index-TV_INT::Axis_Vector(arrays_axis))){
                        sum+=face_values_slip(2,arrays_axis,face_index-TV_INT::Axis_Vector(arrays_axis));
                        count++;}
                    if(side==2 && face_velocities_valid_mask(1,arrays_axis,face_index+TV_INT::Axis_Vector(arrays_axis))){
                        sum+=face_values_slip(1,arrays_axis,face_index+TV_INT::Axis_Vector(arrays_axis));
                        count++;}
                    TV_INT reference_cell=(side==1?face_index-TV_INT::Axis_Vector(arrays_axis):face_index);
                    for(int axis=1;axis<=T_GRID::dimension;axis++){
                        if(axis!=arrays_axis){
                            TV_INT min_face=face_index-TV_INT::Axis_Vector(axis),max_face=face_index+TV_INT::Axis_Vector(axis);
                            // only "right" neighbor visibility is coherent?
                            if(cell_neighbors_visible(reference_cell)(axis) && face_velocities_valid_mask(side,arrays_axis,max_face)){sum+=face_values_slip(side,arrays_axis,max_face);count++;}
                            if(cell_neighbors_visible(reference_cell-TV_INT::Axis_Vector(axis))(axis) && face_velocities_valid_mask(side,arrays_axis,min_face)){
                                sum+=face_values_slip(side,arrays_axis,min_face);count++;}}}
                    if(count){
                        face_values_slip(side,arrays_axis,invalid_indices(k).y)=sum/(T)count;invalid_indices(k).z=true;done=false;
                        if(!face_velocities_valid_mask(other_side,arrays_axis,face_index) && cell_neighbors_visible(face_index-TV_INT::Axis_Vector(arrays_axis))(arrays_axis)){
                            // find the other version of this face and set it too if it's invalid.  This is expensive; try not to have it happen much.
                            //LOG::cout<<"Doing expensive check at side "<<side<<" axis "<<arrays_axis<<" face index "<<face_index<<std::endl;
                            face_values_slip(other_side,arrays_axis,face_index)=face_values_slip(side,arrays_axis,face_index);
                            expensive_count++;
                            for(int kk=k+1;kk<=invalid_indices.m;kk++)
                                if(invalid_indices(kk).x==other_side && invalid_indices(kk).y==face_index)
                                    invalid_indices(kk).z=true;}}}}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).z){
                face_velocities_valid_mask(invalid_indices(k).x,arrays_axis,invalid_indices(k).y)=true;invalid_indices.Remove_Index_Lazy(k);}}

        // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
        done=false;
        
        while(!done){
            done=true;
            for(int k=1;k<=invalid_indices.m;k++){ 
                T sum=0;int count=0;
                int side=invalid_indices(k).x; // across cell always visible, neighboring cell visible iff cells can see each other
                TV_INT face_index=invalid_indices(k).y;
                int other_side=side%2+1;
                if(face_velocities_valid_mask(other_side,arrays_axis,face_index) && cell_neighbors_visible(face_index-TV_INT::Axis_Vector(arrays_axis))(arrays_axis)){
                    face_values_slip(side,arrays_axis,face_index)=face_values_slip(other_side,arrays_axis,face_index);
                }else{
                    if(side==1 && face_velocities_valid_mask(2,arrays_axis,face_index-TV_INT::Axis_Vector(arrays_axis))){
                        sum+=face_values_slip(2,arrays_axis,face_index-TV_INT::Axis_Vector(arrays_axis));
                        count++;}
                    if(side==2 && face_velocities_valid_mask(1,arrays_axis,face_index+TV_INT::Axis_Vector(arrays_axis))){
                        sum+=face_values_slip(1,arrays_axis,face_index+TV_INT::Axis_Vector(arrays_axis));
                        count++;}
                    TV_INT reference_cell=(side==1?face_index-TV_INT::Axis_Vector(arrays_axis):face_index);
                    for(int axis=1;axis<=T_GRID::dimension;axis++){
                        if(axis!=arrays_axis){
                            TV_INT min_face=face_index-TV_INT::Axis_Vector(axis),max_face=face_index+TV_INT::Axis_Vector(axis);
                            // only "right" neighbor visibility is coherent?
                            if(cell_neighbors_visible(reference_cell)(axis) && face_velocities_valid_mask(side,arrays_axis,max_face)){sum+=face_values_slip(side,arrays_axis,max_face);count++;}
                            if(cell_neighbors_visible(reference_cell-TV_INT::Axis_Vector(axis))(axis) && face_velocities_valid_mask(side,arrays_axis,min_face)){
                                sum+=face_values_slip(side,arrays_axis,min_face);count++;}}}
                    if(!count){
                        sum=face_values_slip(side,arrays_axis,face_index);
                        count=1;}
                    if(count){face_values_slip(side,arrays_axis,invalid_indices(k).y)=sum/(T)count;invalid_indices(k).z=true;done=false;}}}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).z){
                face_velocities_valid_mask(invalid_indices(k).x,arrays_axis,invalid_indices(k).y)=true;invalid_indices.Remove_Index_Lazy(k);}}}

    // set valid for future advection
    for(int side=1;side<=2;side++) for(int axis=1;axis<=TV::dimension;axis++) grid.Put_Ghost(true,face_velocities_valid_mask.Component(side,axis),3);
    
    // deal with faces which may have become the same (e.g. an object moves away from them)
    /*for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();const TV_INT& face_index=iterator.Face_Index();
        if(face_velocities_valid_mask(1,axis,face_index) && face_velocities_valid_mask(2,axis,face_index)
            && face_values(1,axis,face_index) != face_values(2,axis,face_index)){
            face_values(1,axis,face_index)=face_values(2,axis,face_index)=(T).5*(face_values(1,axis,face_index)+face_values(2,axis,face_index));
            }}*/
}
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
#endif
#endif
