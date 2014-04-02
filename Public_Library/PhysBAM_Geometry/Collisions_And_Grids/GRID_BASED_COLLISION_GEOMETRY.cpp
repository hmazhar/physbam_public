//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/RIGID_GEOMETRY_RASTERIZATION_DYADIC.h>
#elif !COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/RIGID_GEOMETRY_RASTERIZATION_DYADIC.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/RIGID_GEOMETRY_RASTERIZATION_RLE.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/RIGID_GEOMETRY_RASTERIZATION_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
GRID_BASED_COLLISION_GEOMETRY(T_GRID& grid_input)
    :grid(grid_input),number_of_ghost_cells(3)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
~GRID_BASED_COLLISION_GEOMETRY()
{}
//##################################################################### 
// Function Add_Bodies
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Add_Bodies(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection)
{
    for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i))
        collision_geometry_collection.Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(rigid_geometry_collection.Rigid_Geometry(i)),i,true);
}
//##################################################################### 
// Function Rasterize_Objects
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Rasterize_Objects()
{
    objects_in_cell.Reset(grid,number_of_ghost_cells);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++)
        if(collision_geometry_collection.bodies(i) && collision_geometry_collection.Is_Active(i) && collision_geometry_collection.bodies(i)->active)
            RASTERIZATION::Rasterize_Object(*collision_geometry_collection.bodies(i),grid,objects_in_cell,i);
}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T min_time=FLT_MAX;bool collision=false;T current_hit_time;T_WEIGHTS current_weights;int current_simplex_id;
    if(!objects){for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Earliest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id) && current_hit_time < min_time){
            min_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;}}
    else for(int k=1;k<=objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Earliest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id) && current_hit_time < min_time){
            min_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;}}
    return collision;
}
//##################################################################### 
// Function Latest_Crossover
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Latest_Crossover(const TV& start_X,const TV& end_X,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,TV& initial_hit_point,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T hit_time;T_WEIGHTS weights;POINT_SIMPLEX_COLLISION_TYPE returned_collision_type;
    bool crossover=Latest_Simplex_Crossover(start_X,end_X,dt,hit_time,weights,body_id,simplex_id,returned_collision_type,objects);
    if(crossover){
        T_SIMPLEX initial_simplex=collision_geometry_collection(body_id).World_Space_Simplex(simplex_id);
        initial_hit_point=initial_simplex.Point_From_Barycentric_Coordinates(weights);}
    return crossover;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type,
    const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;POINT_SIMPLEX_COLLISION_TYPE collision_type;
    T max_time=-FLT_MAX;bool collision=false;T current_hit_time;T_WEIGHTS current_weights;int current_simplex_id;
    if(!objects){for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Latest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id,collision_type) && current_hit_time > max_time){
            max_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;
            returned_collision_type=collision_type;}}
    else for(int k=1;k<=objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && 
        collision_geometry_collection.bodies(i)->Latest_Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,current_simplex_id,collision_type) && current_hit_time > max_time){
            max_time=hit_time=current_hit_time;weights=current_weights;body_id=i;simplex_id=current_simplex_id;collision=true;
            returned_collision_type=collision_type;}}
    return collision;
}
//##################################################################### 
// Function Any_Simplex_Crossover
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    if(!objects){for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && collision_geometry_collection.bodies(i)->Any_Simplex_Crossover(start_X,end_X,dt)) return true;}
    else for(int k=1;k<=objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);if(Is_Active(i) && collision_geometry_collection.bodies(i)->active && collision_geometry_collection.bodies(i)->Any_Simplex_Crossover(start_X,end_X,dt)) return true;}
    return false;
}
//##################################################################### 
// Function Update_Intersection_Acceleration_Structures
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1,const int state2)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Update_Intersection_Acceleration_Structures(use_swept_simplex_hierarchy,state1,state2);
}
//##################################################################### 
// Function Get_Body_Penetration
//##################################################################### 
// normal is flipped to ensure that start_phi is positive
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,TV& body_velocity,
    const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T hit_time=FLT_MAX;T current_hit_time;int current_simplex_id;T current_start_phi,current_end_phi;TV current_end_body_normal,current_body_velocity;
    if(objects) for(int k=1;k<=objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);
        if(collision_geometry_collection.bodies(i)->Get_Body_Penetration(start_X,end_X,contour_value,dt,current_hit_time,current_simplex_id,current_start_phi,current_end_phi,current_end_body_normal,
             current_body_velocity) && current_hit_time<hit_time){
             body_id=i;hit_time=current_hit_time;
             simplex_id=current_simplex_id;start_phi=current_start_phi;end_phi=current_end_phi;end_body_normal=current_end_body_normal;body_velocity=current_body_velocity;}}
    else for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.Size();i++) 
        if(Is_Active(i) && collision_geometry_collection.bodies(i)->Get_Body_Penetration(start_X,end_X,contour_value,dt,current_hit_time,current_simplex_id,current_start_phi,current_end_phi,current_end_body_normal,
                current_body_velocity) && current_hit_time<hit_time){
            body_id=i;hit_time=current_hit_time;
            simplex_id=current_simplex_id;start_phi=current_start_phi;end_phi=current_end_phi;end_body_normal=current_end_body_normal;body_velocity=current_body_velocity;}
    return hit_time<FLT_MAX;
}
//##################################################################### 
// Function Push_Out_Point
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Push_Out_Point(TV& X,const T collision_distance,const bool check_particle_crossover,bool& particle_crossover,const ARRAY<COLLISION_GEOMETRY_ID>* objects) const
{
    T distance=FLT_MAX,current_distance;TV X_old=X;
    if(objects) for(int k=1;k<=objects->m;k++){COLLISION_GEOMETRY_ID i=(*objects)(k);TV current_X=X_old;
        if(collision_geometry_collection.bodies(i)->Push_Out_Point(current_X,collision_distance,current_distance) && current_distance<distance){X=current_X;distance=current_distance;}}
    else for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.Size();i++){TV current_X=X_old;
        if(Is_Active(i) && collision_geometry_collection.bodies(i)->Push_Out_Point(current_X,collision_distance,current_distance) && current_distance<distance){X=current_X;distance=current_distance;}}
    if(distance<FLT_MAX){COLLISION_GEOMETRY_ID body_id;int simplex_id;TV intersection_point;
        if(check_particle_crossover) particle_crossover=collision_geometry_collection.Intersection_Between_Points(X_old,X,body_id,simplex_id,intersection_point,objects);
        return true;}
    return false;
}
//##################################################################### 
// Function Occupied_Block
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Occupied_Block(const T_BLOCK& block) const
{
    return occupied_blocks(block.Block());
}
//##################################################################### 
// Function Swept_Occupied_Block
//##################################################################### 
template<class T_GRID> bool GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Swept_Occupied_Block(const T_BLOCK& block) const
{
    return swept_occupied_blocks(block.Block());
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//##################################################################### 
// Function Read_State
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Read_State(TYPED_ISTREAM& input,const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Read_State(input,state_index);
}
//##################################################################### 
// Function Write_State
//##################################################################### 
template<class T_GRID> void GRID_BASED_COLLISION_GEOMETRY<T_GRID>::
Write_State(TYPED_OSTREAM& output,const int state_index) const
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Write_State(output,state_index);
}
#endif
//#####################################################################
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<float,1> > >;
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<float,2> > >;
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<OCTREE_GRID<float> >;
template class GRID_BASED_COLLISION_GEOMETRY<QUADTREE_GRID<float> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<RLE_GRID_2D<float> >;
template class GRID_BASED_COLLISION_GEOMETRY<RLE_GRID_3D<float> >;
#endif
#if COMPILE_WITH_BINTREE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<BINTREE_GRID<float> >;
#endif

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<double,1> > >;
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<double,2> > >;
template class GRID_BASED_COLLISION_GEOMETRY<GRID<VECTOR<double,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<OCTREE_GRID<double> >;
template class GRID_BASED_COLLISION_GEOMETRY<QUADTREE_GRID<double> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<RLE_GRID_2D<double> >;
template class GRID_BASED_COLLISION_GEOMETRY<RLE_GRID_3D<double> >;
#endif
#if COMPILE_WITH_BINTREE_SUPPORT
template class GRID_BASED_COLLISION_GEOMETRY<BINTREE_GRID<double> >;
#endif
#endif
