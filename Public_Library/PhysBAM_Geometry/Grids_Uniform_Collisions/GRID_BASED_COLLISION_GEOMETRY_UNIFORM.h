//#####################################################################
// Copyright 2005-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY_UNIFORM
//#####################################################################
#ifndef __GRID_BASED_COLLISION_GEOMETRY_UNIFORM__
#define __GRID_BASED_COLLISION_GEOMETRY_UNIFORM__

#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
namespace PhysBAM{

template <class T_GRID>
class GRID_BASED_COLLISION_GEOMETRY_UNIFORM:public GRID_BASED_COLLISION_GEOMETRY<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<bool,TV::dimension> TV_BOOL;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:

    typedef GRID_BASED_COLLISION_GEOMETRY<T_GRID> BASE;
    using BASE::collision_geometry_collection;using BASE::grid;using BASE::collision_thickness;using BASE::objects_in_cell;using BASE::cell_neighbors_visible;
    using BASE::face_neighbors_visible;using BASE::Get_Body_Penetration;using BASE::occupied_blocks;using BASE::swept_occupied_blocks;using BASE::Is_Active;

    const T_ARRAYS_BOOL* outside_fluid;
    bool use_collision_face_neighbors;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM(T_GRID& grid_input);
    ~GRID_BASED_COLLISION_GEOMETRY_UNIFORM();

    bool Occupied_Cell_Center(const TV_INT& cell_index) const
    {return occupied_blocks(cell_index);} // this will check the occupied block that is left,bottom,front of the current cell.

    bool Swept_Occupied_Cell_Center(const TV_INT& cell_index) const
    {return swept_occupied_blocks(cell_index);} // this will check the occupied block that is left,bottom,front of the current cell.

    bool Occupied_Face_Center(FACE_ITERATOR& iterator) const
    {return occupied_blocks(iterator.Face_Node_Index(1));}

    bool Swept_Occupied_Face_Center(FACE_ITERATOR& iterator) const
    {return swept_occupied_blocks(iterator.Face_Node_Index(1));}

    bool Inside_Any_Simplex_Of_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;TV_INT cell_index=grid.Clamp_To_Cell(location);objects_in_cell.Get_Objects_For_Cell(cell_index,objects);if(!objects.m) return false;
    return collision_geometry_collection.Inside_Any_Simplex_Of_Any_Body(location,body_id,aggregate_id,&objects);}

    // TODO(jontg): Slow, but it works...
    bool Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const
    {return collision_geometry_collection.Inside_Any_Body(location,collision_thickness*(T).5,body_id);}

    bool Implicit_Geometry_Lazy_Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const
    {return collision_geometry_collection.Implicit_Geometry_Lazy_Inside_Any_Body(location,body_id);}

    bool Closest_Non_Intersecting_Point_Of_Any_Body(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id) const
    {// since the ray is a general ray, we don't get specific id's.
    return collision_geometry_collection.Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(ray,body_id);}

    bool Cell_Center_Visible_From_Face(const TV_INT& cell,const int axis,const TV_INT& face_index) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cell(cell,objects);if(!objects.m) return true;
    return !collision_geometry_collection.Intersection_Between_Points(grid.Center(cell),grid.Face(axis,face_index),&objects);}

    bool Face_Velocity(const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cells,number_of_cells,collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    COLLISION_GEOMETRY_ID body_id;int triangle_id;TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,grid.Face(axis,face_index),body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];return true;}
    return false;}

    bool Face_Velocity(const int side,const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cells,number_of_cells,collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    COLLISION_GEOMETRY_ID body_id;int triangle_id;TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,grid.X(side==1?face_index-TV_INT::Axis_Vector(axis):face_index),body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];return true;}
    return false;}

    bool Latest_Cell_Crossover(const TV_INT& cell_index,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Center(cell_index);
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    bool Latest_Cell_Crossover_And_Velocity(const TV_INT& cell_index,const T dt,TV& velocity) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Center(cell_index);
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point);return true;}
    return false;}

    bool Latest_Face_Crossover(const TV_INT& face_index,const int axis,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Axis_X_Face(face_index,axis);
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    bool Latest_Velocity_Crossover(const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Face(axis,face_index);
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[axis];return true;}
    return false;}

    bool Latest_Velocity_Crossover(const int side,const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Face(axis,face_index);
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[axis];return true;}
    ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(face_index,face_index-TV_INT::Axis_Vector(axis),collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    int triangle_id;
    if(collision_geometry_collection.Intersection_Between_Points(grid.Center((side==2)?face_index:(face_index-TV_INT::Axis_Vector(axis))),X,body_id,triangle_id,initial_hit_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,initial_hit_point)[axis];
        return true;}
    return false;}

    bool Any_Crossover(const TV& start_X,const TV& end_X,const T dt) const
    {return Any_Simplex_Crossover(start_X,end_X,dt);} // TODO: use object ids

    TV Object_Velocity(const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,const TV& X) const
    {return collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,X);}

    bool Cell_Center_Intersection(const TV_INT& cell_index,const TV_INT* cell_indices_for_body_id,const int number_of_cells_for_body_id,const TV& X,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id,
        TV& intersection_point) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell_indices_for_body_id,number_of_cells_for_body_id,collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    return collision_geometry_collection.Intersection_Between_Points(grid.Center(cell_index),X,body_id,aggregate_id,intersection_point,&objects);}

    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,int ghost_cells,T thickness) const
    {Compute_Simplices_In_Cell(simplices_in_cell,collision_geometry_collection.bodies,ghost_cells,thickness,false);}

//#####################################################################
    void Initialize_Grids();
    virtual void Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor);
    void Compute_Grid_Visibility();
    void Compute_Psi_N(T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR* face_velocities=0) const;
    void Compute_Psi_N_Zero_Velocity(T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR* face_velocities=0) const;
    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
        int ghost_cells,T thickness,bool assume_active) const;
//#####################################################################
};
}
#endif
