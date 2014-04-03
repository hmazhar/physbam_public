//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __GRID_BASED_COLLISION_GEOMETRY_DYADIC__
#define __GRID_BASED_COLLISION_GEOMETRY_DYADIC__

#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
namespace PhysBAM{

template <class T_GRID>
class GRID_BASED_COLLISION_GEOMETRY_DYADIC:public GRID_BASED_COLLISION_GEOMETRY<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::INDEX INDEX;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT T_SEGMENT;
    typedef typename T_GRID::CELL T_CELL;typedef typename T_GRID::UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    typedef GRID_BASED_COLLISION_GEOMETRY<T_GRID> BASE;
    using BASE::collision_geometry_collection;using BASE::grid;using BASE::objects_in_cell;using BASE::cell_neighbors_visible;using BASE::face_neighbors_visible;
    using BASE::Get_Body_Penetration;using BASE::occupied_blocks;using BASE::swept_occupied_blocks;using BASE::Is_Active;

    const T_ARRAYS_BOOL* outside_fluid;
    bool use_collision_face_neighbors;

    GRID_BASED_COLLISION_GEOMETRY_DYADIC(T_GRID& grid_input);
    ~GRID_BASED_COLLISION_GEOMETRY_DYADIC();

    bool Occupied_Cell_Center(const int cell_index)
    {return occupied_blocks(cell_index);} // this will check the occupied block that is right,top,behind the current cell.

    bool Swept_Occupied_Cell_Center(const int cell_index) const
    {return swept_occupied_blocks(cell_index);} // this will check the occupied block that is right,top,behind of the current cell.

    bool Occupied_Face_Center(FACE_ITERATOR& iterator)
    {return occupied_blocks(iterator.First_Cell_Index());}

    bool Swept_Occupied_Face_Center(FACE_ITERATOR& iterator) const
    {return swept_occupied_blocks(iterator.First_Cell_Index());}

    bool Inside_Any_Simplex_Of_Any_Body(const BLOCK_DYADIC<T_GRID>& block,const TV& location,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id) const
    {int cell_indices[T_GRID::number_of_cells_per_block];block.All_Cell_Indices(cell_indices);
    ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell_indices,T_GRID::number_of_cells_per_block,collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    return collision_geometry_collection.Inside_Any_Simplex_Of_Any_Body(location,body_id,aggregate_id,&objects);}

    bool Closest_Non_Intersecting_Point_Of_Any_Body(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id) const
    {// since the ray is a general ray, we don't get specific id's.
    return collision_geometry_collection.Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(ray,body_id);}

    bool Cell_Center_Visible_From_Face(const T_CELL* cell,const int face_index_in_cell) const
    {assert(1<=face_index_in_cell&&face_index_in_cell<=T_GRID::number_of_faces_per_cell);
    ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cell(cell->Cell(),objects);
    return !collision_geometry_collection.Intersection_Between_Points(cell->Center(),grid.Face_Location(face_index_in_cell,cell),&objects);}

    bool Face_Velocity(const int axis,const int& face_index,const int* cell_indices_for_body_id,const int number_of_cells_for_body_id,const TV& X,T& face_velocity) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell_indices_for_body_id,number_of_cells_for_body_id,collision_geometry_collection.bodies.m,objects);
    COLLISION_GEOMETRY_ID body_id;int triangle_id;TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,grid.Face_Location(face_index),body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];return true;}
    return false;}

    bool Latest_Cell_Crossover(const int cell_index,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Cell_Location(cell_index);
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    bool Latest_Velocity_Crossover(const int face_index,const T dt,T& face_velocity) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Face_Location(face_index);
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[grid.Face_Axis(face_index)+1];return true;}
    return false;}

    bool Any_Crossover(const TV& start_X,const TV& end_X,const T dt) const
    {return Any_Simplex_Crossover(start_X,end_X,dt);} // TODO: use object ids

    TV Object_Velocity(const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,const TV& X) const
    {return collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,X);}

    bool Cell_Center_Intersection(const int& cell_index,const int* cell_indices_for_body_id,const int number_of_cells_for_body_id,const TV& X,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id,
        TV& intersection_point) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell_indices_for_body_id,number_of_cells_for_body_id,collision_geometry_collection.bodies.m,objects);
    return collision_geometry_collection.Intersection_Between_Points(grid.Cell_Location(cell_index),X,body_id,aggregate_id,intersection_point,&objects);}

    void Transfer_Occupied_Cell_Values_To_Parents()
    {PHYSBAM_NOT_IMPLEMENTED();} // TODO: currently just a place holder
    
    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,INDEX>& simplices_in_cell,int ghost_cells,T thickness) const
    {Compute_Simplices_In_Cell(simplices_in_cell,collision_geometry_collection.bodies,ghost_cells,thickness,false);}

//#####################################################################
    void Initialize_Grids();
    void Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor);
    void Compute_Grid_Visibility();
    void Compute_Psi_N(ARRAY<bool>& psi_N,ARRAY<T>* face_velocities) const;
    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,INDEX>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
        int ghost_cells,T thickness,bool assume_active) const;
//#####################################################################
};
}
#endif
#endif
