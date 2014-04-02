//#####################################################################
// Copyright 2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __GRID_BASED_COLLISION_GEOMETRY_RLE__
#define __GRID_BASED_COLLISION_GEOMETRY_RLE__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/RLE_LEVELSET_ON_A_RAY.h>
namespace PhysBAM{

template <class T_GRID>
class GRID_BASED_COLLISION_GEOMETRY_RLE:public GRID_BASED_COLLISION_GEOMETRY<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::RUN T_RUN;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::HORIZONTAL_GRID::NODE_ITERATOR HORIZONTAL_NODE_ITERATOR;
public:

    typedef GRID_BASED_COLLISION_GEOMETRY<T_GRID> BASE;
    using BASE::collision_geometry_collection;using BASE::grid;using BASE::objects_in_cell;using BASE::cell_neighbors_visible;using BASE::face_neighbors_visible;
    using BASE::occupied_blocks;using BASE::swept_occupied_blocks;

    GRID_BASED_COLLISION_GEOMETRY_RLE(T_GRID& grid_input);
    ~GRID_BASED_COLLISION_GEOMETRY_RLE();

    bool Occupied_Cell_Center(const CELL_ITERATOR& cell)
    {T_BLOCK block(grid,cell.I());return block?Occupied_Block(block):0;}

    bool Swept_Occupied_Cell_Center(const CELL_ITERATOR& cell) const
    {T_BLOCK block(grid,cell.I());return block?Swept_Occupied_Block(block):0;}

    template<class T_FACE>
    bool Occupied_Face_Center(const T_FACE& face)
    {assert(face.Both_Cells_Short());T_BLOCK block(grid,face.cell1.I());return block?Occupied_Block(block):0;}

    template<class T_FACE>
    bool Swept_Occupied_Face_Center(const T_FACE& face) const
    {assert(face.Both_Cells_Short());T_BLOCK block(grid,face.cell1.I());return block?Swept_Occupied_Block(block):0;}

    bool Inside_Any_Simplex_Of_Any_Body(const T_BLOCK& block,const TV& location,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(block,collision_geometry_collection.bodies.m,objects);if(!objects.m) return false;
    return collision_geometry_collection.Inside_Any_Simplex_Of_Any_Body(location,body_id,aggregate_id,&objects);}

    bool Closest_Non_Intersecting_Point_Of_Any_Body(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id) const
    {// since the ray is a general ray, we don't get specific id's.
    return collision_geometry_collection.Closest_Non_Intersecting_Point_Of_Any_Simplicial_Object(ray,body_id);}

/*
    bool Cell_Center_Visible_From_Face(const T_CELL* cell,const int face_index_in_cell) const
    {assert(1<=face_index_in_cell&&face_index_in_cell<=T_GRID::number_of_faces_per_cell);
    ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cell(cell->Cell(),objects);
    return !Intersection_Between_Points(cell->Center(),grid.Face_Location(face_index_in_cell,cell),&objects);}
*/

    bool Face_Velocity(const T_BLOCK& block,const int axis,const TV& face_X,const TV& X,T& face_velocity) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(block,collision_geometry_collection.bodies.m,objects);
    COLLISION_GEOMETRY_ID body_id;int triangle_id;TV intersection_point;
    if(collision_geometry_collection.Intersection_Between_Points(X,face_X,body_id,triangle_id,intersection_point,&objects)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(triangle_id,intersection_point)[axis];return true;}
    return false;}

    bool Latest_Cell_Crossover(const CELL_ITERATOR& cell_iterator,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=cell_iterator.X();
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    template<class T_FACE>
    bool Latest_Velocity_Crossover(const T_FACE& face_iterator,const T dt,T& face_velocity) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=face_iterator.X();
    if(Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point)){
        face_velocity=collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point)[face_iterator.Axis()];return true;}
    return false;}

    bool Any_Crossover(const TV& start_X,const TV& end_X,const T dt) const
    {return Any_Simplex_Crossover(start_X,end_X,dt);} // TODO: use object ids

    TV Object_Velocity(const T_BLOCK& block,const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,const TV& X) const
    {return collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,X);}

    bool Cell_Center_Intersection(const T_BLOCK& block,const int cell_in_block,const TV& X,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id,TV& intersection_point) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(block,collision_geometry_collection.bodies.m,objects);
    return collision_geometry_collection.Intersection_Between_Points(block.Cell_X(cell_in_block),X,body_id,aggregate_id,intersection_point,&objects);}

    bool Refine_Due_To_Objects(const int cell1,const int cell2) const
    {ARRAY<COLLISION_GEOMETRY_ID> objects;objects_in_cell.Get_Objects_For_Cells(cell1,cell2,collision_geometry_collection.bodies.m,objects);
    for(int i=1;i<=objects.m;i++)if(collision_geometry_collection.bodies(objects(i))->refine_nearby_fluid) return true;return false;}

//#####################################################################
    void Initialize_Grids();
    void Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor);
    void Compute_Grid_Visibility();
    void Compute_Psi_N(ARRAY<bool>& psi_N,ARRAY<T>* face_velocities) const;
private:
    struct Rasterize_Objects_Accurate{template<class T_FACE> static void Apply(GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list);};
    struct Compute_Cell_Neighbors_Visible_Helper{template<class T_FACE> static void Apply(GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list);};
    struct Compute_Face_Neighbors_Visible_Helper{template<class T_FACE> static void Apply(GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list,ARRAY<int>& face_cell);};
    struct Compute_Psi_N_Helper{template<class T_FACE> static void Apply(const GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list,ARRAY<bool>& psi_N,ARRAY<T>* face_velocities);};
//#####################################################################
};
}
#endif
#endif
