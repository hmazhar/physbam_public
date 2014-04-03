//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Rasterization/RIGID_GEOMETRY_RASTERIZATION_HELPER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
namespace RASTERIZATION{
//#####################################################################
// Function Rasterize_Object
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id)
{
    Rasterize_Object_Generic(collision_geometry,grid,objects,id);
}
//#####################################################################
// Function Compute_Occupied_Blocks
//#####################################################################
template<class T,class TV,class T_GRID> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,
    typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor)
{
    if(collision_geometry.Number_Of_Simplices()) Compute_Occupied_Blocks_Generic(collision_geometry,grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Rasterize_Box_Onto_Blocks
//#####################################################################
template<class T_GRID> void Rasterize_Box_Onto_Blocks_Helper(ARRAY<bool>& occupied,const typename T_GRID::CELL& cell,const RANGE<typename T_GRID::VECTOR_T>& box)
{
    if(!box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(cell.Has_Children()) for(int i=0;i<T_GRID::number_of_children_per_cell;i++) Rasterize_Box_Onto_Blocks_Helper<T_GRID>(occupied,*cell.Child(i),box);
    else occupied(cell.Cell())=true;
}
template<class TV,class T_GRID> void Rasterize_Box_Onto_Blocks(const T_GRID& grid,ARRAY<bool>& occupied,const RANGE<TV>& box)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    typedef typename T_GRID::UNIFORM_GRID T_UNIFORM_GRID;
    TV DX_over_two=(typename TV::SCALAR).5*grid.Minimum_Cell_DX();
    TV_INT min_index=grid.uniform_grid.Clamp_To_Cell(box.Minimum_Corner()-DX_over_two,grid.number_of_ghost_cells);
    TV_INT max_index=grid.uniform_grid.Clamp_To_Cell(box.Maximum_Corner()-DX_over_two,grid.number_of_ghost_cells);
    for(typename T_UNIFORM_GRID::CELL_ITERATOR iterator(grid.uniform_grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
        Rasterize_Box_Onto_Blocks_Helper<T_GRID>(occupied,*grid.cells(iterator.Cell_Index()),box);
}
//#####################################################################
// Function Rasterize_Box
//#####################################################################
template<class T_GRID> static void Rasterize_Box_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const typename T_GRID::CELL& cell,
    const RANGE<typename T_GRID::VECTOR_T>& box,const COLLISION_GEOMETRY_ID id)
{
    if(!box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(cell.Has_Children()) for(int i=0;i<T_GRID::number_of_children_per_cell;i++) Rasterize_Box_Helper<T_GRID>(objects_in_cell,*cell.Child(i),box,id);
    else objects_in_cell.Add_Object_To_Cell(cell.Cell(),id);
}
template<class TV,class T_GRID> void Rasterize_Box(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    typedef typename T_GRID::UNIFORM_GRID T_UNIFORM_GRID;
    TV_INT min_index=grid.uniform_grid.Clamp_To_Cell(box.Minimum_Corner(),grid.number_of_ghost_cells);
    TV_INT max_index=grid.uniform_grid.Clamp_To_Cell(box.Maximum_Corner(),grid.number_of_ghost_cells);
    for(typename T_UNIFORM_GRID::CELL_ITERATOR iterator(grid.uniform_grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
        Rasterize_Box_Helper<T_GRID>(objects_in_cell,*grid.cells(iterator.Cell_Index()),box,id);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void Rasterize_Object(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID&,OBJECTS_IN_CELL<DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID, \
        COLLISION_GEOMETRY_ID>&,const COLLISION_GEOMETRY_ID&); \
    template void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID&, \
        GRID_ARRAYS_POLICY<DYADIC_GRID_POLICY<VECTOR<bool,d> >::DYADIC_GRID>::ARRAYS_SCALAR&,const bool,const T&,const T&);
#define INSTANTIATION_HELPER_2(T,d) \
    INSTANTIATION_HELPER(T,d); \
    template void Rasterize_Box_Onto_Blocks(const DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID&,ARRAY<bool>&,const RANGE<VECTOR<T,d> >&); \
    template void Rasterize_Box(const DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID&,OBJECTS_IN_CELL<DYADIC_GRID_POLICY<VECTOR<T,d> >::DYADIC_GRID,COLLISION_GEOMETRY_ID>&,const RANGE<VECTOR<T,d> >&,const COLLISION_GEOMETRY_ID);

INSTANTIATION_HELPER_2(float,1);
// INSTANTIATION_HELPER_2(float,2);
// INSTANTIATION_HELPER_2(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER_2(double,1);
// INSTANTIATION_HELPER_2(double,2);
// INSTANTIATION_HELPER_2(double,3);
#endif
};
};
#endif
