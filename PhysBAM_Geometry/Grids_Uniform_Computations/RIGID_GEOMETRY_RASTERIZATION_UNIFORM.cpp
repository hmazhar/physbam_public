//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Rasterization/RIGID_GEOMETRY_RASTERIZATION_HELPER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV> class GRID;

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
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    if(collision_geometry.Number_Of_Simplices()) Compute_Occupied_Blocks_Generic(collision_geometry,grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
    else{
        for(typename T_GRID::NODE_ITERATOR node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()){
            TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<T_GRID> block(grid,block_index);
            for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++){TV_INT cell=block.Cell(cell_index);
                if(collision_geometry.Implicit_Geometry_Extended_Value(grid.X(cell))<=extra_thickness){occupied(block_index)=true;break;}}}}
}
//#####################################################################
// Function Rasterize_Box_Onto_Blocks
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box_Onto_Blocks(const T_GRID& grid,typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE& occupied,const RANGE<TV>& box)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    TV DX_over_two=(typename TV::SCALAR).5*grid.dX;
    TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner()-DX_over_two,3);
    TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner()-DX_over_two,3);
    for(typename T_GRID::CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()) occupied(iterator.Cell_Index()+TV_INT::All_Ones_Vector())=true;
}
//#####################################################################
// Function Rasterize_Box
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),3);
    TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),3);
    for(typename T_GRID::CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()) objects_in_cell.Add_Object_To_Cell(iterator.Cell_Index(),id);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void Rasterize_Object(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const GRID<VECTOR<T,d> >&,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>&,const COLLISION_GEOMETRY_ID&); \
    template void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const GRID<VECTOR<T,d> >&,GRID_ARRAYS_POLICY<GRID<VECTOR<bool,d> > >::ARRAYS_SCALAR&,const bool, \
        const T&,const T&); \
    template void Rasterize_Box_Onto_Blocks(const GRID<VECTOR<T,d> >&,GRID_ARRAYS_POLICY<GRID<VECTOR<bool,d> > >::ARRAYS_SCALAR&,const RANGE<VECTOR<T,d> >&); \
    template void Rasterize_Box(const GRID<VECTOR<T,d> >&,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>&,const RANGE<VECTOR<T,d> >&,const COLLISION_GEOMETRY_ID);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
};
};
