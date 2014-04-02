//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_1D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Rasterization/RIGID_GEOMETRY_RASTERIZATION_HELPER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
namespace RASTERIZATION{
//#####################################################################
// Function Specializations for grids which we don't support
//#####################################################################
SPECIALIZE_RASTERIZE_OBJECT(RLE_GRID_1D<float>)
SPECIALIZE_COMPUTE_OCCUPIED(RLE_GRID_1D<float>)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
SPECIALIZE_RASTERIZE_OBJECT(RLE_GRID_1D<double>)
SPECIALIZE_COMPUTE_OCCUPIED(RLE_GRID_1D<double>)
#endif
//#####################################################################
// Function Rasterize_Object
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id)
{
    Rasterize_Object_Generic(collision_geometry,grid,objects,id);
}
template<> void Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<VECTOR<float,1> >&,const RLE_GRID_1D<float>&,ARRAY<bool>&,const bool,const float,const float){PHYSBAM_NOT_IMPLEMENTED();}
template<> void Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<VECTOR<double,1> >&,const RLE_GRID_1D<double>&,ARRAY<bool>&,const bool,const double,const double){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Compute_Occupied_Blocks
//#####################################################################
template<class T,class TV,class T_GRID> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const T_GRID& grid,
    typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor)
{
    if(collision_geometry.Number_Of_Simplices()) Compute_Occupied_Blocks_Generic(collision_geometry,grid,occupied,with_body_motion,extra_thickness,body_thickness_factor);
    else PHYSBAM_FATAL_ERROR();
/*        for(typename T_GRID::BLOCK_ITERATOR block(grid,grid.number_of_ghost_cells-1);block;block++){int b=block.Block();
            for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
            if(collision_geometry.Implicit_Geometry_Extended_Value(grid.X(block.Cell(i)))<=extra_thickness){occupied(b)=true;break;}}}*/
}
//#####################################################################
// Function Rasterize_Box_Onto_Blocks
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box_Onto_Blocks(const T_GRID& grid,ARRAY<bool>& occupied,const RANGE<TV>& box)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    GRID<TV> block_grid=grid.uniform_grid.Get_MAC_Grid(); // TODO: make Block_Index work for regular grids
    TV_INT min_index=block_grid.Block_Index(box.Minimum_Corner(),grid.number_of_ghost_cells-1)-TV_INT::All_Ones_Vector(); // rle indexes block columns by first cell
    TV_INT max_index=block_grid.Block_Index(box.Maximum_Corner(),grid.number_of_ghost_cells-1)-TV_INT::All_Ones_Vector();
    typename T_GRID::BOX_HORIZONTAL_INT horizontal_box=RANGE<TV_INT>(min_index,max_index).Get_Horizontal_Box();
    for(typename T_GRID::HORIZONTAL_GRID::NODE_ITERATOR iterator(grid.horizontal_grid,horizontal_box);iterator.Valid();iterator.Next()){
        RANGE<VECTOR<int,1> > block_range=grid.Block_Range(iterator.Node_Index(),min_index.y,max_index.y);
        for(int b=block_range.min_corner.x;b<=block_range.max_corner.x;b++)occupied(b)=true;}
}
//#####################################################################
// Function Rasterize_Box
//#####################################################################
template<class TV,class T_GRID> void Rasterize_Box(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id)
{
    typedef typename REBIND<TV,int>::TYPE TV_INT;
    typedef typename T_GRID::RUN T_RUN;
    TV_INT min_index=grid.uniform_grid.Clamp_To_Cell(box.Minimum_Corner(),grid.number_of_ghost_cells);
    TV_INT max_index=grid.uniform_grid.Clamp_To_Cell(box.Maximum_Corner(),grid.number_of_ghost_cells);
    typename T_GRID::BOX_HORIZONTAL_INT horizontal_box=RANGE<TV_INT>(min_index,max_index).Get_Horizontal_Box();
    for(typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR iterator(grid.horizontal_grid,horizontal_box);iterator.Valid();iterator.Next()){
        const ARRAY<T_RUN>& column=grid.columns(iterator.Cell_Index());
        const T_RUN* start_run=grid.Clamped_Run_In_Column(column,min_index.y);int start_cell=start_run->is_long?start_run->cell+2:start_run->cell+min_index.y-start_run->jmin;
        const T_RUN* end_run=grid.Clamped_Run_In_Column(column,max_index.y);int end_cell=end_run->is_long?end_run->cell-1:end_run->cell+max_index.y-end_run->jmin;
        for(int c=start_cell;c<=end_cell;c++)objects_in_cell.Add_Object_To_Cell(c,id);}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void Rasterize_Object(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID&,OBJECTS_IN_CELL<RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID, \
        COLLISION_GEOMETRY_ID>&,const COLLISION_GEOMETRY_ID&);
#define INSTANTIATION_HELPER_2(T,d) \
    INSTANTIATION_HELPER(T,d); \
    template void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<VECTOR<T,d> >&,const RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID&, \
        GRID_ARRAYS_POLICY<RLE_GRID_POLICY<VECTOR<bool,d> >::RLE_GRID>::ARRAYS_SCALAR&,const bool,const T&,const T&); \
    template void Rasterize_Box_Onto_Blocks(const RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID&,ARRAY<bool>&,const RANGE<VECTOR<T,d> >&); \
    template void Rasterize_Box(const RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID&,OBJECTS_IN_CELL<RLE_GRID_POLICY<VECTOR<T,d> >::RLE_GRID,COLLISION_GEOMETRY_ID>&,const RANGE<VECTOR<T,d> >&,const COLLISION_GEOMETRY_ID);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER_2(float,2);
INSTANTIATION_HELPER_2(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER_2(double,2);
INSTANTIATION_HELPER_2(double,3);
#endif
};
};
#endif
