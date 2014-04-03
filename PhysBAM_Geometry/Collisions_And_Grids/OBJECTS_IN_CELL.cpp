//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
using namespace PhysBAM;
template <class T_GRID,class ID> OBJECTS_IN_CELL<T_GRID,ID>::
OBJECTS_IN_CELL()
{
}
template <class T_GRID,class ID> OBJECTS_IN_CELL<T_GRID,ID>::
~OBJECTS_IN_CELL() 
{
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Reset(const T_GRID& grid,const int number_of_ghost_cells)
{
    object_in_cell.Resize(grid.Cell_Indices(number_of_ghost_cells),false,false);
    object_in_cell.Fill(ID());
    object_list.Clean_Memory();
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Add_Object_To_Cell(const T_INDEX& cell_index,const ID object_id)
{
    assert(object_id>ID());
    if(!object_in_cell(cell_index)) object_in_cell(cell_index)=object_id;
    else if(object_in_cell(cell_index)<ID()) object_list(ID(-Value(object_in_cell(cell_index)))).Append_Unique(object_id);
    else if(object_in_cell(cell_index)!=object_id){
        object_list.Resize(object_list.m+1);object_list.Last().Resize(2);
        object_list.Last()(1)=object_in_cell(cell_index);
        object_list.Last()(2)=object_id;
        object_in_cell(cell_index)=ID(-Value(object_list.m));}
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells_Cell(const T_INDEX& cell) const
{
    if(object_in_cell(cell)>ID()){
        ID object=object_in_cell(cell);
        if(!operation_hash.Is_Marked_Current(object)){operation_hash.Mark(object);merge.Append(object);}}
    else if(object_in_cell(cell)<ID())for(int item=1;item<=object_list(ID(-Value(object_in_cell(cell)))).m;item++){
        ID object=object_list(ID(-Value(object_in_cell(cell))))(item);
        if(!operation_hash.Is_Marked_Current(object)){operation_hash.Mark(object);merge.Append(object);}}
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_INDEX* cells,const int number_of_cells,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    for(int i=0;i<number_of_cells;i++) Get_Objects_For_Cells_Cell(cells[i]);
    Get_Objects_For_Cells_End(objects);
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_BLOCK& block,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    for(int i=1;i<T_GRID::number_of_cells_per_block;i++) Get_Objects_For_Cells_Cell(block.Cell(i));
    Get_Objects_For_Cells_End(objects);
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cell(const T_INDEX& cell_index,ARRAY<ID>& objects) const
{
    assert(!objects.m);
    if(object_in_cell(cell_index)>ID()) objects.Append(object_in_cell(cell_index));
    else if(object_in_cell(cell_index)<ID()) objects=object_list(ID(-Value(object_in_cell(cell_index))));
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_INDEX& cell_index1,const T_INDEX& cell_index2,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    Get_Objects_For_Cells_Cell(cell_index1);Get_Objects_For_Cells_Cell(cell_index2);
    Get_Objects_For_Cells_End(objects);
}
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells_Start(const ID number_of_collision_bodies) const
{
    if(operation_hash.operations.m!=number_of_collision_bodies)operation_hash.Initialize(number_of_collision_bodies);
    merge.Remove_All();operation_hash.Next_Operation();
}
template class OBJECTS_IN_CELL<GRID<VECTOR<float,1> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<float,2> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<float,3> >,COLLISION_GEOMETRY_ID>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OBJECTS_IN_CELL<GRID<VECTOR<double,1> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<double,2> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<double,3> >,COLLISION_GEOMETRY_ID>;
#endif

#if COMPILE_WITH_BINTREE_SUPPORT
template class OBJECTS_IN_CELL<BINTREE_GRID<float>,COLLISION_GEOMETRY_ID>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OBJECTS_IN_CELL<BINTREE_GRID<double>,COLLISION_GEOMETRY_ID>;
#endif
#endif
