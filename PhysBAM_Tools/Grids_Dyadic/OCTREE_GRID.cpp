//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_GRID
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
template<class T> const char* OCTREE_GRID<T>::name="octree";
//#####################################################################
// Constructor
//#####################################################################
template<class T> OCTREE_GRID<T>::
OCTREE_GRID()
    :uniform_grid(2,2,2,0,1,0,1,0,1),minimum_cell_size(0),number_of_ghost_cells(0),number_of_cells(0),number_of_nodes(0),number_of_faces(0),maximum_depth(0)
{
    Tree_Topology_Changed();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OCTREE_GRID<T>::
~OCTREE_GRID()
{
    for(int i=1;i<=allocated_children.m;i++)delete allocated_children(i);allocated_children.Resize(0);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Initialize(const GRID<TV> uniform_grid_input,const int maximum_depth_input,const int number_of_ghost_cells_input,const bool use_nodes,const bool use_faces)
{
    int i,j,ij;uniform_grid=uniform_grid_input;Set_Maximum_Depth(maximum_depth_input);number_of_ghost_cells=number_of_ghost_cells_input;
    if(uniform_grid.numbers_of_cells.x%2 || uniform_grid.numbers_of_cells.y%2 || uniform_grid.numbers_of_cells.z%2){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Only even sized grids are supported (grid = "<<uniform_grid.numbers_of_cells<<")"<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
    Tree_Topology_Changed();
    number_of_cells=number_of_nodes=number_of_faces=0;
    for(i=1;i<=allocated_children.m;i++)delete allocated_children(i);allocated_children.Resize(0);
    cells.Resize(1-number_of_ghost_cells,uniform_grid.numbers_of_cells.x+number_of_ghost_cells,1-number_of_ghost_cells,uniform_grid.numbers_of_cells.y+number_of_ghost_cells,1-number_of_ghost_cells,uniform_grid.numbers_of_cells.z+number_of_ghost_cells);
    ARRAY<int,TV_INT> nodes(cells.domain.min_corner.x,cells.domain.max_corner.x+1,cells.domain.min_corner.y,cells.domain.max_corner.y+1,cells.domain.min_corner.z,cells.domain.max_corner.z+1);
    ARRAY<int,TV_INT> faces_u(cells.domain.min_corner.x,cells.domain.max_corner.x+1,cells.domain.min_corner.y,cells.domain.max_corner.y,cells.domain.min_corner.z,cells.domain.max_corner.z),
        faces_v(cells.domain.min_corner.x,cells.domain.max_corner.x,cells.domain.min_corner.y,cells.domain.max_corner.y+1,cells.domain.min_corner.z,cells.domain.max_corner.z),
        faces_w(cells.domain.min_corner.x,cells.domain.max_corner.x,cells.domain.min_corner.y,cells.domain.max_corner.y,cells.domain.min_corner.z,cells.domain.max_corner.z+1);
    if(use_nodes) for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y+1;j++)for(ij=cells.domain.min_corner.z;ij<=cells.domain.max_corner.z+1;ij++)nodes(i,j,ij)=++number_of_nodes;
    if(use_faces){
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y;j++)for(ij=cells.domain.min_corner.z;ij<=cells.domain.max_corner.z;ij++)faces_u(i,j,ij)=++number_of_faces;
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y+1;j++)for(ij=cells.domain.min_corner.z;ij<=cells.domain.max_corner.z;ij++)faces_v(i,j,ij)=++number_of_faces;
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y;j++)for(ij=cells.domain.min_corner.z;ij<=cells.domain.max_corner.z+1;ij++)faces_w(i,j,ij)=++number_of_faces;}
    ARRAY<int,VECTOR<int,1> > nodes_to_input(0,26);ARRAY<int,VECTOR<int,1> > faces_to_input(0,35);
    for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i+=2)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y;j+=2)for(ij=cells.domain.min_corner.z;ij<=cells.domain.max_corner.z;ij+=2){
        OCTREE_CHILDREN<T>* owner=new OCTREE_CHILDREN<T>(use_nodes,use_faces);
        allocated_children.Append(owner);
        TV center=uniform_grid.Node(i+1,j+1,ij+1),DX=uniform_grid.dX;
        if(use_nodes) for(int ii=i;ii<=i+2;ii++)for(int jj=j;jj<=j+2;jj++)for(int iijj=ij;iijj<=ij+2;iijj++)
            nodes_to_input((ii-i)+(jj-j)*3+(iijj-ij)*9)=nodes(ii,jj,iijj); // CHECK THIS
        if(use_faces){
            for(int ii=i;ii<=i+2;ii++)for(int jj=j;jj<=j+1;jj++)for(int iijj=ij;iijj<=ij+1;iijj++)
            {faces_to_input((ii-i)+(jj-j)*3+(iijj-ij)*6)=faces_u(ii,jj,iijj);}
            for(int ii=i;ii<=i+1;ii++)for(int jj=j;jj<=j+2;jj++)for(int iijj=ij;iijj<=ij+1;iijj++)
            {faces_to_input(12+(ii-i)*3+(jj-j)+(iijj-ij)*6)=faces_v(ii,jj,iijj);}
            for(int ii=i;ii<=i+1;ii++)for(int jj=j;jj<=j+1;jj++)for(int iijj=ij;iijj<=ij+2;iijj++)
            {faces_to_input(24+(ii-i)*3+(jj-j)*6+(iijj-ij))=faces_w(ii,jj,iijj);}}
        owner->Initialize_Pseudo_Root_Cells(number_of_cells,nodes_to_input,faces_to_input,center,DX);
        cells(i,j,ij)=&owner->children[0];cells(i+1,j,ij)=&owner->children[1];cells(i,j+1,ij)=&owner->children[2];cells(i+1,j+1,ij)=&owner->children[3];
        cells(i,j,ij+1)=&owner->children[4];cells(i+1,j,ij+1)=&owner->children[5];cells(i,j+1,ij+1)=&owner->children[6];cells(i+1,j+1,ij+1)=&owner->children[7];}
}
//#####################################################################
// Function Tree_Topology_Changed
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Tree_Topology_Changed()
{
    cell_pointer_from_index.Resize(0);cell_pointer_from_index_up_to_date=false;
    neighbors.Resize(0);neighbors_up_to_date=false;
    node_neighbors.Resize(0);node_neighbors_up_to_date=false;
    fully_refined_block.Resize(0);fully_refined_block_up_to_date=false;
    node_iterator_data_up_to_date=false;face_iterator_data_up_to_date=false;
}
//#####################################################################
// Function Leaf_Cell
//#####################################################################
template<class T> OCTREE_CELL<T>* OCTREE_GRID<T>::
Leaf_Cell(const TV& location,const T thickness) const
{
    TV_INT i=uniform_grid.Cell(location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    if(cells.Domain_Indices().Lazy_Outside(i)) return 0;
    return Inside_Offspring(cells(i),location,thickness);
}
//#####################################################################
// Function Clamped_Leaf_Cell
//#####################################################################
template<class T> OCTREE_CELL<T>* OCTREE_GRID<T>::
Clamped_Leaf_Cell(const TV& location,const T thickness) const
{
    TV clamped_location=uniform_grid.Clamp(location,number_of_ghost_cells);
    TV_INT i=uniform_grid.Cell(clamped_location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    i=TV_INT::Componentwise_Min(i,TV_INT(cells.domain.max_corner.x,cells.domain.max_corner.y,cells.domain.max_corner.z)); // need this clamp because Cell could potentially return i, j or ij beyond last index
    return Inside_Offspring(cells(i),clamped_location,thickness);
}
//#####################################################################
// Function Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> const OCTREE_CELL<T>* OCTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(const OCTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<OCTREE_CELL<T>*,6> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell,tolerance_times_small_cell,tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else if(location.y<center.y-DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(3);
        else if(location.y>center.y+DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(4);
        else if(location.z<center.z-DX_over_two.z) starting_cell=neighbors(starting_cell->Cell())(5);
        else if(location.z>center.z+DX_over_two.z) starting_cell=neighbors(starting_cell->Cell())(6);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Function Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> OCTREE_CELL<T>* OCTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(OCTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<OCTREE_CELL<T>*,6> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell,tolerance_times_small_cell,tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else if(location.y<center.y-DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(3);
        else if(location.y>center.y+DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(4);
        else if(location.z<center.z-DX_over_two.z) starting_cell=neighbors(starting_cell->Cell())(5);
        else if(location.z>center.z+DX_over_two.z) starting_cell=neighbors(starting_cell->Cell())(6);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Function Inside_Cell
//#####################################################################
template<class T> bool OCTREE_GRID<T>::
Inside_Cell(const OCTREE_CELL<T>* cell,const TV& location) const
{
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    if(location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x && location.y>=center.y-DX_over_two.y && location.y<=center.y+DX_over_two.y && 
        location.z>=center.z-DX_over_two.z && location.z<=center.z+DX_over_two.z) return true;
    else return false;
}
//#####################################################################
// Function Inside_Thickened_Cell
//#####################################################################
template<class T> bool OCTREE_GRID<T>::
Inside_Thickened_Cell(const OCTREE_CELL<T>* cell,const TV& location,const T thickness) const
{
    TV center(cell->Center()),DX_over_two(((T).5+thickness)*cell->DX());
    if(location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x && location.y>=center.y-DX_over_two.y && 
        location.y<=center.y+DX_over_two.y && location.z>=center.z-DX_over_two.z && location.z<=center.z+DX_over_two.z) return true;
    else return false;
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> OCTREE_CELL<T>* OCTREE_GRID<T>::
Inside_Offspring(OCTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        int child_to_use=0;if(location.x>center.x)child_to_use+=1;if(location.y>center.y)child_to_use+=2;if(location.z>center.z)child_to_use+=4;
        cell=cell->Child(child_to_use);}
    return cell;
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> const OCTREE_CELL<T>* OCTREE_GRID<T>::
Inside_Offspring(const OCTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        int child_to_use=0;if(location.x>center.x)child_to_use+=1;if(location.y>center.y)child_to_use+=2;if(location.z>center.z)child_to_use+=4;
        cell=cell->Child(child_to_use);}
    return cell;
}
//#####################################################################
// Function Refine_Cell
//#####################################################################
// refines the given cell so that the node at the given location is touching a refined cell instead of a coarse one
// the location needs to be inside or on one of the faces of the cell
template<class T> void OCTREE_GRID<T>::
Refine_Cell(const int max_depth,OCTREE_CELL<T>* cell,const TV& location,ARRAY<OCTREE_CELL<T>*>* new_cells,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes, 
            ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const OCTREE_GRID<T>* grid)
{
    int depth_of_cell=cell->Depth_Of_This_Cell();assert(depth_of_cell<=maximum_depth);
    if(depth_of_cell>=max_depth) return;
    if(!Inside_Thickened_Cell(cell,location)) return; // make sure that the node is actually in the cell
    if(!cell->Has_Children()) cell->Create_Children(number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,grid);
    for(int i=0;i<number_of_children_per_cell;i++) Refine_Cell(max_depth,cell->Child(i),location, new_cells,new_nodes,new_faces,grid); // recur on the children
}
//#####################################################################
// Function Compact_Array_Indices
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array,ARRAY<int>* face_mapping_array)
{
    if(cell_mapping_array){cell_mapping_array->Resize(number_of_cells,false);ARRAYS_COMPUTATIONS::Fill(*cell_mapping_array,0);number_of_cells=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j++)
                for(int ij=1-number_of_ghost_cells;ij<=uniform_grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
                    cells(i,j,ij)->Create_Cell_Compaction_Mapping(*cell_mapping_array,number_of_cells);}
    if(face_mapping_array){face_mapping_array->Resize(number_of_faces,false);ARRAYS_COMPUTATIONS::Fill(*face_mapping_array,0);number_of_faces=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j+=2)
                for(int ij=1-number_of_ghost_cells;ij<=uniform_grid.numbers_of_cells.z+number_of_ghost_cells;ij+=2)
                    cells(i,j,ij)->owner->Create_Face_Compaction_Mapping_Helper(*face_mapping_array,number_of_faces);}
    if(node_mapping_array){node_mapping_array->Resize(number_of_nodes,false);ARRAYS_COMPUTATIONS::Fill(*node_mapping_array,0);number_of_nodes=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j+=2)
                for(int ij=1-number_of_ghost_cells;ij<=uniform_grid.numbers_of_cells.z+number_of_ghost_cells;ij+=2)
                    cells(i,j,ij)->owner->Create_Node_Compaction_Mapping_Helper(*node_mapping_array,number_of_nodes);}
}
//#####################################################################
//#####################################################################
// Function Refine_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Refine_Cells_Intersecting_Box_Helper(OCTREE_GRID<T>& grid,OCTREE_CELL<T>& cell,const RANGE<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    if(cell.Depth_Of_This_Cell()>=refinement_depth || !box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){ // refine
        refined_cells.Append(&cell);cell.Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}
    for(int i=0;i<8;i++) Refine_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,refined_cells,refinement_depth);
}
template<class T> void OCTREE_GRID<T>::
Refine_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<OCTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,OCTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,OCTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index+=TV_INT(1,1,1);
    int depth=maximum_depth;if(refinement_depth) depth=refinement_depth;
    for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) for(int ij=min_index.z;ij<=max_index.z;ij++) Refine_Cells_Intersecting_Box_Helper(*this,*cells(i,j,ij),box,refined_cells,depth);
}
//#####################################################################
// Function Get_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Get_Cells_Intersecting_Box_Helper(OCTREE_GRID<T>& grid,OCTREE_CELL<T>& cell,const RANGE<VECTOR<T,3> >& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells)
{
    if(!box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){
        intersecting_cells.Append(&cell);}
    else for(int i=0;i<8;i++) Get_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,intersecting_cells);
}
template<class T> void OCTREE_GRID<T>::
Get_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells)
{
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,OCTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,OCTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index+=TV_INT(1,1,1);
    for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) for(int ij=min_index.z;ij<=max_index.z;ij++)
        Get_Cells_Intersecting_Box_Helper(*this,*cells(i,j,ij),box,intersecting_cells);
}
//#####################################################################
// Function Node_Iterator_Data
//#####################################################################
namespace PhysBAM{template<class T> struct NODE_ITERATOR_DATA_HELPER{const OCTREE_GRID<T>* grid;int* number_of_nodes;};}
template<class T> static void Node_Iterator_Data_Helper(void* data,const OCTREE_CELL<T>* cell1,const OCTREE_CELL<T>* cell2,const OCTREE_CELL<T>* cell3,const OCTREE_CELL<T>* cell4,
                                                        const OCTREE_CELL<T>* cell5,const OCTREE_CELL<T>* cell6,const OCTREE_CELL<T>* cell7,const OCTREE_CELL<T>* cell8)
{
    NODE_ITERATOR_DATA_HELPER<T>* helper=(NODE_ITERATOR_DATA_HELPER<T>*)data;
    const OCTREE_CELL<T>* cells[8]={cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8};int maximum_depth,deepest_cell;MAP_OCTREE_MESH<T>::Get_Deepest_Cell(8,cells,deepest_cell,maximum_depth);
    int node=7-deepest_cell;int node_index=cells[deepest_cell]->Node(node);(*helper->number_of_nodes)++;
    helper->grid->node_iterator_deepest_cells(node_index)=cells[deepest_cell]->Cell();helper->grid->nodes(node_index)=node;
}
template<class T> void OCTREE_GRID<T>::
Node_Iterator_Data() const
{
    if(!node_iterator_data_up_to_date){
        node_iterator_deepest_cells.Resize(0);node_iterator_deepest_cells.Resize(number_of_nodes);
        nodes.Resize(number_of_nodes);
        
        NODE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;
        
        int number_of_internal_nodes=0;
        helper.number_of_nodes=&number_of_internal_nodes;
        MAP_OCTREE_MESH<T>::Map_Nodes(uniform_grid,cells,0,&helper,Node_Iterator_Data_Helper);
        internal_nodes.Resize(number_of_internal_nodes);
        
        int number_of_boundary_nodes=0;
        helper.number_of_nodes=&number_of_boundary_nodes;
        MAP_OCTREE_MESH<T>::Map_Boundary_Nodes(uniform_grid,cells,&helper,Node_Iterator_Data_Helper);
        boundary_nodes.Resize(number_of_boundary_nodes);
        individual_side_boundary_nodes.Resize(6);
        for(int i=1;i<=6;i++)individual_side_boundary_nodes(i).Resize(number_of_boundary_nodes);
        
        int number_of_ghost_nodes=0;
        helper.number_of_nodes=&number_of_ghost_nodes;
        MAP_OCTREE_MESH<T>::Map_Ghost_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
        MAP_OCTREE_MESH<T>::Map_Exterior_Ghost_Cell_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
        ghost_nodes.Resize(number_of_ghost_nodes);
        individual_side_ghost_nodes.Resize(6);
        for(int i=1;i<=6;i++)individual_side_ghost_nodes(i).Resize(number_of_ghost_nodes);
        
        number_of_internal_nodes=0;number_of_boundary_nodes=0;number_of_ghost_nodes=0;
        ARRAY<int> number_of_individual_boundary_nodes(6);ARRAY<int> number_of_individual_ghost_nodes(6);
        
        RANGE<TV> interior_domain=Domain();interior_domain.Change_Size((T)-.5*Minimum_Edge_Length());
        RANGE<TV> regular_domain=Domain();regular_domain.Change_Size((T).5*Minimum_Edge_Length());
        Cell_Pointer_From_Index();
        for(int i=1;i<=node_iterator_deepest_cells.m;i++)if(node_iterator_deepest_cells(i)){
            TV node_location=Node_Location(nodes(i),cell_pointer_from_index(node_iterator_deepest_cells(i)));
            if(interior_domain.Lazy_Inside(node_location)) internal_nodes(++number_of_internal_nodes)=i;
            else{
                if(regular_domain.Lazy_Inside(node_location)){
                    if(node_location.x<interior_domain.min_corner.x) individual_side_boundary_nodes(1)(++number_of_individual_boundary_nodes(1))=i;
                    else if(node_location.x>interior_domain.max_corner.x) individual_side_boundary_nodes(2)(++number_of_individual_boundary_nodes(2))=i;
                    if(node_location.y<interior_domain.min_corner.y) individual_side_boundary_nodes(3)(++number_of_individual_boundary_nodes(3))=i;
                    else if(node_location.y>interior_domain.max_corner.y) individual_side_boundary_nodes(4)(++number_of_individual_boundary_nodes(4))=i;
                    if(node_location.z<interior_domain.min_corner.z) individual_side_boundary_nodes(5)(++number_of_individual_boundary_nodes(5))=i;
                    else if(node_location.z>interior_domain.max_corner.z) individual_side_boundary_nodes(6)(++number_of_individual_boundary_nodes(6))=i;
                    boundary_nodes(++number_of_boundary_nodes)=i;}
                else{
                    if(node_location.x<regular_domain.min_corner.x) individual_side_ghost_nodes(1)(++number_of_individual_ghost_nodes(1))=i;
                    else if(node_location.x>regular_domain.max_corner.x) individual_side_ghost_nodes(2)(++number_of_individual_ghost_nodes(2))=i;
                    if(node_location.y<regular_domain.min_corner.y) individual_side_ghost_nodes(3)(++number_of_individual_ghost_nodes(3))=i;
                    else if(node_location.y>regular_domain.max_corner.y) individual_side_ghost_nodes(4)(++number_of_individual_ghost_nodes(4))=i;
                    if(node_location.z<regular_domain.min_corner.z) individual_side_ghost_nodes(5)(++number_of_individual_ghost_nodes(5))=i;
                    else if(node_location.z>regular_domain.max_corner.z) individual_side_ghost_nodes(6)(++number_of_individual_ghost_nodes(6))=i;
                    ghost_nodes(++number_of_ghost_nodes)=i;}}}
                    
        for(int i=1;i<=6;i++)individual_side_boundary_nodes(i).Resize(number_of_individual_boundary_nodes(i));
        for(int i=1;i<=6;i++)individual_side_ghost_nodes(i).Resize(number_of_individual_ghost_nodes(i));
        ghost_nodes.Resize(number_of_ghost_nodes);
        node_iterator_data_up_to_date=true;}
}
//#####################################################################
// Function Face_Iterator_Data
//#####################################################################
namespace PhysBAM{template<class T> struct FACE_ITERATOR_DATA_HELPER{const OCTREE_GRID<T>* grid;int* number_of_faces;};}
template<class T> static void Face_Iterator_Data_Helper(void* data,const OCTREE_CELL<T>* cell1,const OCTREE_CELL<T>* cell2,const int axis)
{
    FACE_ITERATOR_DATA_HELPER<T>* helper=(FACE_ITERATOR_DATA_HELPER<T>*)data;
    const OCTREE_CELL<T>* cells[]={cell1,cell2};int deepest_cell,deepest_depth;MAP_OCTREE_MESH<T>::Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    const OCTREE_CELL<T>* smaller_cell=cells[deepest_cell];
    int face=MAP_OCTREE_MESH<T>::face_by_axis[axis][deepest_cell];int face_index=smaller_cell->Face(face);(*helper->number_of_faces)++;
    helper->grid->face_iterator_deepest_cells(face_index)=smaller_cell->Cell();helper->grid->faces(face_index)=face;
}
template<class T> void OCTREE_GRID<T>::
Face_Iterator_Data() const
{
    if(!face_iterator_data_up_to_date){
        face_iterator_deepest_cells.Resize(0);face_iterator_deepest_cells.Resize(number_of_faces);
        faces.Resize(number_of_faces);
        
        FACE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;
        
        int number_of_internal_faces=0;
        helper.number_of_faces=&number_of_internal_faces;
        MAP_OCTREE_MESH<T>::Map_Faces(uniform_grid,cells,0,&helper,Face_Iterator_Data_Helper);
        internal_faces.Resize(number_of_internal_faces);
        
        int number_of_boundary_faces=0;
        helper.number_of_faces=&number_of_boundary_faces;
        MAP_OCTREE_MESH<T>::Map_Boundary_Faces(uniform_grid,cells,&helper,Face_Iterator_Data_Helper);
        boundary_faces.Resize(number_of_boundary_faces);
        individual_side_boundary_faces.Resize(6);
        for(int i=1;i<=6;i++)individual_side_boundary_faces(i).Resize(number_of_boundary_faces);
        
        int number_of_ghost_faces=0;
        helper.number_of_faces=&number_of_ghost_faces;
        MAP_OCTREE_MESH<T>::Map_Ghost_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        MAP_OCTREE_MESH<T>::Map_Exterior_Ghost_Cell_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        ghost_faces.Resize(number_of_ghost_faces);
        individual_side_ghost_faces.Resize(6);
        for(int i=1;i<=6;i++)individual_side_ghost_faces(i).Resize(number_of_ghost_faces);
        individual_side_domain_ghost_faces.Resize(6);
        for(int i=1;i<=6;i++)individual_side_domain_ghost_faces(i).Resize(number_of_ghost_faces);
        
        number_of_internal_faces=0;number_of_boundary_faces=0;number_of_ghost_faces=0;
        ARRAY<int> number_of_individual_boundary_faces(6);ARRAY<int> number_of_individual_domain_ghost_faces(6);ARRAY<int> number_of_individual_ghost_faces(6);
        
        T tolerance=(T).25*Minimum_Edge_Length();
        RANGE<TV> domain=Domain();
        RANGE<TV> interior_domain=domain;interior_domain.Change_Size(-tolerance);
        RANGE<TV> regular_domain=domain;regular_domain.Change_Size(tolerance);
        Cell_Pointer_From_Index();
        for(int i=1;i<=face_iterator_deepest_cells.m;i++)if(face_iterator_deepest_cells(i)){
            TV face_location=Face_Location(faces(i),cell_pointer_from_index(face_iterator_deepest_cells(i)));
            if(interior_domain.Lazy_Inside(face_location)){internal_faces(++number_of_internal_faces)=i;}
            else{
                if(regular_domain.Lazy_Inside(face_location)){
                    if(face_location.x<interior_domain.min_corner.x) individual_side_boundary_faces(1)(++number_of_individual_boundary_faces(1))=i;
                    else if(face_location.x>interior_domain.max_corner.x) individual_side_boundary_faces(2)(++number_of_individual_boundary_faces(2))=i;
                    if(face_location.y<interior_domain.min_corner.y) individual_side_boundary_faces(3)(++number_of_individual_boundary_faces(3))=i;
                    else if(face_location.y>interior_domain.max_corner.y) individual_side_boundary_faces(4)(++number_of_individual_boundary_faces(4))=i;
                    if(face_location.z<interior_domain.min_corner.z) individual_side_boundary_faces(5)(++number_of_individual_boundary_faces(5))=i;
                    else if(face_location.z>interior_domain.max_corner.z) individual_side_boundary_faces(6)(++number_of_individual_boundary_faces(6))=i;
                    boundary_faces(++number_of_boundary_faces)=i;}
                else{
                    if(face_location.x<regular_domain.min_corner.x) individual_side_ghost_faces(1)(++number_of_individual_ghost_faces(1))=i;
                    else if(face_location.x>regular_domain.max_corner.x) individual_side_ghost_faces(2)(++number_of_individual_ghost_faces(2))=i;
                    if(face_location.y<regular_domain.min_corner.y) individual_side_ghost_faces(3)(++number_of_individual_ghost_faces(3))=i;
                    else if(face_location.y>regular_domain.max_corner.y) individual_side_ghost_faces(4)(++number_of_individual_ghost_faces(4))=i;
                    if(face_location.z<regular_domain.min_corner.z) individual_side_ghost_faces(5)(++number_of_individual_ghost_faces(5))=i;
                    else if(face_location.z>regular_domain.max_corner.z) individual_side_ghost_faces(6)(++number_of_individual_ghost_faces(6))=i;
                    if(abs(face_location.x-domain.min_corner.x)<tolerance) individual_side_domain_ghost_faces(1)(++number_of_individual_domain_ghost_faces(1))=i;
                    else if(abs(face_location.x-domain.max_corner.x)<tolerance) individual_side_domain_ghost_faces(2)(++number_of_individual_domain_ghost_faces(2))=i;
                    if(abs(face_location.y-domain.min_corner.y)<tolerance) individual_side_domain_ghost_faces(3)(++number_of_individual_domain_ghost_faces(3))=i;
                    else if(abs(face_location.y-domain.max_corner.y)<tolerance) individual_side_domain_ghost_faces(4)(++number_of_individual_domain_ghost_faces(4))=i;
                    if(abs(face_location.z-domain.min_corner.z)<tolerance) individual_side_domain_ghost_faces(5)(++number_of_individual_domain_ghost_faces(5))=i;
                    else if(abs(face_location.z-domain.max_corner.z)<tolerance) individual_side_domain_ghost_faces(6)(++number_of_individual_domain_ghost_faces(6))=i;
                    ghost_faces(++number_of_ghost_faces)=i;}}}
                    
        for(int i=1;i<=6;i++)individual_side_boundary_faces(i).Resize(number_of_individual_boundary_faces(i));
        for(int i=1;i<=6;i++)individual_side_domain_ghost_faces(i).Resize(number_of_individual_domain_ghost_faces(i));
        for(int i=1;i<=6;i++)individual_side_ghost_faces(i).Resize(number_of_individual_ghost_faces(i));
        ghost_faces.Resize(number_of_ghost_faces);
        face_iterator_data_up_to_date=true;}
}
//#####################################################################
// Function Calculate_Cell_Pointer_From_Index_Array
//#####################################################################
template<class T> static void Calculate_Cell_Pointer_From_Index_Array_Helper(OCTREE_CELL<T>* cell,ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index)
{
    assert(cell_pointer_from_index(cell->Cell())==0);cell_pointer_from_index(cell->Cell())=cell;
    if(cell->Has_Children()) for(int i=0;i<8;i++) Calculate_Cell_Pointer_From_Index_Array_Helper(cell->Child(i),cell_pointer_from_index);
}
template<class T> void OCTREE_GRID<T>::
Calculate_Cell_Pointer_From_Index_Array(ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index)const 
{
    cell_pointer_from_index.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_pointer_from_index,0);
    for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
        for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j++)
            for(int ij=1-number_of_ghost_cells;ij<=uniform_grid.numbers_of_cells.z+number_of_ghost_cells;ij++)
                Calculate_Cell_Pointer_From_Index_Array_Helper(cells(i,j,ij),cell_pointer_from_index);
}
//#####################################################################
// Function Calculate_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Neighbors_Array_Helper(void* data,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int axis)
{
    ARRAY<VECTOR<OCTREE_CELL<T>*,6> >* neighbors=(ARRAY<VECTOR<OCTREE_CELL<T>*,6> >*)data;
    const OCTREE_CELL<T>* cells[2]={c1,c2};
    int smaller_cell=0;int larger_cell=1;if(cells[0]->Depth_Of_This_Cell()<cells[1]->Depth_Of_This_Cell()){smaller_cell=1;larger_cell=0;}
    int min_depth=cells[larger_cell]->Depth_Of_This_Cell();
    while(cells[smaller_cell]->Depth_Of_This_Cell()>min_depth){ // mark the large cell as the neighbor of all of the cells on the path up to the cell of the same level as the large cell
        (*neighbors)(cells[smaller_cell]->Cell())(axis*2+1+(1-smaller_cell))=(OCTREE_CELL<T>*)cells[larger_cell];
        cells[smaller_cell]=cells[smaller_cell]->Parent();}
    (*neighbors)(cells[0]->Cell())(axis*2+2)=(OCTREE_CELL<T>*)cells[1]; // mark the two cells that are the same size as neighbors of each other
    (*neighbors)(cells[1]->Cell())(axis*2+1)=(OCTREE_CELL<T>*)cells[0];
}
template<class T> void OCTREE_GRID<T>::
Calculate_Neighbors_Array(ARRAY<VECTOR<OCTREE_CELL<T>*,6> >& neighbors)const
{
    neighbors.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(neighbors.Flattened(),0);
    MAP_OCTREE_MESH<T>::Map_Faces(uniform_grid,cells,number_of_ghost_cells,&neighbors,Calculate_Neighbors_Array_Helper);
}
//#####################################################################
// Function Calculate_Node_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Node_Neighbors_Array_Helper(void* data,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const int axis)
{
    ARRAY<VECTOR<int,6> >* node_neighbors=(ARRAY<VECTOR<int,6> >*)data;
    const OCTREE_CELL<T>* cells[4]={c1,c2,c3,c4};
    int deepest_cell,deepest_depth;MAP_OCTREE_MESH<T>::Get_Deepest_Cell(4,cells,deepest_cell,deepest_depth);
    int node_index1=cells[deepest_cell]->Node(MAP_OCTREE_MESH<T>::edge_nodes_by_axis[axis][deepest_cell][0]),
        node_index2=cells[deepest_cell]->Node(MAP_OCTREE_MESH<T>::edge_nodes_by_axis[axis][deepest_cell][1]);
    // add each other as neighbors (in the correct directions)
    int neighbor_position1=2*axis+1,neighbor_position2=2*axis+2;
    assert((*node_neighbors)(node_index1)(neighbor_position2)==0);assert((*node_neighbors)(node_index2)(neighbor_position1)==0);
    (*node_neighbors)(node_index1)(neighbor_position2)=node_index2;(*node_neighbors)(node_index2)(neighbor_position1)=node_index1;
}
template<class T> void OCTREE_GRID<T>::
Calculate_Node_Neighbors_Array(ARRAY<VECTOR<int,6> >& node_neighbors)const
{
    node_neighbors.Resize(number_of_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(node_neighbors.Flattened(),0);
    MAP_OCTREE_MESH<T>::Map_Edges(uniform_grid,cells,number_of_ghost_cells,&node_neighbors,Calculate_Node_Neighbors_Array_Helper);
}
//#####################################################################
// Function Calculate_Fully_Refined_Block
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Calculate_Fully_Refined_Block_Array(ARRAY<bool>& fully_refined_block) const
{
    fully_refined_block.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(fully_refined_block,false);Neighbors();
    OCTREE_CELL<T>* cells[number_of_children_per_cell];int lookup[][2]={{0,0},{2,0},{4,0},{4,1},{6,0},{6,1},{6,2},{6,3}};
    for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(*this,3);iterator.Valid();iterator.Next()){
        if(iterator.Cell_Pointer()->Depth_Of_This_Cell()!=maximum_depth)continue;
        bool block_fully_refined=true;cells[0]=iterator.Cell_Pointer();
        for(int i=1;i<number_of_cells_per_block;i++){
            cells[i]=neighbors(cells[lookup[i][1]]->Cell())(lookup[i][0]);
            if(!cells[i]||cells[i]->Depth_Of_This_Cell()!=maximum_depth){block_fully_refined=false;break;}}
        fully_refined_block(iterator.Cell_Index())=block_fully_refined;}
}
//#####################################################################
// Function Enslave_T_Junction_Nodes
//#####################################################################
// enslaves nodes between levels to be linear combinations of the corner nodes of the largest cell touching the given node
namespace PhysBAM{template<class T,class T2> struct ENSLAVE_T_JUNCTION_NODES_HELPER{OCTREE_GRID<T>* grid;ARRAY<T2>* nodes;int depth_to_enforce;};}
template<class T,class T2> static void Enslave_T_Junction_Nodes_Helper(void* data,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2, 
    const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8)
{
    ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2>* helper_struct=(ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2>*)data;
    const OCTREE_CELL<T>* cells[]={c1,c2,c3,c4,c5,c6,c7,c8};
    int deepest_cell,maximum_depth,shallowest_cell,minimum_depth;
    MAP_OCTREE_MESH<T>::Get_Deepest_And_Shallowest_Cell(8,cells,deepest_cell,maximum_depth,shallowest_cell,minimum_depth);
    if(helper_struct->depth_to_enforce>=0&&minimum_depth!=helper_struct->depth_to_enforce) return;
    int current_node=cells[deepest_cell]->Node(7-deepest_cell);
    // check if we are in the corner of the shallowest cell -  if so, we are guaranteed not to be a node that needs to be a linear combination of the corner nodes
    if(cells[shallowest_cell]->Node(7-shallowest_cell)==current_node) return;
    // figure out what our value should be based on our position in the cell
    VECTOR<T,3> cell_location=helper_struct->grid->Node_Location(7-deepest_cell,cells[deepest_cell]);
    (*(helper_struct->nodes))(current_node)=LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(*helper_struct->grid,cells[shallowest_cell],*helper_struct->nodes,cell_location);
}
template<class T> template<class T2> void OCTREE_GRID<T>::
Enslave_T_Junction_Nodes(ARRAY<T2>* nodes,int depth_to_enforce)
{
    ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2> helper;helper.nodes=nodes;helper.grid=this;helper.depth_to_enforce=depth_to_enforce;
    MAP_OCTREE_MESH<T>::Map_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Enslave_T_Junction_Nodes_Helper<T,T2>);
}
//#####################################################################
// Function Check_Tree_Consistency
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Check_Tree_Consistency(bool check_cells,bool check_nodes,bool check_faces,bool check_neighbor_structures)
{
    Tree_Topology_Changed();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(check_cells){
        LOG::cout<<"Checking cell indices"<<std::endl;
        ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index=Cell_Pointer_From_Index();
        for(int i=1;i<=number_of_cells;i++)for(int j=1;j<=number_of_cells;j++)if(i!=j){
            if(cell_pointer_from_index(i)==cell_pointer_from_index(j)) LOG::cout<<"Different index pointing to the same cell: "<<i<<", "<<j<<std::endl;
            if((cell_pointer_from_index(i)->Center()-cell_pointer_from_index(j)->Center()).Magnitude_Squared()<1E-5)LOG::cout<<"Different cells on top of eachother: "<<i<<", "<<j<<std::endl;}}
    if(check_nodes){
        LOG::cout<<"Checking nodes"<<std::endl;
        for(int i=1;i<=number_of_nodes;i++)for(int j=1;j<=number_of_nodes;j++)if(i!=j){
        if((Node_Location(i)-Node_Location(j)).Magnitude_Squared()<1E-5) LOG::cout<<"Different nodes sitting on top of eachother: "<<i<<", "<<j<<". Location: "<<Node_Location(i)<<std::endl;}}
    if(check_faces){
        LOG::cout<<"Checking faces"<<std::endl;
        for(int i=1;i<=number_of_faces;i++)for(int j=1;j<=number_of_faces;j++)if(i!=j){
        if((Face_Location(i)-Face_Location(j)).Magnitude_Squared()<1E-5) LOG::cout<<"Different faces sitting on top of eachother: "<<i<<", "<<j<<". Location: "<<Face_Location(i)<<std::endl;}}
#endif
}
//#####################################################################
// Function Print_Storage_Requirement
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Print_Storage_Requirement()
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    int node_iterator_size=sizeof(int)*(internal_nodes.m+boundary_nodes.m+ghost_nodes.m+node_iterator_deepest_cells.m);
    for(int i=1;i<=individual_side_boundary_nodes.m;i++)node_iterator_size+=sizeof(int)+individual_side_boundary_nodes(i).m;
    for(int i=1;i<=individual_side_ghost_nodes.m;i++)node_iterator_size+=sizeof(int)+individual_side_ghost_nodes(i).m;
    node_iterator_size+=sizeof(unsigned char)*nodes.m;
    LOG::cout<<"Node Iterator: "<<node_iterator_size/(1024*1024)<<" mb"<<std::endl;

    int face_iterator_size=sizeof(int)*(internal_faces.m+boundary_faces.m+ghost_faces.m+face_iterator_deepest_cells.m);
    for(int i=1;i<=individual_side_boundary_faces.m;i++)face_iterator_size+=sizeof(int)+individual_side_boundary_faces(i).m;
    for(int i=1;i<=individual_side_ghost_faces.m;i++)face_iterator_size+=sizeof(int)+individual_side_ghost_faces(i).m;
    face_iterator_size+=sizeof(unsigned char)*faces.m;
    LOG::cout<<"Face Iterator: "<<face_iterator_size/(1024*1024)<<" mb"<<std::endl;

    int allocated_children_size=sizeof(OCTREE_CHILDREN<T>*)*allocated_children.m;
    LOG::cout<<"Allocated Children Buffer: "<<allocated_children_size/(1024*1024)<<" mb"<<std::endl;

    int cell_pointer_from_index_size=sizeof(OCTREE_CHILDREN<T>*)*cell_pointer_from_index.m;
    LOG::cout<<"Cell Pointer From Index: "<<cell_pointer_from_index_size/(1024*1024)<<" mb"<<std::endl;

    int neighbors_size=sizeof(OCTREE_CHILDREN<T>*)*neighbors.m;
    LOG::cout<<"Neighbors Size: "<<neighbors_size/(1024*1024)<<" mb"<<std::endl;

    int node_neighbors_size=sizeof(int)*node_neighbors.m;
    LOG::cout<<"Node Neighbors: "<<node_neighbors_size/(1024*1024)<<" mb"<<std::endl;

    int tree_size=sizeof(OCTREE_CELL<T>*)*cells.counts.Product();
    for(int i=1;i<=allocated_children.m;i++)tree_size+=allocated_children(i)->Storage_Requirement();
    LOG::cout<<"Tree_Size: "<<tree_size/(1024*1024)<<" mb"<<std::endl;

    int total_size=node_iterator_size+face_iterator_size+allocated_children_size+cell_pointer_from_index_size+neighbors_size+node_neighbors_size+tree_size;
    LOG::cout<<"Total octree storage size: "<<total_size/(1<<20)<<"MB"<<std::endl;
#endif
}
//#####################################################################
// Function Move_Contents_Left_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Left_Two()
{
    int start=1-number_of_ghost_cells;TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(start,start+1,start,end.y,start,end.z);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(-2,0,0));temp_data.domain.min_corner.x=end.x-1;temp_data.domain.max_corner.x=end.x;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int j=start;j<=end.y;j++)for(int ij=start;ij<=end.z;ij++){assert(!cells(end.x-1,j,ij)->Has_Children()&&!cells(end.x,j,ij)->Has_Children());
        cells(end.x,j,ij)->Node(1)=cells(end.x-1,j,ij)->Node(0);cells(end.x,j,ij)->Node(3)=cells(end.x-1,j,ij)->Node(2);
        cells(end.x,j,ij)->Node(5)=cells(end.x-1,j,ij)->Node(4);cells(end.x,j,ij)->Node(7)=cells(end.x-1,j,ij)->Node(6);}
    for(int j=start;j<=end.y;j++)for(int ij=start;ij<=end.z;ij++){
        cells(end.x-1,j,ij)->Node(0)=cells(end.x-2,j,ij)->Node(1);cells(end.x-1,j,ij)->Node(2)=cells(end.x-2,j,ij)->Node(3);
        cells(end.x-1,j,ij)->Node(4)=cells(end.x-2,j,ij)->Node(5);cells(end.x-1,j,ij)->Node(6)=cells(end.x-2,j,ij)->Node(7);
        cells(end.x,j,ij)->Face(1)=cells(end.x-1,j,ij)->Face(0);cells(end.x-1,j,ij)->Face(0)=cells(end.x-2,j,ij)->Face(1);}
    for(int j=start;j<=end.y;j+=2)for(int ij=start;ij<=end.z;ij+=2)cells(end.x,j,ij)->owner->parents_center=uniform_grid.Node(end.x+2,j+1,ij+1);
}
//#####################################################################
// Function Move_Contents_Right_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Right_Two()
{
    int start=1-number_of_ghost_cells; TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(end.x-1,end.x,start,end.y,start,end.z);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(2,0,0));temp_data.domain.min_corner.x=start;temp_data.domain.max_corner.x=start+1;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int j=start;j<=end.y;j++)for(int ij=start;ij<=end.z;ij++){assert(!cells(start,j,ij)->Has_Children()&&!cells(start+1,j,ij)->Has_Children());
        cells(start,j,ij)->Node(0)=cells(start+1,j,ij)->Node(1);cells(start,j,ij)->Node(2)=cells(start+1,j,ij)->Node(3);
        cells(start,j,ij)->Node(4)=cells(start+1,j,ij)->Node(5);cells(start,j,ij)->Node(6)=cells(start+1,j,ij)->Node(7);}
    for(int j=start;j<=end.y;j++)for(int ij=start;ij<=end.z;ij++){
        cells(start+1,j,ij)->Node(1)=cells(start+2,j,ij)->Node(0);cells(start+1,j,ij)->Node(3)=cells(start+2,j,ij)->Node(2);
        cells(start+1,j,ij)->Node(5)=cells(start+2,j,ij)->Node(4);cells(start+1,j,ij)->Node(7)=cells(start+2,j,ij)->Node(6);
        cells(start,j,ij)->Face(0)=cells(start+1,j,ij)->Face(1);cells(start+1,j,ij)->Face(1)=cells(start+2,j,ij)->Face(0);}
    for(int j=start;j<=end.y;j+=2)for(int ij=start;ij<=end.z;ij+=2)cells(start,j,ij)->owner->parents_center=uniform_grid.Node(start-1,j+1,ij+1);
}
//#####################################################################
// Function Move_Contents_Down_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Down_Two()
{
    int start=1-number_of_ghost_cells; TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(start,end.x,start,start+1,start,end.z);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,-2,0));temp_data.domain.min_corner.y=end.y-1;temp_data.domain.max_corner.y=end.y;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int i=start;i<=end.x;i++)for(int ij=start;ij<=end.z;ij++){assert(!cells(i,end.y-1,ij)->Has_Children());assert(!cells(i,end.y,ij)->Has_Children());
        cells(i,end.y,ij)->Node(2)=cells(i,end.y-1,ij)->Node(0);cells(i,end.y,ij)->Node(3)=cells(i,end.y-1,ij)->Node(1);
        cells(i,end.y,ij)->Node(6)=cells(i,end.y-1,ij)->Node(4);cells(i,end.y,ij)->Node(7)=cells(i,end.y-1,ij)->Node(5);}
    for(int i=start;i<=end.x;i++)for(int ij=start;ij<=end.z;ij++){
        cells(i,end.y-1,ij)->Node(0)=cells(i,end.y-2,ij)->Node(2);cells(i,end.y-1,ij)->Node(1)=cells(i,end.y-2,ij)->Node(3);
        cells(i,end.y-1,ij)->Node(4)=cells(i,end.y-2,ij)->Node(6);cells(i,end.y-1,ij)->Node(5)=cells(i,end.y-2,ij)->Node(7);
        cells(i,end.y,ij)->Face(3)=cells(i,end.y-1,ij)->Face(2);cells(i,end.y-1,ij)->Face(2)=cells(i,end.y-2,ij)->Face(3);}
    for(int i=start;i<=end.x;i+=2)for(int ij=start;ij<=end.z;ij+=2)cells(i,end.y,ij)->owner->parents_center=uniform_grid.Node(i+1,end.y+2,ij+1);
}
//#####################################################################
// Function Move_Contents_Up_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Up_Two()
{
    int start=1-number_of_ghost_cells; TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(start,end.x,end.y-1,end.y,start,end.z);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,2,0));temp_data.domain.min_corner.y=start;temp_data.domain.max_corner.y=start+1;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int i=start;i<=end.x;i++)for(int ij=start;ij<=end.z;ij++){assert(!cells(i,start,ij)->Has_Children());assert(!cells(i,start+1,ij)->Has_Children());
        cells(i,start,ij)->Node(0)=cells(i,start+1,ij)->Node(2);cells(i,start,ij)->Node(1)=cells(i,start+1,ij)->Node(3);
        cells(i,start,ij)->Node(4)=cells(i,start+1,ij)->Node(6);cells(i,start,ij)->Node(5)=cells(i,start+1,ij)->Node(7);}
    for(int i=start;i<=end.x;i++)for(int ij=start;ij<=end.z;ij++){
        cells(i,start+1,ij)->Node(2)=cells(i,start+2,ij)->Node(0);cells(i,start+1,ij)->Node(3)=cells(i,start+2,ij)->Node(1);
        cells(i,start+1,ij)->Node(6)=cells(i,start+2,ij)->Node(4);cells(i,start+1,ij)->Node(7)=cells(i,start+2,ij)->Node(5);
        cells(i,start,ij)->Face(2)=cells(i,start+1,ij)->Face(3);cells(i,start+1,ij)->Face(3)=cells(i,start+2,ij)->Face(2);}
    for(int i=start;i<=end.x;i+=2)for(int ij=start;ij<=end.z;ij+=2)cells(i,start,ij)->owner->parents_center=uniform_grid.Node(i+1,start-1,ij+1);
}
//#####################################################################
// Function Move_Contents_Forward_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Forward_Two()
{
    int start=1-number_of_ghost_cells; TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(start,end.x,start,end.y,start,start+1);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,0,-2));temp_data.domain.min_corner.z=end.z-1;temp_data.domain.max_corner.z=end.z;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int i=start;i<=end.x;i++)for(int j=start;j<=end.y;j++){assert(!cells(i,j,end.z-1)->Has_Children());assert(!cells(i,j,end.z)->Has_Children());
        cells(i,j,end.z)->Node(4)=cells(i,j,end.z-1)->Node(0);cells(i,j,end.z)->Node(5)=cells(i,j,end.z-1)->Node(1);
        cells(i,j,end.z)->Node(6)=cells(i,j,end.z-1)->Node(2);cells(i,j,end.z)->Node(7)=cells(i,j,end.z-1)->Node(3);}
    for(int i=start;i<=end.x;i++)for(int j=start;j<=end.y;j++){
        cells(i,j,end.z-1)->Node(0)=cells(i,j,end.z-2)->Node(4);cells(i,j,end.z-1)->Node(1)=cells(i,j,end.z-2)->Node(5);
        cells(i,j,end.z-1)->Node(2)=cells(i,j,end.z-2)->Node(6);cells(i,j,end.z-1)->Node(3)=cells(i,j,end.z-2)->Node(7);
        cells(i,j,end.z)->Face(5)=cells(i,j,end.z-1)->Face(4);cells(i,j,end.z-1)->Face(4)=cells(i,j,end.z-2)->Face(5);}
    for(int i=start;i<=end.x;i+=2)for(int j=start;j<=end.y;j+=2)cells(i,j,end.z)->owner->parents_center=uniform_grid.Node(i+1,j+1,end.z+2);
}
//#####################################################################
// Function Move_Contents_Backward_Two
//#####################################################################
template<class T> void OCTREE_GRID<T>::
Move_Contents_Backward_Two()
{
    int start=1-number_of_ghost_cells; TV_INT end=uniform_grid.numbers_of_cells+number_of_ghost_cells;
    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > temp_data(start,end.x,start,end.y,end.z-1,end.z);ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,0,2));temp_data.domain.min_corner.z=start;temp_data.domain.max_corner.z=start+1;temp_data.Calculate_Acceleration_Constants();ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >::Put(temp_data,cells);
    for(int i=start;i<=end.x;i++)for(int j=start;j<=end.y;j++){assert(!cells(i,j,start)->Has_Children());assert(!cells(i,j,start+1)->Has_Children());
        cells(i,j,start)->Node(0)=cells(i,j,start+1)->Node(4);cells(i,j,start)->Node(1)=cells(i,j,start+1)->Node(5);
        cells(i,j,start)->Node(2)=cells(i,j,start+1)->Node(6);cells(i,j,start)->Node(3)=cells(i,j,start+1)->Node(7);}
    for(int i=start;i<=end.x;i++)for(int j=start;j<=end.y;j++){
        cells(i,j,start+1)->Node(4)=cells(i,j,start+2)->Node(0);cells(i,j,start+1)->Node(5)=cells(i,j,start+2)->Node(1);
        cells(i,j,start+1)->Node(6)=cells(i,j,start+2)->Node(2);cells(i,j,start+1)->Node(7)=cells(i,j,start+2)->Node(3);
        cells(i,j,start)->Face(4)=cells(i,j,start+1)->Face(5);cells(i,j,start+1)->Face(5)=cells(i,j,start+2)->Face(4);}
    for(int i=start;i<=end.x;i+=2)for(int j=start;j<=end.y;j+=2)cells(i,j,start)->owner->parents_center=uniform_grid.Node(i+1,j+1,start-1);
}
//#####################################################################
template class OCTREE_GRID<float>;
template void OCTREE_GRID<float>::Enslave_T_Junction_Nodes(ARRAY<float>* nodes,int depth_to_enforce);
template void OCTREE_GRID<float>::Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTREE_GRID<double>;
template void OCTREE_GRID<double>::Enslave_T_Junction_Nodes(ARRAY<double>* nodes,int depth_to_enforce);
template void OCTREE_GRID<double>::Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce);
#endif
#endif
