//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
using namespace PhysBAM;
//#####################################################################
template<class T> const char* QUADTREE_GRID<T>::name="quadtree";
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUADTREE_GRID<T>::
QUADTREE_GRID()
    :uniform_grid(2,2,0,1,0,1),minimum_cell_size(0),number_of_ghost_cells(0),number_of_cells(0),number_of_nodes(0),number_of_faces(0),maximum_depth(1)
{
    Tree_Topology_Changed();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> QUADTREE_GRID<T>::
~QUADTREE_GRID()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Initialize(const GRID<TV> uniform_grid_input,const int maximum_depth_input,const int number_of_ghost_cells_input,const bool use_nodes,const bool use_faces)
{
    int i,j;uniform_grid=uniform_grid_input;Set_Maximum_Depth(maximum_depth_input);number_of_ghost_cells=number_of_ghost_cells_input;
    if(uniform_grid.numbers_of_cells.x%2 || uniform_grid.numbers_of_cells.y%2){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Only even sized grids are supported (grid = "<<uniform_grid.numbers_of_cells<<")"<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
    number_of_cells=number_of_nodes=number_of_faces=0;
    for(i=1;i<=allocated_children.m;i++)delete allocated_children(i);allocated_children.Resize(0);
    cells.Resize(RANGE<TV_INT>(TV_INT()+(1-number_of_ghost_cells),uniform_grid.numbers_of_cells+number_of_ghost_cells));
    ARRAY<int,VECTOR<int,2> > nodes(cells.domain.min_corner.x,cells.domain.max_corner.x+1,cells.domain.min_corner.y,cells.domain.max_corner.y+1);
    ARRAY<int,VECTOR<int,2> > faces_u(cells.domain.min_corner.x,cells.domain.max_corner.x+1,cells.domain.min_corner.y,cells.domain.max_corner.y),faces_v(cells.domain.min_corner.x,cells.domain.max_corner.x,cells.domain.min_corner.y,cells.domain.max_corner.y+1);
    if(use_nodes) for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y+1;j++)nodes(i,j)=++number_of_nodes;
    if(use_faces){
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y;j++)faces_u(i,j)=++number_of_faces;
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i++)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y+1;j++)faces_v(i,j)=++number_of_faces;}
    ARRAY<int,VECTOR<int,1> > nodes_to_input(0,8);ARRAY<int,VECTOR<int,1> > faces_to_input(0,11);
    for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i+=2)for(j=cells.domain.min_corner.y;j<=cells.domain.max_corner.y;j+=2){
        QUADTREE_CHILDREN<T>* owner=new QUADTREE_CHILDREN<T>(use_nodes,use_faces);
        allocated_children.Append(owner);
        TV center=uniform_grid.Node(i+1,j+1),DX=uniform_grid.dX;
        if(use_nodes) for(int ii=i;ii<=i+2;ii++)for(int jj=j;jj<=j+2;jj++)nodes_to_input((ii-i)+(jj-j)*3)=nodes(ii,jj);
        if(use_faces){
            for(int ii=i;ii<=i+2;ii++)for(int jj=j;jj<=j+1;jj++){faces_to_input((ii-i)+(jj-j)*3)=faces_u(ii,jj);}
            for(int ii=i;ii<=i+1;ii++)for(int jj=j;jj<=j+2;jj++){faces_to_input(6+(jj-j)+(ii-i)*3)=faces_v(ii,jj);}}
        owner->Initialize_Pseudo_Root_Cells(number_of_cells,nodes_to_input,faces_to_input,center,DX);
        cells(i,j)=&owner->children[0];cells(i+1,j)=&owner->children[1];cells(i,j+1)=&owner->children[2];cells(i+1,j+1)=&owner->children[3];}
}
//#####################################################################
// Function Tree_Topology_Changed
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
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
template<class T> QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Leaf_Cell(const TV& location,const T thickness) const
{
    TV_INT i=uniform_grid.Cell(location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    if(cells.Domain_Indices().Lazy_Outside(i)) return 0;
    return Inside_Offspring(cells(i),location,thickness);
}
//#####################################################################
// Function Clamped_Leaf_Cell
//#####################################################################
template<class T> QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Clamped_Leaf_Cell(const TV& location,const T thickness) const
{
    TV clamped_location=uniform_grid.Clamp(location,number_of_ghost_cells);
    TV_INT i=uniform_grid.Cell(clamped_location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    i=TV_INT::Componentwise_Min(i,TV_INT(cells.domain.max_corner.x,cells.domain.max_corner.y)); // need this clamp because Cell could potentially return i or ij beyond last index
    return Inside_Offspring(cells(i),clamped_location,thickness);
}
//#####################################################################
// Function Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> const QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(const QUADTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell,tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else if(location.y<center.y-DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(3);
        else if(location.y>center.y+DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(4);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Function Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(QUADTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell,tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else if(location.y<center.y-DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(3);
        else if(location.y>center.y+DX_over_two.y) starting_cell=neighbors(starting_cell->Cell())(4);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Function Inside_Cell
//#####################################################################
template<class T> bool QUADTREE_GRID<T>::
Inside_Cell(const QUADTREE_CELL<T>* cell,const TV& location) const
{
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    if(location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x && location.y>=center.y-DX_over_two.y && location.y<=center.y+DX_over_two.y) return true;
    else return false;
}
//#####################################################################
// Function Inside_Thickened_Cell
//#####################################################################
template<class T> bool QUADTREE_GRID<T>::
Inside_Thickened_Cell(const QUADTREE_CELL<T>* cell,const TV& location,const T thickness) const
{
    TV center(cell->Center()),DX_over_two(((T).5+thickness)*cell->DX());
    if(location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x && location.y>=center.y-DX_over_two.y && location.y<=center.y+DX_over_two.y) return true;
    else return false;
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Inside_Offspring(QUADTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        int child_to_use=0;if(location.x>center.x)child_to_use+=1;if(location.y>center.y)child_to_use+=2;
        cell=cell->Child(child_to_use);}
    return cell;
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> const QUADTREE_CELL<T>* QUADTREE_GRID<T>::
Inside_Offspring(const QUADTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        int child_to_use=0;if(location.x>center.x)child_to_use+=1;if(location.y>center.y)child_to_use+=2;
        cell=cell->Child(child_to_use);}
    return cell;
}
//#####################################################################
// Function Refine_Cell
//#####################################################################
// refines the given cell so that the node at the given location is touching a refined cell instead of a coarse one
// the location needs to be inside or on one of the faces of the cell
template<class T> void QUADTREE_GRID<T>::
Refine_Cell(const int max_depth,QUADTREE_CELL<T>* cell,const TV& location,ARRAY<QUADTREE_CELL<T>*>* new_cells,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_nodes, 
            ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_faces,const QUADTREE_GRID<T>* grid)
{
    int depth_of_cell=cell->Depth_Of_This_Cell();assert(depth_of_cell<=maximum_depth);
    if(depth_of_cell>=max_depth) return;
    if(!Inside_Thickened_Cell(cell,location)) return; // make sure that the node is actually in the cell
    if(!cell->Has_Children()) cell->Create_Children(number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,grid);
    for(int i=0;i<4;i++) Refine_Cell(max_depth,cell->Child(i),location, new_cells,new_nodes,new_faces,grid); // recur on the children
}
//#####################################################################
// Function Compact_Array_Indices
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array,ARRAY<int>* face_mapping_array)
{
    if(cell_mapping_array){cell_mapping_array->Resize(number_of_cells,false);ARRAYS_COMPUTATIONS::Fill(*cell_mapping_array,0);number_of_cells=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j++)
                cells(i,j)->Create_Cell_Compaction_Mapping(*cell_mapping_array,number_of_cells);}
    if(face_mapping_array){face_mapping_array->Resize(number_of_faces,false);ARRAYS_COMPUTATIONS::Fill(*face_mapping_array,0);number_of_faces=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j+=2)
                cells(i,j)->owner->Create_Face_Compaction_Mapping_Helper(*face_mapping_array,number_of_faces);}
    if(node_mapping_array){node_mapping_array->Resize(number_of_nodes,false);ARRAYS_COMPUTATIONS::Fill(*node_mapping_array,0);number_of_nodes=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
            for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j+=2)
                cells(i,j)->owner->Create_Node_Compaction_Mapping_Helper(*node_mapping_array,number_of_nodes);}
}
//#####################################################################
// Function Refine_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Refine_Cells_Intersecting_Box_Helper(QUADTREE_GRID<T>& grid,QUADTREE_CELL<T>& cell,const RANGE<VECTOR<T,2> >& box,ARRAY<QUADTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    if(cell.Depth_Of_This_Cell()>=refinement_depth || !box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){ // refine
        refined_cells.Append(&cell);cell.Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}
    for(int i=0;i<4;i++) Refine_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,refined_cells,refinement_depth);
}
template<class T> void QUADTREE_GRID<T>::
Refine_Cells_Intersecting_Box(const RANGE<VECTOR<T,2> >& box,ARRAY<QUADTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,QUADTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,QUADTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index+=TV_INT(1,1);
    int depth=maximum_depth;if(refinement_depth) depth=refinement_depth;
    for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) Refine_Cells_Intersecting_Box_Helper(*this,*cells(i,j),box,refined_cells,refinement_depth);
}
//#####################################################################
// Function Get_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Get_Cells_Intersecting_Box_Helper(QUADTREE_GRID<T>& grid,QUADTREE_CELL<T>& cell,const RANGE<VECTOR<T,2> >& box,ARRAY<QUADTREE_CELL<T>*>& intersecting_cells)
{   
    if(!box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){
        intersecting_cells.Append(&cell);}
    else for(int i=0;i<4;i++) Get_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,intersecting_cells);
}   
template<class T> void QUADTREE_GRID<T>::
Get_Cells_Intersecting_Box(const RANGE<VECTOR<T,2> >& box,ARRAY<QUADTREE_CELL<T>*>& intersecting_cells)
{   
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,QUADTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,QUADTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index+=TV_INT(1,1);
    for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++)
        Get_Cells_Intersecting_Box_Helper(*this,*cells(i,j),box,intersecting_cells);
}
//#####################################################################
// Function Node_Iterator_Data
//#####################################################################
template<class T> struct NODE_ITERATOR_DATA_HELPER{const QUADTREE_GRID<T>* grid;int* number_of_nodes;};
template<class T> static void Node_Iterator_Data_Helper(void* data,const QUADTREE_CELL<T>* cell1,const QUADTREE_CELL<T>* cell2,const QUADTREE_CELL<T>* cell3,const QUADTREE_CELL<T>* cell4)
{
    NODE_ITERATOR_DATA_HELPER<T>* helper=(NODE_ITERATOR_DATA_HELPER<T>*)data;
    const QUADTREE_CELL<T>* cells[4]={cell1,cell2,cell3,cell4};int maximum_depth,deepest_cell;MAP_QUADTREE_MESH<T>::Get_Deepest_Cell(4,cells,deepest_cell,maximum_depth);
    int node=3-deepest_cell;int node_index=cells[deepest_cell]->Node(node);(*helper->number_of_nodes)++;
    helper->grid->node_iterator_deepest_cells(node_index)=cells[deepest_cell]->Cell();helper->grid->nodes(node_index)=node;
}
template<class T> void QUADTREE_GRID<T>::
Node_Iterator_Data() const
{
    if(!node_iterator_data_up_to_date){
        node_iterator_deepest_cells.Resize(0);node_iterator_deepest_cells.Resize(number_of_nodes);
        nodes.Resize(number_of_nodes);
        
        NODE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;
        
        int number_of_internal_nodes=0;
        helper.number_of_nodes=&number_of_internal_nodes;
        MAP_QUADTREE_MESH<T>::Map_Nodes(uniform_grid,cells,0,&helper,Node_Iterator_Data_Helper);
        internal_nodes.Resize(number_of_internal_nodes);
        
        int number_of_boundary_nodes=0;
        helper.number_of_nodes=&number_of_boundary_nodes;
        MAP_QUADTREE_MESH<T>::Map_Boundary_Nodes(uniform_grid,cells,&helper,Node_Iterator_Data_Helper);
        boundary_nodes.Resize(number_of_boundary_nodes);
        individual_side_boundary_nodes.Resize(4);
        for(int i=1;i<=4;i++)individual_side_boundary_nodes(i).Resize(number_of_boundary_nodes);
        
        int number_of_ghost_nodes=0;
        helper.number_of_nodes=&number_of_ghost_nodes;
        MAP_QUADTREE_MESH<T>::Map_Ghost_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
        MAP_QUADTREE_MESH<T>::Map_Exterior_Ghost_Cell_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
        ghost_nodes.Resize(number_of_ghost_nodes);
        individual_side_ghost_nodes.Resize(4);
        for(int i=1;i<=4;i++)individual_side_ghost_nodes(i).Resize(number_of_ghost_nodes);
        
        number_of_internal_nodes=0;number_of_boundary_nodes=0;number_of_ghost_nodes=0;
        ARRAY<int> number_of_individual_boundary_nodes(4);ARRAY<int> number_of_individual_ghost_nodes(4);
        
        RANGE<VECTOR<T,2> > interior_domain=Domain();interior_domain.Change_Size((T)-.5*Minimum_Edge_Length());
        RANGE<VECTOR<T,2> > regular_domain=Domain();regular_domain.Change_Size((T).5*Minimum_Edge_Length());
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
                    boundary_nodes(++number_of_boundary_nodes)=i;}
                else{
                    if(node_location.x<regular_domain.min_corner.x) individual_side_ghost_nodes(1)(++number_of_individual_ghost_nodes(1))=i;
                    else if(node_location.x>regular_domain.max_corner.x) individual_side_ghost_nodes(2)(++number_of_individual_ghost_nodes(2))=i;
                    if(node_location.y<regular_domain.min_corner.y) individual_side_ghost_nodes(3)(++number_of_individual_ghost_nodes(3))=i;
                    else if(node_location.y>regular_domain.max_corner.y) individual_side_ghost_nodes(4)(++number_of_individual_ghost_nodes(4))=i;
                    ghost_nodes(++number_of_ghost_nodes)=i;}}}
                    
        for(int i=1;i<=4;i++)individual_side_boundary_nodes(i).Resize(number_of_individual_boundary_nodes(i));
        for(int i=1;i<=4;i++)individual_side_ghost_nodes(i).Resize(number_of_individual_ghost_nodes(i));
        ghost_nodes.Resize(number_of_ghost_nodes);
        node_iterator_data_up_to_date=true;}
}
//#####################################################################
// Function Face_Iterator_Data
//#####################################################################
template<class T> struct FACE_ITERATOR_DATA_HELPER{const QUADTREE_GRID<T>* grid;int* number_of_faces;};
template<class T> static void Face_Iterator_Data_Helper(void* data,const QUADTREE_CELL<T>* cell1,const QUADTREE_CELL<T>* cell2,const int axis)
{
    FACE_ITERATOR_DATA_HELPER<T>* helper=(FACE_ITERATOR_DATA_HELPER<T>*)data;
    const QUADTREE_CELL<T>* cells[]={cell1,cell2};int deepest_cell,deepest_depth;MAP_QUADTREE_MESH<T>::Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    const QUADTREE_CELL<T>* smaller_cell=cells[deepest_cell];
    int face=MAP_QUADTREE_MESH<T>::face_by_axis[axis][deepest_cell];int face_index=smaller_cell->Face(face);(*helper->number_of_faces)++;
    helper->grid->face_iterator_deepest_cells(face_index)=smaller_cell->Cell();helper->grid->faces(face_index)=face;
}
template<class T> void QUADTREE_GRID<T>::
Face_Iterator_Data() const
{
    if(!face_iterator_data_up_to_date){
        face_iterator_deepest_cells.Resize(0);face_iterator_deepest_cells.Resize(number_of_faces);
        faces.Resize(number_of_faces);
        
        FACE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;
        
        int number_of_internal_faces=0;
        helper.number_of_faces=&number_of_internal_faces;
        MAP_QUADTREE_MESH<T>::Map_Faces(uniform_grid,cells,0,&helper,Face_Iterator_Data_Helper);
        internal_faces.Resize(number_of_internal_faces);
        
        int number_of_boundary_faces=0;
        helper.number_of_faces=&number_of_boundary_faces;
        MAP_QUADTREE_MESH<T>::Map_Boundary_Faces(uniform_grid,cells,&helper,Face_Iterator_Data_Helper<T>);
        boundary_faces.Resize(number_of_boundary_faces);
        individual_side_boundary_faces.Resize(4);
        for(int i=1;i<=4;i++)individual_side_boundary_faces(i).Resize(number_of_boundary_faces);
        
        int number_of_ghost_faces=0;
        helper.number_of_faces=&number_of_ghost_faces;
        MAP_QUADTREE_MESH<T>::Map_Ghost_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        MAP_QUADTREE_MESH<T>::Map_Exterior_Ghost_Cell_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        ghost_faces.Resize(number_of_ghost_faces);
        individual_side_ghost_faces.Resize(4);
        for(int i=1;i<=4;i++)individual_side_ghost_faces(i).Resize(number_of_ghost_faces);
        individual_side_domain_ghost_faces.Resize(4);
        for(int i=1;i<=4;i++)individual_side_domain_ghost_faces(i).Resize(number_of_ghost_faces);
        
        number_of_internal_faces=0;number_of_boundary_faces=0;number_of_ghost_faces=0;
        ARRAY<int> number_of_individual_boundary_faces(4);ARRAY<int> number_of_individual_domain_ghost_faces(4);ARRAY<int> number_of_individual_ghost_faces(4);
        
        T tolerance=(T).25*Minimum_Edge_Length();
        RANGE<VECTOR<T,2> > domain=Domain();
        RANGE<VECTOR<T,2> > interior_domain=domain;interior_domain.Change_Size(-tolerance);
        RANGE<VECTOR<T,2> > regular_domain=domain;regular_domain.Change_Size(tolerance);
        Cell_Pointer_From_Index();
        for(int i=1;i<=face_iterator_deepest_cells.m;i++)if(face_iterator_deepest_cells(i)){
            TV face_location=Face_Location(faces(i),cell_pointer_from_index(face_iterator_deepest_cells(i)));
            if(interior_domain.Lazy_Inside(face_location)) internal_faces(++number_of_internal_faces)=i;
            else{
                if(regular_domain.Lazy_Inside(face_location)){
                    if(face_location.x<interior_domain.min_corner.x) individual_side_boundary_faces(1)(++number_of_individual_boundary_faces(1))=i;
                    else if(face_location.x>interior_domain.max_corner.x) individual_side_boundary_faces(2)(++number_of_individual_boundary_faces(2))=i;
                    if(face_location.y<interior_domain.min_corner.y) individual_side_boundary_faces(3)(++number_of_individual_boundary_faces(3))=i;
                    else if(face_location.y>interior_domain.max_corner.y) individual_side_boundary_faces(4)(++number_of_individual_boundary_faces(4))=i;
                    boundary_faces(++number_of_boundary_faces)=i;}
                else{
                    if(face_location.x<regular_domain.min_corner.x) individual_side_ghost_faces(1)(++number_of_individual_ghost_faces(1))=i;
                    else if(face_location.x>regular_domain.max_corner.x) individual_side_ghost_faces(2)(++number_of_individual_ghost_faces(2))=i;
                    if(face_location.y<regular_domain.min_corner.y) individual_side_ghost_faces(3)(++number_of_individual_ghost_faces(3))=i;
                    else if(face_location.y>regular_domain.max_corner.y) individual_side_ghost_faces(4)(++number_of_individual_ghost_faces(4))=i;
                    if(abs(face_location.x-domain.min_corner.x)<tolerance) individual_side_domain_ghost_faces(1)(++number_of_individual_domain_ghost_faces(1))=i;
                    else if(abs(face_location.x-domain.max_corner.x)<tolerance) individual_side_domain_ghost_faces(2)(++number_of_individual_domain_ghost_faces(2))=i;
                    if(abs(face_location.y-domain.min_corner.y)<tolerance) individual_side_domain_ghost_faces(3)(++number_of_individual_domain_ghost_faces(3))=i;
                    else if(abs(face_location.y-domain.max_corner.y)<tolerance) individual_side_domain_ghost_faces(4)(++number_of_individual_domain_ghost_faces(4))=i;
                    ghost_faces(++number_of_ghost_faces)=i;}}}
                    
        for(int i=1;i<=4;i++)individual_side_boundary_faces(i).Resize(number_of_individual_boundary_faces(i));
        for(int i=1;i<=4;i++)individual_side_domain_ghost_faces(i).Resize(number_of_individual_domain_ghost_faces(i));
        for(int i=1;i<=4;i++)individual_side_ghost_faces(i).Resize(number_of_individual_ghost_faces(i));
        ghost_faces.Resize(number_of_ghost_faces);
        face_iterator_data_up_to_date=true;}
}
//#####################################################################
// Function Calculate_Cell_Pointer_From_Index_Array
//#####################################################################
template<class T> static void Calculate_Cell_Pointer_From_Index_Array_Helper(QUADTREE_CELL<T>* cell,ARRAY<QUADTREE_CELL<T>*>& cell_pointer_from_index)
{
    assert(cell_pointer_from_index(cell->Cell())==0);cell_pointer_from_index(cell->Cell())=cell;
    if(cell->Has_Children()) for(int i=0;i<4;i++) Calculate_Cell_Pointer_From_Index_Array_Helper(cell->Child(i),cell_pointer_from_index);
}
template<class T> void QUADTREE_GRID<T>::
Calculate_Cell_Pointer_From_Index_Array(ARRAY<QUADTREE_CELL<T>*>& cell_pointer_from_index) const 
{
    cell_pointer_from_index.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_pointer_from_index,0);
    for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
        for(int j=1-number_of_ghost_cells;j<=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;j++)
            Calculate_Cell_Pointer_From_Index_Array_Helper(cells(i,j),cell_pointer_from_index);
}
//#####################################################################
// Function Calculate_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Neighbors_Array_Helper(void* data,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const int axis)
{
    ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >* neighbors=(ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >*)data;
    const QUADTREE_CELL<T>* cells[2]={c1,c2};
    int smaller_cell=0;int larger_cell=1;if(cells[0]->Depth_Of_This_Cell()<cells[1]->Depth_Of_This_Cell()){smaller_cell=1;larger_cell=0;}
    int min_depth=cells[larger_cell]->Depth_Of_This_Cell();
    while(cells[smaller_cell]->Depth_Of_This_Cell()>min_depth){ // mark the large cell as the neighbor of all of the cells on the path up to the cell of the same level as the large cell
        (*neighbors)(cells[smaller_cell]->Cell())(axis*2+1+(1-smaller_cell))=(QUADTREE_CELL<T>*)cells[larger_cell];
        cells[smaller_cell]=cells[smaller_cell]->Parent();}
    (*neighbors)(cells[0]->Cell())(axis*2+2)=(QUADTREE_CELL<T>*)cells[1]; // mark the two cells that are the same size as neighbors of each other
    (*neighbors)(cells[1]->Cell())(axis*2+1)=(QUADTREE_CELL<T>*)cells[0];
}
template<class T> void QUADTREE_GRID<T>::
Calculate_Neighbors_Array(ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >& neighbors)const
{
    neighbors.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(neighbors.Flattened(),0);
    MAP_QUADTREE_MESH<T>::Map_Faces(uniform_grid,cells,number_of_ghost_cells,&neighbors,Calculate_Neighbors_Array_Helper);
}
//#####################################################################
// Function Calculate_Node_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Node_Neighbors_Array_Helper(void* data,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const int axis)
{
    ARRAY<VECTOR<int,4> >* node_neighbors=(ARRAY<VECTOR<int,4> >*)data;
    const QUADTREE_CELL<T>* cells[2]={c1,c2};
    int deepest_tree,deepest_depth;MAP_QUADTREE_MESH<T>::Get_Deepest_Cell(2,cells,deepest_tree,deepest_depth);
    int node_index1=cells[deepest_tree]->Node(MAP_QUADTREE_MESH<T>::nodes_by_axis[axis][deepest_tree][0]),
        node_index2=cells[deepest_tree]->Node(MAP_QUADTREE_MESH<T>::nodes_by_axis[axis][deepest_tree][1]);
    // add each other as neighbors (in the correct directions)
    int neighbor_position1=2*(1-axis)+1,neighbor_position2=2*(1-axis)+2; // in 2D, since a map face is really a map edge, the axis needs to be reversed
    assert((*node_neighbors)(node_index1)(neighbor_position2)==0);assert((*node_neighbors)(node_index2)(neighbor_position1)==0);
    (*node_neighbors)(node_index1)(neighbor_position2)=node_index2;(*node_neighbors)(node_index2)(neighbor_position1)=node_index1;
}
template<class T> void QUADTREE_GRID<T>::
Calculate_Node_Neighbors_Array(ARRAY<VECTOR<int,4> >& node_neighbors)const
{
    node_neighbors.Resize(number_of_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(node_neighbors.Flattened(),0);
    MAP_QUADTREE_MESH<T>::Map_Faces(uniform_grid,cells,number_of_ghost_cells,&node_neighbors,Calculate_Node_Neighbors_Array_Helper);
}
//#####################################################################
// Function Calculate_Fully_Refined_Block
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Calculate_Fully_Refined_Block_Array(ARRAY<bool>& fully_refined_block) const
{
    fully_refined_block.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(fully_refined_block,false);Neighbors();
    QUADTREE_CELL<T>* cells[number_of_children_per_cell];int lookup[][2]={{0,0},{2,0},{4,0},{4,1}};
    for(DYADIC_GRID_ITERATOR_CELL<QUADTREE_GRID<T> > iterator(*this,3);iterator.Valid();iterator.Next()){
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
template<class T,class T2> struct ENSLAVE_T_JUNCTION_NODES_HELPER{QUADTREE_GRID<T>* quadtree_volume;ARRAY<T2>* nodes;int depth_to_enforce;};
template<class T,class T2> static void Enslave_T_Junction_Nodes_Helper(void* data,const QUADTREE_CELL<T>* c1,const QUADTREE_CELL<T>* c2,const QUADTREE_CELL<T>* c3,const QUADTREE_CELL<T>* c4)
{
    ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2>* helper=(ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2>*)data;
    const QUADTREE_CELL<T>* cells[]={c1,c2,c3,c4};int deepest_cell,maximum_depth,shallowest_cell,minimum_depth;
    MAP_QUADTREE_MESH<T>::Get_Deepest_And_Shallowest_Cell(4,cells,deepest_cell,maximum_depth,shallowest_cell,minimum_depth);
    if(helper->depth_to_enforce>=0&&minimum_depth!=helper->depth_to_enforce) return;
    int current_node=cells[deepest_cell]->Node(3-deepest_cell);
    // check if we are in the corner of the shallowest cell -  if so, we are guaranteed not to be a node that needs to be a linear combination of the corner nodes
    if(cells[shallowest_cell]->Node(3-shallowest_cell)==current_node) return;
    // figure out what our value should be based on our position in the cell
    VECTOR<T,2> cell_location=helper->quadtree_volume->Node_Location(3-deepest_cell,cells[deepest_cell]);
    (*(helper->nodes))(current_node)=LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(*helper->quadtree_volume,cells[shallowest_cell],*helper->nodes,cell_location);
}
template<class T> template<class T2> void QUADTREE_GRID<T>::
Enslave_T_Junction_Nodes(ARRAY<T2>* nodes,int depth_to_enforce)
{
    ENSLAVE_T_JUNCTION_NODES_HELPER<T,T2> helper;helper.nodes=nodes;helper.quadtree_volume=this;helper.depth_to_enforce=depth_to_enforce;
    MAP_QUADTREE_MESH<T>::Map_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Enslave_T_Junction_Nodes_Helper<T,T2>);
}
//#####################################################################
// Function Check_Tree_Consistency
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Check_Tree_Consistency(bool check_cells,bool check_nodes,bool check_faces,bool check_neighbor_structures)
{
    Tree_Topology_Changed();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(check_cells){
        LOG::cout<<"Checking cell indices"<<std::endl;
        ARRAY<QUADTREE_CELL<T>*>& cell_pointer_from_index=Cell_Pointer_From_Index();
        for(int i=1;i<=number_of_cells;i++)for(int j=1;j<=number_of_cells;j++)if(i!=j){
            if(cell_pointer_from_index(i)==cell_pointer_from_index(j)) LOG::cout<<"Different index pointing to the same cell: "<<i<<", "<<j<<std::endl;
            if((cell_pointer_from_index(i)->Center()-cell_pointer_from_index(j)->Center()).Magnitude_Squared()<1E-5) LOG::cout<<"Different cells on top of eachother: "<<i<<", "<<j<<std::endl;}}
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
// Function Move_Contents_Left_Two
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Move_Contents_Left_Two()
{
    int start=1-number_of_ghost_cells,m_end=uniform_grid.numbers_of_cells.x+number_of_ghost_cells,n_end=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;
    ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> > temp_data(start,start+1,start,n_end);ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(-2,0));temp_data.domain.min_corner.x=m_end-1;temp_data.domain.max_corner.x=m_end;temp_data.Calculate_Acceleration_Constants();ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Put(temp_data,cells);
    for(int j=start;j<=n_end;j++){assert(!cells(m_end-1,j)->Has_Children()&&!cells(m_end,j)->Has_Children());
        cells(m_end,j)->Node(1)=cells(m_end-1,j)->Node(0);cells(m_end,j)->Node(3)=cells(m_end-1,j)->Node(2);}
    for(int j=start;j<=n_end;j++){cells(m_end-1,j)->Node(0)=cells(m_end-2,j)->Node(1);cells(m_end-1,j)->Node(2)=cells(m_end-2,j)->Node(3);
        cells(m_end,j)->Face(1)=cells(m_end-1,j)->Face(0);cells(m_end-1,j)->Face(0)=cells(m_end-2,j)->Face(1);}
    for(int j=start;j<=n_end;j+=2)cells(m_end,j)->owner->parents_center=uniform_grid.Node(m_end+2,j+1);
}
//#####################################################################
// Function Move_Contents_Right_Two
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Move_Contents_Right_Two()
{
    int start=1-number_of_ghost_cells,m_end=uniform_grid.numbers_of_cells.x+number_of_ghost_cells,n_end=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;
    ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> > temp_data(m_end-1,m_end,start,n_end);ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(2,0));temp_data.domain.min_corner.x=start;temp_data.domain.max_corner.x=start+1;temp_data.Calculate_Acceleration_Constants();ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Put(temp_data,cells);
    for(int j=start;j<=n_end;j++){assert(!cells(start,j)->Has_Children()&&!cells(start+1,j)->Has_Children());
        cells(start,j)->Node(0)=cells(start+1,j)->Node(1);cells(start,j)->Node(2)=cells(start+1,j)->Node(3);}
    for(int j=start;j<=n_end;j++){cells(start+1,j)->Node(1)=cells(start+2,j)->Node(0);cells(start+1,j)->Node(3)=cells(start+2,j)->Node(2);
        cells(start,j)->Face(0)=cells(start+1,j)->Face(1);cells(start+1,j)->Face(1)=cells(start+2,j)->Face(0);}
    for(int j=start;j<=n_end;j+=2)cells(start,j)->owner->parents_center=uniform_grid.Node(start-1,j+1);
}
//#####################################################################
// Function Move_Contents_Down_Two
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Move_Contents_Down_Two()
{
    int start=1-number_of_ghost_cells,m_end=uniform_grid.numbers_of_cells.x+number_of_ghost_cells,n_end=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;
    ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> > temp_data(start,m_end,start,start+1);ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,-2));temp_data.domain.min_corner.y=n_end-1;temp_data.domain.max_corner.y=n_end;temp_data.Calculate_Acceleration_Constants();ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Put(temp_data,cells);
    for(int i=start;i<=m_end;i++){assert(!cells(i,n_end-1)->Has_Children());assert(!cells(i,n_end)->Has_Children());
        cells(i,n_end)->Node(2)=cells(i,n_end-1)->Node(0);cells(i,n_end)->Node(3)=cells(i,n_end-1)->Node(1);}
    for(int i=start;i<=m_end;i++){cells(i,n_end-1)->Node(0)=cells(i,n_end-2)->Node(2);cells(i,n_end-1)->Node(1)=cells(i,n_end-2)->Node(3);
        cells(i,n_end)->Face(3)=cells(i,n_end-1)->Face(2);cells(i,n_end-1)->Face(2)=cells(i,n_end-2)->Face(3);}
    for(int i=start;i<=m_end;i+=2)cells(i,n_end)->owner->parents_center=uniform_grid.Node(i+1,n_end+2);
}
//#####################################################################
// Function Move_Contents_Up_Two
//#####################################################################
template<class T> void QUADTREE_GRID<T>::
Move_Contents_Up_Two()
{
    int start=1-number_of_ghost_cells,m_end=uniform_grid.numbers_of_cells.x+number_of_ghost_cells,n_end=uniform_grid.numbers_of_cells.y+number_of_ghost_cells;
    ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> > temp_data(start,m_end,n_end-1,n_end);ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Get(temp_data,cells);
    cells.Move_Contents_By_Offset(TV_INT(0,2));temp_data.domain.min_corner.y=start;temp_data.domain.max_corner.y=start+1;temp_data.Calculate_Acceleration_Constants();ARRAY<QUADTREE_CELL<T>*,VECTOR<int,2> >::Put(temp_data,cells); 
    for(int i=start;i<=m_end;i++){assert(!cells(i,start)->Has_Children());assert(!cells(i,start+1)->Has_Children());
        cells(i,start)->Node(0)=cells(i,start+1)->Node(2);cells(i,start)->Node(1)=cells(i,start+1)->Node(3);}
    for(int i=start;i<=m_end;i++){cells(i,start+1)->Node(2)=cells(i,start+2)->Node(0);cells(i,start+1)->Node(3)=cells(i,start+2)->Node(1);
        cells(i,start)->Face(2)=cells(i,start+1)->Face(3);cells(i,start+1)->Face(3)=cells(i,start+2)->Face(2);}
    for(int i=start;i<=m_end;i+=2)cells(i,start)->owner->parents_center=uniform_grid.Node(i+1,start-1);
}
//#####################################################################
template class QUADTREE_GRID<float>;
template void QUADTREE_GRID<float>::Enslave_T_Junction_Nodes(ARRAY<float>* nodes,int depth_to_enforce);
template void QUADTREE_GRID<float>::Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUADTREE_GRID<double>;
template void QUADTREE_GRID<double>::Enslave_T_Junction_Nodes(ARRAY<double>* nodes,int depth_to_enforce);
template void QUADTREE_GRID<double>::Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce);
#endif
#endif
