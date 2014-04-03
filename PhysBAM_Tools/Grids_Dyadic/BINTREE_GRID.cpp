//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINTREE_GRID
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> BINTREE_GRID<T>::
BINTREE_GRID()
    :uniform_grid(2,0,1),minimum_cell_size(0),number_of_ghost_cells(0),number_of_cells(0),number_of_nodes(0),number_of_faces(0),maximum_depth(1)
{
    Tree_Topology_Changed();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> BINTREE_GRID<T>::
~BINTREE_GRID()
{for(int i=1;i<=allocated_children.m;i++)delete allocated_children(i);allocated_children.Resize(0);}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void BINTREE_GRID<T>::
Initialize(const GRID<TV> uniform_grid_input,const int maximum_depth_input,const int number_of_ghost_cells_input,const bool use_nodes,const bool use_faces)
{
    int i;uniform_grid=uniform_grid_input;Set_Maximum_Depth(maximum_depth_input);number_of_ghost_cells=number_of_ghost_cells_input;
    if(uniform_grid.numbers_of_cells.x%2){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Only even sized grids are supported (grid = "<<uniform_grid.numbers_of_cells<<")"<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
    number_of_cells=number_of_nodes=number_of_faces=0;
    for(i=1;i<=allocated_children.m;i++)delete allocated_children(i);allocated_children.Resize(0);
    cells.Resize(RANGE<TV_INT>(TV_INT()+(1-number_of_ghost_cells),uniform_grid.numbers_of_cells+number_of_ghost_cells));
    ARRAY<int,VECTOR<int,1> > nodes(cells.domain.min_corner.x,cells.domain.max_corner.x+1);
    ARRAY<int,VECTOR<int,1> > faces_u(cells.domain.min_corner.x,cells.domain.max_corner.x+1);
    if(use_nodes) for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++) nodes(i)=++number_of_nodes;
    if(use_faces){
        for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x+1;i++) faces_u(i)=++number_of_faces;}
    ARRAY<int,VECTOR<int,1> > nodes_to_input(0,2);ARRAY<int,VECTOR<int,1> > faces_to_input(0,2);
    for(i=cells.domain.min_corner.x;i<=cells.domain.max_corner.x;i+=2){
        BINTREE_CHILDREN<T>* owner=new BINTREE_CHILDREN<T>(use_nodes,use_faces);
        allocated_children.Append(owner);
        TV center=uniform_grid.Node(i+1),DX=uniform_grid.dX;
        if(use_nodes) for(int ii=i;ii<=i+2;ii++) nodes_to_input(ii-i)=nodes(ii);
        if(use_faces) for(int ii=i;ii<=i+2;ii++) faces_to_input(ii-i)=faces_u(ii);
        owner->Initialize_Pseudo_Root_Cells(number_of_cells,nodes_to_input,faces_to_input,center,DX);
        cells(i)=&owner->children[0];cells(i+1)=&owner->children[1];}
}
//#####################################################################
// Function Tree_Topology_Changed
//#####################################################################
template<class T> void BINTREE_GRID<T>::
Tree_Topology_Changed()
{
    cell_pointer_from_index.Resize(0);cell_pointer_from_index_up_to_date=false;
    neighbors.Resize(0);neighbors_up_to_date=false;
    node_neighbors.Resize(0);node_neighbors_up_to_date=false;
    fully_refined_block.Resize(0);fully_refined_block_up_to_date=false;
    node_iterator_data_up_to_date=false;face_iterator_data_up_to_date=false;
}
//#####################################################################
// Leaf_Cell
//#####################################################################
template<class T> BINTREE_CELL<T>* BINTREE_GRID<T>::
Leaf_Cell(const TV& location,const T thickness) const
{
    TV_INT i=uniform_grid.Cell(location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    if(cells.Domain_Indices().Lazy_Outside(i)) return 0;
    return Inside_Offspring(cells(i),location,thickness);
}
//#####################################################################
// Function Clamped_Leaf_Cell
//#####################################################################
template<class T> BINTREE_CELL<T>* BINTREE_GRID<T>::
Clamped_Leaf_Cell(const TV& location,const T thickness) const
{
    TV clamped_location=uniform_grid.Clamp(location,number_of_ghost_cells);
    TV_INT i=uniform_grid.Cell(clamped_location,number_of_ghost_cells+2); // the +2 is to make sure the behavior is correct when we pass a negative into the cast
    i=TV_INT::Componentwise_Min(i,TV_INT(cells.domain.max_corner.x)); // need this clamp because Cell could potentially return i or ij beyond last index
    return Inside_Offspring(cells(i),clamped_location,thickness);
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> BINTREE_CELL<T>* BINTREE_GRID<T>::
Inside_Offspring(BINTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        cell=cell->Child(location.x>center.x?1:0);}
    return cell;
}
//#####################################################################
// Function Inside_Offspring
//#####################################################################
template<class T> const BINTREE_CELL<T>* BINTREE_GRID<T>::
Inside_Offspring(const BINTREE_CELL<T>* cell,const TV& location,T tolerance) const
{
    assert(Inside_Thickened_Cell(cell,location,tolerance*10));
    while(cell->Has_Children()){
        const TV& center=cell->Center();
        cell=cell->Child(location.x>center.x?1:0);}
    return cell;
}
//#####################################################################
// Function Refine_Cell
//#####################################################################
// refines the given cell so that the node at the given location is touching a refined cell instead of a coarse one
// the location needs to be inside or on one of the faces of the cell
template<class T> void BINTREE_GRID<T>::
Refine_Cell(const int max_depth,BINTREE_CELL<T>* cell,const TV& location,ARRAY<BINTREE_CELL<T>*>* new_cells,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes, 
            ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const BINTREE_GRID<T>* grid)
{
    int depth_of_cell=cell->Depth_Of_This_Cell();assert(depth_of_cell<=maximum_depth);
    if(depth_of_cell>=max_depth) return;
    if(!Inside_Thickened_Cell(cell,location)) return; // make sure that the node is actually in the cell
    if(!cell->Has_Children()) cell->Create_Children(number_of_cells,new_cells,number_of_nodes,new_nodes,number_of_faces,new_faces,grid);
    for(int i=0;i<2;i++) Refine_Cell(max_depth,cell->Child(i),location, new_cells,new_nodes,new_faces,grid); // recur on the children
    Tree_Topology_Changed();
}
//#####################################################################
// Function Refine_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Refine_Cells_Intersecting_Box_Helper(BINTREE_GRID<T>& grid,BINTREE_CELL<T>& cell,const RANGE<VECTOR<T,1> >& box,ARRAY<BINTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    if(cell.Depth_Of_This_Cell()>=refinement_depth || !box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){ // refine
        refined_cells.Append(&cell);cell.Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}
    for(int i=0;i<2;i++) Refine_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,refined_cells,refinement_depth);
}
template<class T> void BINTREE_GRID<T>::
Refine_Cells_Intersecting_Box(const RANGE<VECTOR<T,1> >& box,ARRAY<BINTREE_CELL<T>*>& refined_cells,const int refinement_depth)
{
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,BINTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,BINTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index.x+=1;
    int depth=maximum_depth;if(refinement_depth) depth=refinement_depth;
    for(int i=min_index.x;i<=max_index.x;i++) Refine_Cells_Intersecting_Box_Helper(*this,*cells(i),box,refined_cells,refinement_depth);
    Tree_Topology_Changed();
}
//#####################################################################
// Function Get_Cells_Intersecting_Box
//#####################################################################
template<class T> static void Get_Cells_Intersecting_Box_Helper(const BINTREE_GRID<T>& grid,const BINTREE_CELL<T>& cell,const RANGE<VECTOR<T,1> >& box,ARRAY<BINTREE_CELL<T>*>& intersecting_cells)
{   
    if(!box.Lazy_Intersection(cell.Bounding_Box())) return;
    if(!cell.Has_Children()){
        intersecting_cells.Append(const_cast<BINTREE_CELL<T>*>(&cell));}
    else for(int i=0;i<2;i++) Get_Cells_Intersecting_Box_Helper(grid,*cell.Child(i),box,intersecting_cells);
}   
template<class T> void BINTREE_GRID<T>::
Get_Cells_Intersecting_Box(const RANGE<VECTOR<T,1> >& box,ARRAY<BINTREE_CELL<T>*>& intersecting_cells) const
{   
    TV_INT min_index=INTERPOLATION_UNIFORM<GRID<TV>,BINTREE_CELL<T>*>::Clamped_Index(uniform_grid,cells,box.Minimum_Corner());
    TV_INT max_index=INTERPOLATION_UNIFORM<GRID<TV>,BINTREE_CELL<T>*>::Clamped_Index_End_Minus_One(uniform_grid,cells,box.Maximum_Corner());
    max_index.x+=1;
    for(int i=min_index.x;i<=max_index.x;i++)
        Get_Cells_Intersecting_Box_Helper(*this,*cells(i),box,intersecting_cells);
}
//#####################################################################
// Function Compact_Array_Indices
//#####################################################################
template<class T> void BINTREE_GRID<T>::
Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array,ARRAY<int>* face_mapping_array)
{
    if(cell_mapping_array){cell_mapping_array->Resize(number_of_cells,false);ARRAYS_COMPUTATIONS::Fill(*cell_mapping_array,0);number_of_cells=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
            cells(i)->Create_Cell_Compaction_Mapping(*cell_mapping_array,number_of_cells);}
    if(face_mapping_array){face_mapping_array->Resize(number_of_faces,false);ARRAYS_COMPUTATIONS::Fill(*face_mapping_array,0);number_of_faces=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
                cells(i)->owner->Create_Face_Compaction_Mapping_Helper(*face_mapping_array,number_of_faces);}
    if(node_mapping_array){node_mapping_array->Resize(number_of_nodes,false);ARRAYS_COMPUTATIONS::Fill(*node_mapping_array,0);number_of_nodes=0;
        for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i+=2)
                cells(i)->owner->Create_Node_Compaction_Mapping_Helper(*node_mapping_array,number_of_nodes);}
}
//#####################################################################
// Function Inside_Cell
//#####################################################################
template<class T> bool BINTREE_GRID<T>::
Inside_Cell(const BINTREE_CELL<T>* cell,const TV& location) const
{
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    return location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x;
}
//#####################################################################
// Function Inside_Thickened_Cell
//#####################################################################
template<class T> bool BINTREE_GRID<T>::
Inside_Thickened_Cell(const BINTREE_CELL<T>* cell,const TV& location,const T thickness) const
{
    TV center(cell->Center()),DX_over_two(((T).5+thickness)*cell->DX());
    return location.x>=center.x-DX_over_two.x && location.x<=center.x+DX_over_two.x;
}
//#####################################################################
// Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> const BINTREE_CELL<T>* BINTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(const BINTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<BINTREE_CELL<T>*,2> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Leaf_Cell_By_Neighbor_Path
//#####################################################################
template<class T> BINTREE_CELL<T>* BINTREE_GRID<T>::
Leaf_Cell_By_Neighbor_Path(BINTREE_CELL<T>* starting_cell,const TV& location,const T tolerance,const int max_iterations) const
{
    ARRAY<VECTOR<BINTREE_CELL<T>*,2> >& neighbors=Neighbors();int iterations=0;T tolerance_times_small_cell=tolerance*Minimum_Edge_Length();
    while(starting_cell!=0 && iterations++<max_iterations){
        const TV& center=starting_cell->Center();TV DX_over_two((T).5*starting_cell->DX()+TV(tolerance_times_small_cell));
        if(location.x<center.x-DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(1);
        else if(location.x>center.x+DX_over_two.x) starting_cell=neighbors(starting_cell->Cell())(2);
        else return Inside_Offspring(starting_cell,location,tolerance);} // if all of the above tests passed, then the point is in the cell, so return it
    return Leaf_Cell(location,tolerance);
}
//#####################################################################
// Node_Iterator_Data
//#####################################################################
template<class T> struct NODE_ITERATOR_DATA_HELPER{const BINTREE_GRID<T>* grid;int* number_of_nodes;};
template<class T> static void Node_Iterator_Data_Helper(void* data,const BINTREE_CELL<T>* cell1,const BINTREE_CELL<T>* cell2)
{
    NODE_ITERATOR_DATA_HELPER<T>* helper=(NODE_ITERATOR_DATA_HELPER<T>*)data;
    const BINTREE_CELL<T>* cells[2]={cell1,cell2};int maximum_depth,deepest_cell;MAP_BINTREE_MESH<T>::Get_Deepest_Cell(2,cells,deepest_cell,maximum_depth);
    int node=1-deepest_cell;int node_index=cells[deepest_cell]->Node(node);(*helper->number_of_nodes)++;
    helper->grid->node_iterator_deepest_cells(node_index)=cells[deepest_cell]->Cell();helper->grid->nodes(node_index)=node;
}
template<class T> void BINTREE_GRID<T>::
Node_Iterator_Data() const
{
    if(node_iterator_data_up_to_date) return;
    node_iterator_deepest_cells.Resize(0);node_iterator_deepest_cells.Resize(number_of_nodes);
    nodes.Resize(number_of_nodes);

    NODE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;

    int number_of_internal_nodes=0;
    helper.number_of_nodes=&number_of_internal_nodes;
    MAP_BINTREE_MESH<T>::Map_Nodes(uniform_grid,cells,0,&helper,Node_Iterator_Data_Helper);
    internal_nodes.Resize(number_of_internal_nodes);

    int number_of_boundary_nodes=0;
    helper.number_of_nodes=&number_of_boundary_nodes;
    MAP_BINTREE_MESH<T>::Map_Boundary_Nodes(uniform_grid,cells,&helper,Node_Iterator_Data_Helper);
    boundary_nodes.Resize(number_of_boundary_nodes);
    individual_side_boundary_nodes.Resize(2);
    for(int i=1;i<=2;i++) individual_side_boundary_nodes(i).Resize(number_of_boundary_nodes);

    int number_of_ghost_nodes=0;
    helper.number_of_nodes=&number_of_ghost_nodes;
    MAP_BINTREE_MESH<T>::Map_Ghost_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
    MAP_BINTREE_MESH<T>::Map_Exterior_Ghost_Cell_Nodes(uniform_grid,cells,number_of_ghost_cells,&helper,Node_Iterator_Data_Helper);
    ghost_nodes.Resize(number_of_ghost_nodes);
    individual_side_ghost_nodes.Resize(2);
    for(int i=1;i<=2;i++) individual_side_ghost_nodes(i).Resize(number_of_ghost_nodes);

    number_of_internal_nodes=0;number_of_boundary_nodes=0;number_of_ghost_nodes=0;
    ARRAY<int> number_of_individual_boundary_nodes(2);ARRAY<int> number_of_individual_ghost_nodes(2);

    RANGE<VECTOR<T,1> > interior_domain=Domain();interior_domain.Change_Size((T)-.5*Minimum_Edge_Length());
    RANGE<VECTOR<T,1> > regular_domain=Domain();regular_domain.Change_Size((T).5*Minimum_Edge_Length());
    Cell_Pointer_From_Index();
    for(int i=1;i<=node_iterator_deepest_cells.m;i++) if(node_iterator_deepest_cells(i)){
        TV node_location=Node_Location(nodes(i),cell_pointer_from_index(node_iterator_deepest_cells(i)));
        if(interior_domain.Lazy_Inside(node_location)) internal_nodes(++number_of_internal_nodes)=i;
        else{
            if(regular_domain.Lazy_Inside(node_location)){
                if(node_location.x<interior_domain.min_corner.x) individual_side_boundary_nodes(1)(++number_of_individual_boundary_nodes(1))=i;
                else if(node_location.x>interior_domain.max_corner.x) individual_side_boundary_nodes(2)(++number_of_individual_boundary_nodes(2))=i;
                boundary_nodes(++number_of_boundary_nodes)=i;}
            else{
                if(node_location.x<regular_domain.min_corner.x) individual_side_ghost_nodes(1)(++number_of_individual_ghost_nodes(1))=i;
                else if(node_location.x>regular_domain.max_corner.x) individual_side_ghost_nodes(2)(++number_of_individual_ghost_nodes(2))=i;
                ghost_nodes(++number_of_ghost_nodes)=i;}}}

    for(int i=1;i<=2;i++) individual_side_boundary_nodes(i).Resize(number_of_individual_boundary_nodes(i));
    for(int i=1;i<=2;i++) individual_side_ghost_nodes(i).Resize(number_of_individual_ghost_nodes(i));
    ghost_nodes.Resize(number_of_ghost_nodes);
    node_iterator_data_up_to_date=true;
}
//#####################################################################
// Face_Iterator_Data
//#####################################################################
template<class T> struct FACE_ITERATOR_DATA_HELPER{const BINTREE_GRID<T>* grid;int* number_of_faces;};
template<class T> static void Face_Iterator_Data_Helper(void* data,const BINTREE_CELL<T>* cell1,const BINTREE_CELL<T>* cell2,const int axis)
{
    FACE_ITERATOR_DATA_HELPER<T>* helper=(FACE_ITERATOR_DATA_HELPER<T>*)data;
    const BINTREE_CELL<T>* cells[]={cell1,cell2};int deepest_cell,deepest_depth;MAP_BINTREE_MESH<T>::Get_Deepest_Cell(2,cells,deepest_cell,deepest_depth);
    const BINTREE_CELL<T>* smaller_cell=cells[deepest_cell];
    int face=MAP_BINTREE_MESH<T>::face_by_axis[deepest_cell];int face_index=smaller_cell->Face(face);(*helper->number_of_faces)++;
    helper->grid->face_iterator_deepest_cells(face_index)=smaller_cell->Cell();helper->grid->faces(face_index)=face;
}
template<class T> void BINTREE_GRID<T>::
Face_Iterator_Data() const
{
    if(!face_iterator_data_up_to_date){
        face_iterator_deepest_cells.Resize(0);face_iterator_deepest_cells.Resize(number_of_faces);
        faces.Resize(number_of_faces);
        
        FACE_ITERATOR_DATA_HELPER<T> helper;helper.grid=this;
        
        int number_of_internal_faces=0;
        helper.number_of_faces=&number_of_internal_faces;
        MAP_BINTREE_MESH<T>::Map_Faces(uniform_grid,cells,0,&helper,Face_Iterator_Data_Helper);
        internal_faces.Resize(number_of_internal_faces);
        
        int number_of_boundary_faces=0;
        helper.number_of_faces=&number_of_boundary_faces;
        MAP_BINTREE_MESH<T>::Map_Boundary_Faces(uniform_grid,cells,&helper,Face_Iterator_Data_Helper<T>);
        boundary_faces.Resize(number_of_boundary_faces);
        individual_side_boundary_faces.Resize(2);
        for(int i=1;i<=2;i++)individual_side_boundary_faces(i).Resize(number_of_boundary_faces);
        
        int number_of_ghost_faces=0;
        helper.number_of_faces=&number_of_ghost_faces;
        MAP_BINTREE_MESH<T>::Map_Ghost_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        MAP_BINTREE_MESH<T>::Map_Exterior_Ghost_Cell_Faces(uniform_grid,cells,number_of_ghost_cells,&helper,Face_Iterator_Data_Helper);
        ghost_faces.Resize(number_of_ghost_faces);
        individual_side_ghost_faces.Resize(2);
        for(int i=1;i<=2;i++)individual_side_ghost_faces(i).Resize(number_of_ghost_faces);
        individual_side_domain_ghost_faces.Resize(2);
        for(int i=1;i<=2;i++)individual_side_domain_ghost_faces(i).Resize(number_of_ghost_faces);
        
        number_of_internal_faces=0;number_of_boundary_faces=0;number_of_ghost_faces=0;
        ARRAY<int> number_of_individual_boundary_faces(2);ARRAY<int> number_of_individual_domain_ghost_faces(2);ARRAY<int> number_of_individual_ghost_faces(2);
        
        T tolerance=(T).25*Minimum_Edge_Length();
        RANGE<VECTOR<T,1> > domain=Domain();
        RANGE<VECTOR<T,1> > interior_domain=domain;interior_domain.Change_Size(-tolerance);
        RANGE<VECTOR<T,1> > regular_domain=domain;regular_domain.Change_Size(tolerance);
        Cell_Pointer_From_Index();
        for(int i=1;i<=face_iterator_deepest_cells.m;i++)if(face_iterator_deepest_cells(i)){
            TV face_location=Face_Location(faces(i),cell_pointer_from_index(face_iterator_deepest_cells(i)));
            if(interior_domain.Lazy_Inside(face_location)) internal_faces(++number_of_internal_faces)=i;
            else{
                if(regular_domain.Lazy_Inside(face_location)){
                    if(face_location.x<interior_domain.min_corner.x) individual_side_boundary_faces(1)(++number_of_individual_boundary_faces(1))=i;
                    else if(face_location.x>interior_domain.max_corner.x) individual_side_boundary_faces(2)(++number_of_individual_boundary_faces(2))=i;
                    boundary_faces(++number_of_boundary_faces)=i;}
                else{
                    if(face_location.x<regular_domain.min_corner.x) individual_side_ghost_faces(1)(++number_of_individual_ghost_faces(1))=i;
                    else if(face_location.x>regular_domain.max_corner.x) individual_side_ghost_faces(2)(++number_of_individual_ghost_faces(2))=i;
                    if(abs(face_location.x-domain.min_corner.x)<tolerance) individual_side_domain_ghost_faces(1)(++number_of_individual_domain_ghost_faces(1))=i;
                    else if(abs(face_location.x-domain.max_corner.x)<tolerance) individual_side_domain_ghost_faces(2)(++number_of_individual_domain_ghost_faces(2))=i;
                    ghost_faces(++number_of_ghost_faces)=i;}}}
                    
        for(int i=1;i<=2;i++)individual_side_boundary_faces(i).Resize(number_of_individual_boundary_faces(i));
        for(int i=1;i<=2;i++)individual_side_domain_ghost_faces(i).Resize(number_of_individual_domain_ghost_faces(i));
        for(int i=1;i<=2;i++)individual_side_ghost_faces(i).Resize(number_of_individual_ghost_faces(i));
        ghost_faces.Resize(number_of_ghost_faces);
        face_iterator_data_up_to_date=true;}
}
//#####################################################################
// Calculate_Cell_Pointer_From_Index_Array
//#####################################################################
template<class T> static void Calculate_Cell_Pointer_From_Index_Array_Helper(BINTREE_CELL<T>* cell,ARRAY<BINTREE_CELL<T>*>& cell_pointer_from_index)
{
    assert(cell_pointer_from_index(cell->Cell())==0);cell_pointer_from_index(cell->Cell())=cell;
    if(cell->Has_Children()) for(int i=0;i<2;i++) Calculate_Cell_Pointer_From_Index_Array_Helper(cell->Child(i),cell_pointer_from_index);
}
template<class T> void BINTREE_GRID<T>::
Calculate_Cell_Pointer_From_Index_Array(ARRAY<BINTREE_CELL<T>*>& cell_pointer_from_index) const
{
    cell_pointer_from_index.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_pointer_from_index,(BINTREE_CELL<T>*)0);
    for(int i=1-number_of_ghost_cells;i<=uniform_grid.numbers_of_cells.x+number_of_ghost_cells;i++)
        Calculate_Cell_Pointer_From_Index_Array_Helper(cells(i),cell_pointer_from_index);
}
//#####################################################################
// Calculate_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Neighbors_Array_Helper(void* data,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,const int axis)
{
    ARRAY<VECTOR<BINTREE_CELL<T>*,2> >* neighbors=(ARRAY<VECTOR<BINTREE_CELL<T>*,2> >*)data;
    const BINTREE_CELL<T>* cells[2]={c1,c2};
    int smaller_cell=0;int larger_cell=1;if(cells[0]->Depth_Of_This_Cell()<cells[1]->Depth_Of_This_Cell()){smaller_cell=1;larger_cell=0;}
    (*neighbors)(cells[1]->Cell())(1)=(BINTREE_CELL<T>*)cells[0];
    (*neighbors)(cells[0]->Cell())(2)=(BINTREE_CELL<T>*)cells[1];
    int min_depth=cells[larger_cell]->Depth_Of_This_Cell();
    while(cells[smaller_cell]->Depth_Of_This_Cell()>min_depth){ // mark the large cell as the neighbor of all of the cells on the path up to the cell of the same level as the large cell
        cells[smaller_cell]=cells[smaller_cell]->Parent();
        (*neighbors)(cells[smaller_cell]->Cell())(1+(1-smaller_cell))=(BINTREE_CELL<T>*)cells[larger_cell];}
}
template<class T> void BINTREE_GRID<T>::
Calculate_Neighbors_Array(ARRAY<VECTOR<BINTREE_CELL<T>*,2> >& neighbors) const
{
    neighbors.Resize(number_of_cells,false,false);
    ARRAY_VIEW<BINTREE_CELL<T>*> neighbor_view=neighbors.Flattened();
    ARRAYS_COMPUTATIONS::Fill(neighbor_view,(BINTREE_CELL<T>*)0);
    MAP_BINTREE_MESH<T>::Map_Faces(uniform_grid,cells,number_of_ghost_cells,&neighbors,Calculate_Neighbors_Array_Helper);
}
//#####################################################################
// Calculate_Node_Neighbors_Array
//#####################################################################
template<class T> static void Calculate_Node_Neighbors_Array_Helper(void* data,const BINTREE_CELL<T>* c)
{
    ARRAY<VECTOR<int,2> >* node_neighbors=(ARRAY<VECTOR<int,2> >*)data;
    (*node_neighbors)(c->Node(0))(2) = c->Node(1);
    (*node_neighbors)(c->Node(1))(1) = c->Node(0);
}
template<class T> void BINTREE_GRID<T>::
Calculate_Node_Neighbors_Array(ARRAY<VECTOR<int,2> >& node_neighbors)const
{
    node_neighbors.Resize(number_of_nodes,false,false);
    ARRAY_VIEW<int> neighbor_view=node_neighbors.Flattened();ARRAYS_COMPUTATIONS::Fill(neighbor_view,0);
    MAP_BINTREE_MESH<T>::Map_Cells(uniform_grid,cells,number_of_ghost_cells,&node_neighbors,Calculate_Node_Neighbors_Array_Helper);
}
//#####################################################################
// Calculate_Fully_Refined_Block_Array
//#####################################################################
template<class T> void BINTREE_GRID<T>::
Calculate_Fully_Refined_Block_Array(ARRAY<bool>& fully_refined_block) const
{
    fully_refined_block.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(fully_refined_block,false);Neighbors();
    for(DYADIC_GRID_ITERATOR_CELL<BINTREE_GRID<T> > iterator(*this,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        BINTREE_CELL<T>* cell_pointer=iterator.Cell_Pointer();
        if(cell_pointer->Depth_Of_This_Cell()!=maximum_depth)continue;
        BINTREE_CELL<T>* neighbor=neighbors(cell_pointer->Cell())(1);
        fully_refined_block(iterator.Cell_Index())=(neighbor && neighbor->Depth_Of_This_Cell()==maximum_depth);}  // TODO(jontg): "fully refined" iff the cell and its RIGHT neighbor are refined?  What about the LEFT one?
}
//#####################################################################
// Function Check_Tree_Consistency
//#####################################################################
template<class T> void BINTREE_GRID<T>::
Check_Tree_Consistency(bool check_cells,bool check_nodes,bool check_faces,bool check_neighbor_structures)
{
    Tree_Topology_Changed();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(check_cells){
        LOG::cout<<"Checking cell indices"<<std::endl;
        ARRAY<BINTREE_CELL<T>*>& cell_pointer_from_index=Cell_Pointer_From_Index();
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
template class BINTREE_GRID<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BINTREE_GRID<double>;
#endif
#endif
