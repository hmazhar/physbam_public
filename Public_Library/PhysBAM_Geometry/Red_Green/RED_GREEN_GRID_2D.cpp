//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Initialize(const GRID<TV>& uniform_grid_input,const int maximum_depth_input)
{
    PHYSBAM_ASSERT(!number_of_cells && !number_of_nodes,"Initialize should never be called twice");
    PHYSBAM_ASSERT(uniform_grid_input.Is_Isotropic(),"Input grid must be isotropic");
    PHYSBAM_ASSERT((uniform_grid_input.counts.x&3)==1 && (uniform_grid_input.counts.y&3)==1,"Input grid must have multiples of 4 cells in each direction");
    uniform_grid=uniform_grid_input;T size=uniform_grid.dX.x;Set_Maximum_Depth(maximum_depth_input);
    elements.Resize(uniform_grid.Domain_Indices());
    ARRAY<int,VECTOR<int,2> > nodes(uniform_grid.Domain_Indices());
    for(int i=1;i<=uniform_grid.counts.x;i+=2) for(int j=1;j<=uniform_grid.counts.y;j+=2)
        if(((i^j)&2)==0)nodes(i,j)=++number_of_nodes;
    for(int i=1;i<uniform_grid.counts.x;i+=4) for(int j=1;j<uniform_grid.counts.y;j+=4){
        RED_CHILDREN_2D<T>* owner[4];for(int a=0;a<4;a++)owner[a]=new RED_CHILDREN_2D<T>;
        owner[0]->Initialize_Pseudo_Root_Triangle(number_of_cells,nodes(i+0,j+0),nodes(i+4,j+0),nodes(i+2,j+2),uniform_grid.Node(i+2,j+1),size,0);
        owner[1]->Initialize_Pseudo_Root_Triangle(number_of_cells,nodes(i+4,j+0),nodes(i+4,j+4),nodes(i+2,j+2),uniform_grid.Node(i+3,j+2),size,1);
        owner[2]->Initialize_Pseudo_Root_Triangle(number_of_cells,nodes(i+4,j+4),nodes(i+0,j+4),nodes(i+2,j+2),uniform_grid.Node(i+2,j+3),size,2);
        owner[3]->Initialize_Pseudo_Root_Triangle(number_of_cells,nodes(i+0,j+4),nodes(i+0,j+0),nodes(i+2,j+2),uniform_grid.Node(i+1,j+2),size,3);
        elements(i+2,j+1)=&owner[0]->children[0];elements(i+3,j+2)=&owner[1]->children[0];
        elements(i+2,j+3)=&owner[2]->children[0];elements(i+1,j+2)=&owner[3]->children[0];}
}
//#####################################################################
// Function Build_Mesh
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Build_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    triangle_mesh.Clean_Memory();
    cell_to_triangle_mapping.Clean_Memory();cell_to_triangle_mapping.Resize(number_of_cells);
    if(node_to_particle_mapping){node_to_particle_mapping->Clean_Memory();node_to_particle_mapping->Resize(number_of_nodes);}
    for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) if(elements(i,j))
        elements(i,j)->Build_Triangle_Mesh(triangle_mesh,phi,cell_to_triangle_mapping,node_to_particle_mapping);
    triangle_mesh.number_nodes=node_to_particle_mapping?node_to_particle_mapping->Max():number_of_nodes;
}
//#####################################################################
// Function Calculate_Cell_Pointer_From_Index
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Calculate_Cell_Pointer_From_Index(ARRAY<RED_TRIANGLE<T>*>& cell_pointer_from_index)const
{
    cell_pointer_from_index.Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(cell_pointer_from_index,0);
    for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) if(elements(i,j))
        elements(i,j)->Calculate_Cell_Pointer_From_Index(cell_pointer_from_index);
}
//#####################################################################
// Function Calculate_Node_Locations
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Calculate_Node_Locations(ARRAY<VECTOR<T,2> >& node_locations) const
{
    node_locations.Resize(number_of_nodes);
    for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) if(elements(i,j)){
        for(int n=0;n<3;n++) node_locations(elements(i,j)->Node(n))=elements(i,j)->Node_Location(n);
        elements(i,j)->Interpolate_Node_Values_To_All_Children(node_locations);}
}
//#####################################################################
// Function Calculate_Node_Neighbors
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Calculate_Node_Neighbors(ARRAY<ARRAY<int> >& node_neighbors) const
{
    ARRAY<ARRAY<int> > neighbors(number_of_nodes),neighbor_links(number_of_nodes);
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j)){
        elements(i,j)->Add_Ordered_Neighbors(neighbors,neighbor_links);}
    node_neighbors.Resize(number_of_nodes);
    for(int i=1;i<=number_of_nodes;i++) if(neighbors(i).m){
        node_neighbors(i).Resize(neighbors(i).m);
        ARRAY<bool> not_first(neighbors(i).m);for(int j=1;j<=neighbors(i).m;j++) if(neighbor_links(i)(j)) not_first(neighbor_links(i)(j))=true;
        int node_index=1;while(node_index <= neighbors(i).m && not_first(node_index)) node_index++; // now find the first node in the linked list
        if(node_index > neighbors(i).m) node_index=1; // if we have a cycle (i is in the interior), just use 1
        for(int j=1;j<=node_neighbors(i).m;j++){
            node_neighbors(i)(j)=neighbors(i)(node_index);
            node_index=neighbor_links(i)(node_index);}}
}
//#####################################################################
// Function Compact_Array_Indices
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array)
{
    if(cell_mapping_array){cell_mapping_array->Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(*cell_mapping_array,0);number_of_cells=0;
        for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j))
            elements(i,j)->Create_Cell_Compaction_Mapping(*cell_mapping_array,number_of_cells);}
    if(node_mapping_array){node_mapping_array->Resize(number_of_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(*node_mapping_array,0);number_of_nodes=0;
        for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j))
            elements(i,j)->owner->Create_Node_Compaction_Mapping_Helper(*node_mapping_array,number_of_nodes);}
}
//#####################################################################
// Function Red_Leaf_Triangle
//#####################################################################
template<class T> const RED_TRIANGLE<T>* RED_GREEN_GRID_2D<T>::
Red_Leaf_Triangle(const VECTOR<T,2>& location) const
{
    VECTOR<T,2> clamped_location=uniform_grid.Clamp(location);
    return Root_Level_Triangle(clamped_location)->Red_Leaf_Triangle(clamped_location);
}
//#####################################################################
// Function Root_Level_Triangle
//#####################################################################
template<class T> const RED_TRIANGLE<T>* RED_GREEN_GRID_2D<T>::
Root_Level_Triangle(const VECTOR<T,2>& location) const
{
    VECTOR<int,2> middle=uniform_grid.Clamp_To_Cell(location,0);middle.x=((middle.x-1)&~3)+3;middle.y=((middle.y-1)&~3)+3;
    VECTOR<T,2> d=location-uniform_grid.Node(middle);
    if(d.x+d.y>0){
        if(d.x-d.y>0) return elements(middle.x+1,middle.y);
        else return elements(middle.x,middle.y+1);}
    else{
        if(d.x-d.y>0) return elements(middle.x,middle.y-1);
        else return elements(middle.x-1,middle.y);}
}
//#####################################################################
// Function Create_Overlay_Grid
//#####################################################################
template<class T> void RED_GREEN_GRID_2D<T>::
Create_Overlay_Grid(QUADTREE_GRID<T>& overlay_grid,const int number_of_ghost_cells,const bool use_nodes,const bool use_faces) const
{
    ARRAY<VECTOR<T,2> >& node_locations=Node_Locations();
    GRID<TV> overlay_uniform_grid((uniform_grid.counts-1)/2+1,uniform_grid.domain);
    overlay_grid.Initialize(overlay_uniform_grid,maximum_depth,number_of_ghost_cells,use_nodes,use_faces);
    T small_number=(T)1e-4*uniform_grid.domain.Edge_Lengths().Max();
    for(int i=1;i<=number_of_nodes;i++){
        QUADTREE_CELL<T>* cell=overlay_grid.Leaf_Cell(node_locations(i));
        while(true){
            VECTOR<T,2> center=cell->Center(),dx_over_two=(T).5*cell->DX();
            if((abs(center.x-dx_over_two.x-node_locations(i).x)<small_number||abs(center.x+dx_over_two.x-node_locations(i).x)<small_number) &&
               (abs(center.y-dx_over_two.y-node_locations(i).y)<small_number||abs(center.y+dx_over_two.y-node_locations(i).y)<small_number))
                break;
            else{
                cell->Create_Children(overlay_grid.number_of_cells,0,overlay_grid.number_of_nodes,0,overlay_grid.number_of_faces,0,&overlay_grid);
                int child_to_use=0;if(node_locations(i).x>center.x)child_to_use+=1;if(node_locations(i).y>center.y)child_to_use+=2;
                cell=cell->Child(child_to_use);}}}
}
//#####################################################################
template class RED_GREEN_GRID_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_GREEN_GRID_2D<double>;
#endif
#endif
