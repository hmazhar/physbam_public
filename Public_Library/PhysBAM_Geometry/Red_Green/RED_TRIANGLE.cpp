//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_TRIANGLE.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    mapping_array(cell)=++cell_count;cell=cell_count;
    if(Has_Red_Children()) for(int i=0;i<4;i++) Red_Child(i)->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
    else if(Has_Green_Children())green_children->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
}
//#####################################################################
// Function Red_Refine
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Red_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid)
{
    assert(!red_children);
    int old_midpoints[3];
    if(green_children){for(int i=0;i<3;i++)old_midpoints[i]=green_children->midpoints[i];delete green_children;green_children=0;}
    else for(int i=0;i<3;i++)old_midpoints[i]=false;
    red_children=new RED_CHILDREN_2D<T>;
    red_children->Initialize_Non_Root_Triangle(this,number_of_cells,number_of_nodes,grid);
    for(int i=0;i<3;i++)if(!old_midpoints[i])Propogate_Refinement(i,number_of_cells,number_of_nodes,grid);
}
//#####################################################################
// Function Propogate_Refinement
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Propogate_Refinement(const int midpoint,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid)
{
    for(;;){
        RED_TRIANGLE<T>* neighbor=Get_Red_Neighbor(midpoint,grid,false);if(!neighbor)return;
        if(neighbor->Depth()<Depth())neighbor->Red_Refine(number_of_cells,number_of_nodes,grid);
        else if(neighbor->Has_Red_Children())return;
        else{neighbor->Simulate_Green_Refine(number_of_cells,number_of_nodes,grid);return;}}
}
//#####################################################################
// Function Simulate_Green_Refine
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid)
{
    assert(!red_children);
    if(!Has_Green_Children())green_children=new GREEN_CHILDREN_2D<T>(this);
    if(!green_children->Simulate_Green_Refine(number_of_nodes,grid))
        Red_Refine(number_of_cells,number_of_nodes,grid);
}
//#####################################################################
// Function Finalize_Green_Refinement
//#####################################################################
template<class T> inline void RED_TRIANGLE<T>::
Finalize_Green_Refinement(int& number_of_cells)
{
    if(Has_Green_Children())green_children->Finalize_Green_Refinement(number_of_cells);
    else if(Has_Red_Children())for(int i=0;i<4;i++)Red_Child(i)->Finalize_Green_Refinement(number_of_cells);
}
//#####################################################################
// Function Delete_Children
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Delete_Children()
{
    if(Has_Red_Children()){for(int i=0;i<4;i++)Red_Child(i)->Delete_Children();delete red_children;red_children=0;}
    else if(Has_Green_Children()){delete green_children;green_children=0;}
}
//#####################################################################
// Function Get_Red_Neighbor
//#####################################################################
template<class T> RED_TRIANGLE<T>* RED_TRIANGLE<T>::
Get_Red_Neighbor(const int edge,const RED_GREEN_GRID_2D<T>& grid,const bool same_depth_only) const
{
    assert(edge>=0 && edge<=2);
    if(!owner->parent){ // we're at the top of the tree, use the uniform grid
        VECTOR<int,2> index=grid.uniform_grid.Closest_Node(Center());
        static const VECTOR<int,2> edge_jump[3]={VECTOR<int,2>(0,-2),VECTOR<int,2>(1,1),VECTOR<int,2>(-1,1)};
        index+=edge_jump[edge].Rotate_Counterclockwise_Multiple_90(owner->orientation); // owner->orientation==Orientation() since we're at the top
        return grid.elements.Valid_Index(index)?grid.elements(index):0;}
    static const int valid_child[4][3]={{-1,3,-1},{-1,-1,2},{3,1,-1},{2,-1,0}}; // lookup table in case we are at the right level
    if(valid_child[child_index][edge]!=-1) return &owner->children[valid_child[child_index][edge]];
    else{ // we need to keep going up the tree to find the node where we can go down the right branch
        static const int parent_direction[4][3]={{0,-1,2},{0,1,-1},{-1,-1,1},{-1,2,-1}}; // transform direction to parent
        RED_TRIANGLE<T>* correct_branch=owner->parent->Get_Red_Neighbor(parent_direction[child_index][edge],grid,same_depth_only);if(!correct_branch)return 0;
        if(!correct_branch->Has_Red_Children())return same_depth_only?0:correct_branch;
        static const int correct_child[4][3]={{1,-1,1},{0,0,-1},{-1,-1,3},{-1,2,-1}}; // lookup table to go back down the tree
        return correct_branch->Red_Child(correct_child[child_index][edge]);}
}
//#####################################################################
// Function Get_Neighbor_Midpoint
//#####################################################################
template<class T> int RED_TRIANGLE<T>::
Get_Neighbor_Midpoint(const int midpoint,const RED_GREEN_GRID_2D<T>& grid) const
{
    static const int reverse_direction[3]={0,2,1};
    RED_TRIANGLE<T>* neighbor=Get_Red_Neighbor(midpoint,grid);
    if(!neighbor)return 0;
    else if(neighbor->Has_Green_Children())return neighbor->green_children->Midpoint(reverse_direction[midpoint]);
    else if(neighbor->Has_Red_Children())return neighbor->red_children->Midpoint(reverse_direction[midpoint]);
    else return 0;
}
//#####################################################################
// Function Build_Triangle_Mesh
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Build_Triangle_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    if(Has_Red_Children())for(int i=0;i<4;i++)Red_Child(i)->Build_Triangle_Mesh(triangle_mesh,phi,cell_to_triangle_mapping,node_to_particle_mapping);
    else if(Has_Green_Children())green_children->Build_Triangle_Mesh(triangle_mesh,phi,cell_to_triangle_mapping,node_to_particle_mapping);
    else if(!phi||(*phi)(Node(0))<=0||(*phi)(Node(1))<=0||(*phi)(Node(2))<=0){
        int p[3];
        if(node_to_particle_mapping) for(int i=0;i<3;i++){if(!(*node_to_particle_mapping)(Node(i))) (*node_to_particle_mapping)(Node(i))=++triangle_mesh.number_nodes;p[i]=(*node_to_particle_mapping)(Node(i));}
        else for(int i=0;i<3;i++) p[i]=Node(i);
        triangle_mesh.elements.Append(VECTOR<int,3>(p[0],p[1],p[2]));
        cell_to_triangle_mapping(cell)=triangle_mesh.elements.m;}
}
//#####################################################################
// Function Node_Location
//#####################################################################
template<class T> VECTOR<T,2> RED_TRIANGLE<T>::
Node_Location(int node_index) const
{
    assert(node_index>=0&&node_index<3);
    static const VECTOR<T,2> node_lookup_table[3]={VECTOR<T,2>(-2,-1),VECTOR<T,2>(2,-1),VECTOR<T,2>(0,1)};
    return Center()+owner->childrens_size*node_lookup_table[node_index].Rotate_Counterclockwise_Multiple_90(Orientation());
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RED_TRIANGLE<T>::
Outside(const VECTOR<T,2>& location,const T thickness_over_2) const
{
    VECTOR<T,2> x[3];
    for(int i=0;i<3;i++)x[i]=Node_Location(i);
    return TRIANGLE_2D<T>::Outside(location,x[0],x[1],x[2],thickness_over_2);
}
//#####################################################################
// Function Red_Leaf_Triangle
//#####################################################################
template<class T> const RED_TRIANGLE<T>* RED_TRIANGLE<T>::
Red_Leaf_Triangle(const VECTOR<T,2>& location) const
{
    if(!Has_Red_Children())return this;
    VECTOR<T,2> d=(location-Center()).Rotate_Counterclockwise_Multiple_90(-Orientation());
    T size=owner->childrens_size;
    if(d.x<0){
        if(d.x+d.y+size>0)return Red_Child(3)->Red_Leaf_Triangle(location);
        else return Red_Child(0)->Red_Leaf_Triangle(location);}
    else{
        if(d.x-d.y-size>0)return Red_Child(1)->Red_Leaf_Triangle(location);
        else return Red_Child(2)->Red_Leaf_Triangle(location);}
}
//#####################################################################
// Function Calculate_Cell_Pointer_From_Index 
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Calculate_Cell_Pointer_From_Index(ARRAY<RED_TRIANGLE<T>*>& cell_pointer_from_index)
{
    cell_pointer_from_index(cell)=this;
    if(Has_Red_Children())for(int i=0;i<4;i++)Red_Child(i)->Calculate_Cell_Pointer_From_Index(cell_pointer_from_index);
}
//#####################################################################
// Function Add_Ordered_Neighbors
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Add_Ordered_Neighbors(ARRAY<ARRAY<int> >& neighbors,ARRAY<ARRAY<int> >& neighbor_links) const
{
    if(Has_Red_Children())for(int i=0;i<4;i++)Red_Child(i)->Add_Ordered_Neighbors(neighbors,neighbor_links);
    else if(Has_Green_Children())green_children->Add_Ordered_Neighbors(neighbors,neighbor_links);
    else{
        int i=Node(0),j=Node(1),k=Node(2);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(i),neighbor_links(i),j,k);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(j),neighbor_links(j),k,i);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(k),neighbor_links(k),i,j);}
}
//#####################################################################
// Function Create_Overlay_Grid
//#####################################################################
template<class T> void RED_TRIANGLE<T>::
Create_Overlay_Grid(QUADTREE_GRID<T>& overlay_grid,QUADTREE_CELL<T>* cell1,QUADTREE_CELL<T>* cell2) const
{
    int orientation=Orientation();QUADTREE_CELL<T>* cells[2]={cell1,cell2};
    if(Has_Red_Children()){
        if(!cell1->Has_Children())cell1->Create_Children(overlay_grid.number_of_cells,0,overlay_grid.number_of_nodes,0,overlay_grid.number_of_faces,0,&overlay_grid);
        if(!cell2->Has_Children())cell2->Create_Children(overlay_grid.number_of_cells,0,overlay_grid.number_of_nodes,0,overlay_grid.number_of_faces,0,&overlay_grid);
        for(int i=0;i<4;i++){
            static const int recursion_cell[4]={0,1,1,0};
            static const int recursion_children[4][4][2]={{{0,1},{0,1},{2,0},{1,3}},{{1,3},{1,3},{0,1},{3,2}},{{3,2},{3,2},{1,3},{2,0}},{{2,0},{2,0},{3,2},{0,1}}};
            int child1=recursion_children[orientation][i][0],child2=recursion_children[orientation][i][1];
            Red_Child(i)->Create_Overlay_Grid(overlay_grid,cells[recursion_cell[i]]->Child(child1),cells[recursion_cell[i]]->Child(child2));}}
}
//#####################################################################
template class RED_TRIANGLE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_TRIANGLE<double>;
#endif
#endif
