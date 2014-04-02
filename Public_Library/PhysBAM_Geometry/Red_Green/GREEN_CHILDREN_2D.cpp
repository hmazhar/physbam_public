//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Red_Green/RED_TRIANGLE.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void GREEN_CHILDREN_2D<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    for(int i=1;i<=cells.m;i++){mapping_array(cells(i))=++cell_count;cells(i)=cell_count;}
}
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void GREEN_CHILDREN_2D<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<3;i++)if(midpoints[i]){
        if(!mapping_array(midpoints[i]))mapping_array(midpoints[i])=++node_count;
        midpoints[i]=mapping_array(midpoints[i]);}
    for(int t=1;t<=cells.m;t++)for(int i=1;i<=3;i++){
        assert(mapping_array(elements(t)(i)));elements(t)(i)=mapping_array(elements(t)(i));}
}
//#####################################################################
// Function Simulate_Green_Refine
//#####################################################################
template<class T> bool GREEN_CHILDREN_2D<T>::
Simulate_Green_Refine(int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid)
{
    bool new_midpoints=false;
    for(int i=0;i<3;i++)if(!midpoints[i]){
        midpoints[i]=parent->Get_Neighbor_Midpoint(i,grid);
        if(midpoints[i])new_midpoints=true;}
    if(new_midpoints){elements.Clean_Memory();cells.Clean_Memory();}else return true;
    int count=(midpoints[0]!=0)+(midpoints[1]!=0)+(midpoints[2]!=0);
    return count<2;
}
//#####################################################################
// Function Finalize_Green_Refinement
//#####################################################################
template<class T> void GREEN_CHILDREN_2D<T>::
Finalize_Green_Refinement(int& number_of_cells)
{
    if(elements.m)return; // old triangles would have been deleted if we need more refinement
    int count=(midpoints[0]!=0)+(midpoints[1]!=0)+(midpoints[2]!=0);assert(count<2);
    if(!count)return;int m=midpoints[0]?0:midpoints[1]?1:2;
    int i=parent->Node((m+2)%3),j=parent->Node(m),k=parent->Node((m+1)%3),jk=midpoints[m];
    elements.Resize(2);elements(1).Set(i,j,jk);elements(2).Set(i,jk,k);
    cells.Resize(2);for(int t=1;t<=2;t++)cells(t)=++number_of_cells;
}
//#####################################################################
// Function Build_Triangle_Mesh
//#####################################################################
template<class T> void GREEN_CHILDREN_2D<T>::
Build_Triangle_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    for(int t=1;t<=cells.m;t++){
        int n[3];elements(t).Get(n[0],n[1],n[2]);
        if(!phi||(*phi)(n[0])<=0||(*phi)(n[1])<=0||(*phi)(n[2])<=0){
            int p[3];
            if(node_to_particle_mapping) for(int i=0;i<3;i++){if(!(*node_to_particle_mapping)(n[i])) (*node_to_particle_mapping)(n[i])=++triangle_mesh.number_nodes;p[i]=(*node_to_particle_mapping)(n[i]);}
            else for(int i=0;i<3;i++) p[i]=n[i];
            triangle_mesh.elements.Append(VECTOR<int,3>(p[0],p[1],p[2]));
            cell_to_triangle_mapping(cells(t))=triangle_mesh.elements.m;}}
}
//#####################################################################
// Function Green_Leaf_Triangle
//#####################################################################
template<class T> int GREEN_CHILDREN_2D<T>::
Green_Leaf_Triangle(const VECTOR<T,2>& location) const
{
    VECTOR<T,2> d=(location-parent->Center()).Rotate_Counterclockwise_Multiple_90(-parent->Orientation());T size=parent->owner->childrens_size;
    if(midpoints[0]){if(d.x<0) return 1;else return 2;}
    else if(midpoints[1]){if(d.x-d.y*3-size<0) return 2;else return 1;}
    else{assert(midpoints[2]);if(d.x+d.y*3+size<0) return 1;else return 2;}
}
//#####################################################################
// Function Add_Ordered_Neighbors
//#####################################################################
template<class T> void GREEN_CHILDREN_2D<T>::
Add_Ordered_Neighbors(ARRAY<ARRAY<int> >& neighbors,ARRAY<ARRAY<int> >& neighbor_links) const
{
    for(int t=1;t<=cells.m;t++){
        int i,j,k;elements(t).Get(i,j,k);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(i),neighbor_links(i),j,k);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(j),neighbor_links(j),k,i);
        TRIANGLE_MESH::Add_Ordered_Neighbors(neighbors(k),neighbor_links(k),i,j);}
}
//#####################################################################
template class GREEN_CHILDREN_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GREEN_CHILDREN_2D<double>;
#endif
#endif
