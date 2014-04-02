//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    mapping_array(cell)=++cell_count;cell=cell_count;
    if(Has_Red_Children()) for(int i=0;i<8;i++) Red_Child(i)->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
    else if(Has_Green_Children())green_children->Create_Cell_Compaction_Mapping(mapping_array,cell_count);
}
//#####################################################################
// Function Red_Refine
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Red_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid)
{
    assert(!red_children);int old_midpoints[6];
    if(green_children){for(int i=0;i<6;i++)old_midpoints[i]=green_children->midpoints[i];delete green_children;green_children=0;}
    else for(int i=0;i<6;i++)old_midpoints[i]=0;
    red_children=new RED_CHILDREN_3D<T>;
    red_children->Initialize_Non_Root_Tetrahedron(this,number_of_cells,number_of_nodes,grid);
    for(int i=0;i<6;i++)if(!old_midpoints[i])Propogate_Refinement(i,number_of_cells,number_of_nodes,grid);
}
//#####################################################################
// Function Propogate_Refinement
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Propogate_Refinement(const int midpoint,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid)
{
    ARRAY<RED_TETRAHEDRON<T>*> neighbors;neighbors.Preallocate(6);int old_number_of_cells;
    do{
        Get_Red_Edge_Neighbors(midpoint,neighbors,grid,false);old_number_of_cells=number_of_cells;
        for(int i=1;i<=neighbors.m;i++){
            RED_TETRAHEDRON<T>* neighbor=neighbors(i);
            if(neighbor->Has_Red_Children()) assert(Depth()==neighbor->Depth() || old_number_of_cells<number_of_cells);
            else if(neighbor->Depth()<Depth()) neighbor->Red_Refine(number_of_cells,number_of_nodes,grid);
            else neighbor->Simulate_Green_Refine(number_of_cells,number_of_nodes,grid);}}
    while(old_number_of_cells<number_of_cells);
}
//#####################################################################
// Function Simulate_Green_Refine
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid)
{
    assert(!red_children);
    if(!Has_Green_Children()){green_children=new GREEN_CHILDREN_3D<T>(this);}
    if(!green_children->Simulate_Green_Refine(number_of_cells,number_of_nodes,grid))
        Red_Refine(number_of_cells,number_of_nodes,grid);
}
//#####################################################################
// Function Finalize_Green_Refinement
//#####################################################################
template<class T> inline void RED_TETRAHEDRON<T>::
Finalize_Green_Refinement(int& number_of_cells)
{
    if(Has_Green_Children())green_children->Finalize_Green_Refinement(number_of_cells);
    else if(Has_Red_Children())for(int i=0;i<8;i++)Red_Child(i)->Finalize_Green_Refinement(number_of_cells);
}
//#####################################################################
// Function Delete_Children
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Delete_Children()
{
    if(Has_Red_Children()){for(int i=0;i<8;i++)Red_Child(i)->Delete_Children();delete red_children;red_children=0;}
    else if(Has_Green_Children()){delete green_children;green_children=0;}
}
//#####################################################################
// Function Get_Red_Neighbor
//#####################################################################
template<class T> RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Get_Red_Neighbor(const int face,const RED_GREEN_GRID_3D<T>& grid,const bool same_depth_only) const
{
    RED_TETRAHEDRON<T>* neighbor=Get_Red_Neighbor_Ancestor(face,grid);if(!neighbor)return 0;
    if(neighbor->Depth()<Depth()){ // go back down as far as we can
        neighbor=neighbor->Red_Descendant(Center(),Depth()); // all half-space tests in Red_Descendant extend out of the tetrahedron to correctly include neighbors 
        if(same_depth_only && neighbor->Depth()!=Depth()) return 0;}
    return neighbor;
}
//#####################################################################
static const VECTOR<int,3> face_jump_table[4]={VECTOR<int,3>(0,-1,-1),VECTOR<int,3>(0,-1,1),VECTOR<int,3>(-1,1,0),VECTOR<int,3>(1,1,0)};
//#####################################################################
// Function Get_Red_Neighbor_Ancestory
//#####################################################################
template<class T> RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Get_Red_Neighbor_Ancestor(const int face,const RED_GREEN_GRID_3D<T>& grid) const
{
    assert(face>=0 && face<4);
    if(!owner->parent){ // we're at the top of the tree, use the uniform grid
        VECTOR<int,3> index=grid.uniform_grid.Closest_Node(Center());
        index+=owner->orientation.Rotate(face_jump_table[face]); // owner->orientation==Orientation() since we're at the top
        return grid.elements.Valid_Index(index)?grid.elements(index):0;}
    static const int valid_child[8][4]={{-1,-1,-1,4},{-1,-1,6,-1},{-1,7,-1,-1},{5,-1,-1,-1},{0,-1,7,5},{4,6,-1,3},{1,-1,5,7},{6,4,-1,2}}; // lookup table in case we are at the right level
    static const int parent_face[8][4]={{0,1,2,-1},{0,1,-1,3},{0,-1,2,3},{-1,1,2,3},{-1,2,-1,-1},{-1,-1,1,-1},{-1,3,-1,-1},{-1,-1,0,-1}}; // transform direction to parent
    if(valid_child[child_index][face]!=-1) return &owner->children[valid_child[child_index][face]];
    else return owner->parent->Get_Red_Neighbor_Ancestor(parent_face[child_index][face],grid); // go further up
}
//#####################################################################
// Function Get_Red_Edge_Neighbors
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Get_Red_Edge_Neighbors(const int edge,ARRAY<RED_TETRAHEDRON<T>*>& neighbors,const RED_GREEN_GRID_3D<T>& grid,const bool same_depth_only)
{
    assert(edge>=0 && edge<6);
    neighbors.Remove_All();
    int depth=Depth();T size=Size();TETRAHEDRAL_GROUP<T> orientation=Orientation();
    int total_steps=edge<2?3:5;
    for(int positivity=0;positivity<2;positivity++){ // search in negative direction and then positive direction if necessary
        static const int first_face[2][6]={{0,2,2,1,3,0},{1,3,0,2,1,3}};
        int face=first_face[positivity][edge];
        static const VECTOR<int,3> edge_direction_table[6]={VECTOR<int,3>(1,0,0),VECTOR<int,3>(0,0,1),VECTOR<int,3>(1,1,-1),VECTOR<int,3>(1,1,1),VECTOR<int,3>(-1,1,1),VECTOR<int,3>(-1,1,-1)};
        VECTOR<int,3> edge_direction=orientation.Rotate(edge_direction_table[edge]);if(!positivity)edge_direction=-edge_direction;
        VECTOR<int,3> face_jump=orientation.Rotate(face_jump_table[face]);
        VECTOR<T,3> edge_center=Midpoint_Location(edge),relative_center=Center()-edge_center;
        RED_TETRAHEDRON<T>* tetrahedron=this;int step;
        for(step=1;step<=total_steps;step++){ // loop around edge in given direction
            assert((VECTOR<int,3>::Dot_Product(edge_direction,face_jump)==0));
            tetrahedron=tetrahedron->Get_Red_Neighbor_Ancestor(face,grid);
            if(!tetrahedron)break;
            if(edge<2){ // edge is axis aligned
                face_jump=VECTOR<int,3>::Cross_Product(edge_direction,face_jump); // rotated 90 degrees around edge since |edge_direction|=1
                relative_center=VECTOR<T,3>::Cross_Product(VECTOR<T,3>(edge_direction),relative_center);} // rotate relative_center as well
            else{ // edge is (+-1,+-1,+-1)
                VECTOR<int,3> face_jump_rotated_90_times_sqrt_3=VECTOR<int,3>::Cross_Product(edge_direction,face_jump); // since |edge_direction|=sqrt(3)
                face_jump=(face_jump+face_jump_rotated_90_times_sqrt_3)/2; // rotated 60 degrees around edge
                VECTOR<T,3> relative_center_rotated_90_times_sqrt_3=VECTOR<T,3>::Cross_Product(VECTOR<T,3>(edge_direction),relative_center); // rotate relative center as well
                relative_center=(relative_center+relative_center_rotated_90_times_sqrt_3)*(T).5;}
            tetrahedron=tetrahedron->Red_Descendant(edge_center+relative_center,depth); // try to go back down to correct depth
            if(!same_depth_only || tetrahedron->Depth()==depth) neighbors.Append(tetrahedron);
            if(tetrahedron->Depth()!=depth){ // check to see if further rotation would remain in this tet
                VECTOR<T,3> next_center=edge_center+relative_center+size*VECTOR<T,3>(face_jump);
                if(!tetrahedron->Outside(next_center)){step+=edge<2?2:3;break;}} // account for additional rotations that would remain in this tet, and switch direction
            face=tetrahedron->Face_From_Face_Jump(face_jump);}
        if(step>total_steps)break;} // if we wrapped all the way around, we don't need to search in other direction
}
//#####################################################################
// Function Get_Neighbor_Midpoint
//#####################################################################
template<class T> int RED_TETRAHEDRON<T>::
Get_Neighbor_Midpoint(const int midpoint,const RED_GREEN_GRID_3D<T>& grid)
{
    ARRAY<RED_TETRAHEDRON<T>*> neighbors;neighbors.Preallocate(6);
    Get_Red_Edge_Neighbors(midpoint,neighbors,grid,true);
    VECTOR<T,3> location=Midpoint_Location(midpoint);
    for(int i=1;i<=neighbors.m;i++){
        RED_TETRAHEDRON<T>* neighbor=neighbors(i);
        if(!(neighbor->Has_Red_Children()||neighbor->Has_Green_Children()))continue;
        int neighbor_midpoint=neighbor->Midpoint_From_Location(location);
        if(neighbor->Has_Red_Children())return neighbor->red_children->Midpoint(neighbor_midpoint);
        int green_midpoint=neighbor->green_children->Midpoint(neighbor_midpoint);if(green_midpoint)return green_midpoint;}
    return 0;
}
//#####################################################################
// Function Build_Tetrahedron_Mesh
//#####################################################################
template<class T> void RED_TETRAHEDRON<T>::
Build_Tetrahedron_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    if(Has_Red_Children())for(int i=0;i<8;i++)Red_Child(i)->Build_Tetrahedron_Mesh(tetrahedron_mesh,phi,cell_to_tetrahedron_mapping,node_to_particle_mapping);
    else if(Has_Green_Children())green_children->Build_Tetrahedron_Mesh(tetrahedron_mesh,phi,cell_to_tetrahedron_mapping,node_to_particle_mapping);
    else if(!phi||(*phi)(Node(0))<=0||(*phi)(Node(1))<=0||(*phi)(Node(2))<=0||(*phi)(Node(3))<=0){
        int p[4];
        if(node_to_particle_mapping) for(int i=0;i<4;i++){
            if(!(*node_to_particle_mapping)(Node(i)))(*node_to_particle_mapping)(Node(i))=++tetrahedron_mesh.number_nodes;
            p[i]=(*node_to_particle_mapping)(Node(i));}
        else for(int i=0;i<4;i++) p[i]=Node(i);
        tetrahedron_mesh.elements.Append(VECTOR<int,4>(p[0],p[1],p[2],p[3]));
        cell_to_tetrahedron_mapping(cell)=tetrahedron_mesh.elements.m;}
}
//#####################################################################
// Function Node_Location
//#####################################################################
template<class T> VECTOR<T,3> RED_TETRAHEDRON<T>::
Node_Location(const int node_index) const
{
    assert(node_index>=0&&node_index<4);
    static const VECTOR<int,3> node_lookup_table[4]={VECTOR<int,3>(-2,-1,0),VECTOR<int,3>(2,-1,0),VECTOR<int,3>(0,1,-2),VECTOR<int,3>(0,1,2)};
    return Center()+Size()*VECTOR<T,3>(Orientation().Rotate(node_lookup_table[node_index]));
}
//#####################################################################
// Function Midpoint_Location
//#####################################################################
template<class T> VECTOR<T,3> RED_TETRAHEDRON<T>::
Midpoint_Location(const int midpoint) const
{
    assert(midpoint>=0&&midpoint<6);
    static const VECTOR<int,3> midpoint_table[6]={VECTOR<int,3>(0,-1,0),VECTOR<int,3>(0,1,0),VECTOR<int,3>(-1,0,-1),VECTOR<int,3>(-1,0,1),VECTOR<int,3>(1,0,1),VECTOR<int,3>(1,0,-1)};
    return Center()+Size()*VECTOR<T,3>(Orientation().Rotate(midpoint_table[midpoint]));
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RED_TETRAHEDRON<T>::
Outside(const VECTOR<T,3>& location,const T thickness_over_2) const
{
    VECTOR<T,3> d=Orientation().Inverse_Rotate(location-Center());
    T threshold=Size()+(T)root_two*thickness_over_2;
    return d.y+d.z<-threshold || d.y-d.z<-threshold || d.y-d.x>threshold || d.y+d.x>threshold;
}
//#####################################################################
// Function Red_Descendant
//#####################################################################
template<class T> const RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Red_Descendant(const VECTOR<T,3>& location,const int depth) const
{
    if(!Has_Red_Children() || Depth()==depth)return this;
    else return Red_Child(Closest_Child(location))->Red_Descendant(location,depth);
}
//#####################################################################
// Function Red_Descendant
//#####################################################################
template<class T> RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Red_Descendant(const VECTOR<T,3>& location,const int depth)
{
    if(!Has_Red_Children() || Depth()==depth)return this;
    else return Red_Child(Closest_Child(location))->Red_Descendant(location,depth);
}
//#####################################################################
// Function Closest_Child
//#####################################################################
template<class T> int RED_TETRAHEDRON<T>::
Closest_Child(const VECTOR<T,3>& location) const
{
    VECTOR<T,3> d=Orientation().Inverse_Rotate(location-Center());T size=Size();
    if(d.x+d.z<0){
        if(d.x-d.z<0)return d.x+d.y+size<0?0:4;
        else return d.z-d.y+size<0?2:7;}
    else{
        if(d.x-d.z<0)return d.z+d.y-size<0?5:3;
        else return d.x-d.y-size<0?6:1;}
}
//#####################################################################
// Function Midpoint_From_Location
//#####################################################################
template<class T> int RED_TETRAHEDRON<T>::
Midpoint_From_Location(const VECTOR<T,3>& location) const
{
    VECTOR<T,3> d=Orientation().Inverse_Rotate(location-Center());int midpoint;
    if(abs(d.y)>(T).5*Size())midpoint=d.y>0?1:0;
    else midpoint=d.x<0?(d.z<0?2:3):(d.z<0?5:4);
    assert((location-Midpoint_Location(midpoint)).Magnitude()<1e-6);
    return midpoint;
}
//#####################################################################
// Function Face_From_Face_Jump
//#####################################################################
template<class T> int RED_TETRAHEDRON<T>::
Face_From_Face_Jump(const VECTOR<int,3>& face_jump) const
{
    VECTOR<int,3> fj=Orientation().Inverse_Rotate(face_jump);int face;
    if(fj.y<0) face=fj.z<0?0:1;
    else face=fj.x<0?2:3;
    assert(fj==face_jump_table[face]);
    return face;
}
//#####################################################################
template class RED_TETRAHEDRON<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_TETRAHEDRON<double>;
#endif
#endif
