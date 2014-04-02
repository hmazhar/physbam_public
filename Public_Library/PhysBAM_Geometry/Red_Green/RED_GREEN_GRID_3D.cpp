//#####################################################################
// Copyright 2004-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Initialize(const GRID<TV>& uniform_grid_input,const int maximum_depth_input)
{
    PHYSBAM_ASSERT(!number_of_cells && !number_of_nodes,"Initialize should never be called twice");
    PHYSBAM_ASSERT(uniform_grid.Is_Isotropic(),"Input grid must be isotropic");
    PHYSBAM_ASSERT((uniform_grid_input.counts.x&3)==1 && (uniform_grid_input.counts.y&3)==1 && (uniform_grid_input.counts.z&3)==1,"Input grid must have multiples of 4 cells in each direction");
    uniform_grid=uniform_grid_input;Set_Maximum_Depth(maximum_depth_input);
    elements.Resize(uniform_grid.Domain_Indices());
    ARRAY<int,VECTOR<int,3> > nodes(uniform_grid.Domain_Indices(2)); // ghost cells to handle elements that poke out of grid
    int m=uniform_grid.counts.x,n=uniform_grid.counts.y,mn=uniform_grid.counts.z;
    for(int i=-1;i<=m+2;i+=2) for(int j=-1;j<=n+2;j+=2) for(int ij=-1;ij<=mn+2;ij+=2)
        if(((i^j)&2)==0 && ((i^ij)&2)==0)nodes(i,j,ij)=++number_of_nodes;
    for(int i=1;i<=m;i+=4) for(int j=3;j<=n;j+=4) for(int ij=3;ij<=mn;ij+=4)Initialize_Tetrahedrons_Across_Primal_Face(1,VECTOR<int,3>(i,j,ij),nodes);
    for(int i=3;i<=m;i+=4) for(int j=1;j<=n;j+=4) for(int ij=3;ij<=mn;ij+=4)Initialize_Tetrahedrons_Across_Primal_Face(2,VECTOR<int,3>(i,j,ij),nodes);
    for(int i=3;i<=m;i+=4) for(int j=3;j<=n;j+=4) for(int ij=1;ij<=mn;ij+=4)Initialize_Tetrahedrons_Across_Primal_Face(3,VECTOR<int,3>(i,j,ij),nodes);
    ARRAY<int> node_mapping_array;Compact_Array_Indices(0,&node_mapping_array); // remove ghost nodes not contained in any tetrahedron
}
//#####################################################################
// Function Initialize_Tetrahedrons_Across_Primal_Face
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Initialize_Tetrahedrons_Across_Primal_Face(const int axis,const VECTOR<int,3>& face,const ARRAY<int,VECTOR<int,3> >& nodes)
{
    TETRAHEDRAL_GROUP<T> g=TETRAHEDRAL_GROUP<T>::Cyclic_Shift_Axes(axis-1); // rotate to push configuration around x axis to given axis
    VECTOR<int,3> dx=g.Rotate(VECTOR<int,3>(1,0,0)),dy=g.Rotate(VECTOR<int,3>(0,1,0)),dz=g.Rotate(VECTOR<int,3>(0,0,1)),dx2=dx*2,dy2=dy*2,dz2=dz*2;
    int node_n=nodes(face-dx2),node_p=nodes(face+dx2),node_nn=nodes(face-dy2-dz2),node_np=nodes(face-dy2+dz2),node_pn=nodes(face+dy2-dz2),node_pp=nodes(face+dy2+dz2); // octahedron vertices
    RED_CHILDREN_3D<T>* owner[4];for(int a=0;a<4;a++)owner[a]=new RED_CHILDREN_3D<T>;
    owner[0]->Initialize_Pseudo_Root_Tetrahedron(number_of_cells,node_n,node_p,node_pn,node_pp,uniform_grid.Node(face+dy),uniform_grid.dX.x,g*TETRAHEDRAL_GROUP<T>::e); elements(face+dy)=&owner[0]->children[0];
    owner[1]->Initialize_Pseudo_Root_Tetrahedron(number_of_cells,node_np,node_pp,node_p,node_n,uniform_grid.Node(face+dz),uniform_grid.dX.x,g*TETRAHEDRAL_GROUP<T>::d2);elements(face+dz)=&owner[1]->children[0];
    owner[2]->Initialize_Pseudo_Root_Tetrahedron(number_of_cells,node_n,node_p,node_np,node_nn,uniform_grid.Node(face-dy),uniform_grid.dX.x,g*TETRAHEDRAL_GROUP<T>::i); elements(face-dy)=&owner[2]->children[0];
    owner[3]->Initialize_Pseudo_Root_Tetrahedron(number_of_cells,node_pn,node_nn,node_p,node_n,uniform_grid.Node(face-dz),uniform_grid.dX.x,g*TETRAHEDRAL_GROUP<T>::b2);elements(face-dz)=&owner[3]->children[0];
}
//#####################################################################
// Function Build_Mesh
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Build_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    tetrahedron_mesh.Clean_Memory();
    cell_to_tetrahedron_mapping.Clean_Memory();cell_to_tetrahedron_mapping.Resize(number_of_cells);
    if(node_to_particle_mapping){node_to_particle_mapping->Clean_Memory();node_to_particle_mapping->Resize(number_of_nodes);}
    for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) for(int ij=1;ij<=uniform_grid.counts.z;ij++) if(elements(i,j,ij))
        elements(i,j,ij)->Build_Tetrahedron_Mesh(tetrahedron_mesh,phi,cell_to_tetrahedron_mapping,node_to_particle_mapping);
    tetrahedron_mesh.number_nodes=node_to_particle_mapping?node_to_particle_mapping->Max():number_of_nodes;
}
//#####################################################################
// Function Calculate_Node_Locations
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Calculate_Node_Locations(ARRAY<VECTOR<T,3> >& node_locations) const
{
    node_locations.Resize(number_of_nodes,false,false);
    for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) for(int ij=1;ij<=uniform_grid.counts.z;ij++) if(elements(i,j,ij)){
        for(int n=0;n<4;n++) node_locations(elements(i,j,ij)->Node(n))=elements(i,j,ij)->Node_Location(n);
        elements(i,j,ij)->Interpolate_Node_Values_To_All_Children(node_locations);}
}
//#####################################################################
// Function Compact_Array_Indices
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array)
{
    if(cell_mapping_array){cell_mapping_array->Resize(number_of_cells,false,false);ARRAYS_COMPUTATIONS::Fill(*cell_mapping_array,0);number_of_cells=0;
        for(int i=1;i<=uniform_grid.counts.x;i++) for(int j=1;j<=uniform_grid.counts.y;j++) for(int ij=1;ij<=uniform_grid.counts.z;ij++) if(elements(i,j,ij))
            elements(i,j,ij)->Create_Cell_Compaction_Mapping(*cell_mapping_array,number_of_cells);}
    if(node_mapping_array){node_mapping_array->Resize(number_of_nodes,false,false);ARRAYS_COMPUTATIONS::Fill(*node_mapping_array,0);number_of_nodes=0;
        for(int i=1;i<=uniform_grid.counts.y;i++) for(int j=1;j<=uniform_grid.counts.y;j++) for(int ij=1;ij<=uniform_grid.counts.z;ij++) if(elements(i,j,ij))
            elements(i,j,ij)->owner->Create_Node_Compaction_Mapping_Helper(*node_mapping_array,number_of_nodes);}
}
//#####################################################################
// Function Red_Leaf_Tetrahedron
//#####################################################################
template<class T> const RED_TETRAHEDRON<T>* RED_GREEN_GRID_3D<T>::
Red_Leaf_Tetrahedron(const VECTOR<T,3>& location) const
{
    VECTOR<T,3> clamped_location=uniform_grid.domain.Thickened((T)-1e-4*uniform_grid.min_dX).Clamp(location);
    return Root_Level_Tetrahedron(clamped_location)->Red_Leaf_Tetrahedron(location);
}
//#####################################################################
// Function Root_Level_Tetrahedron
//#####################################################################
template<class T> const RED_TETRAHEDRON<T>* RED_GREEN_GRID_3D<T>::
Root_Level_Tetrahedron(const VECTOR<T,3>& location) const
{
    if(Outside(location))return 0;
    VECTOR<T,3> v=(T).25*uniform_grid.one_over_dX.x*(location-uniform_grid.domain.Minimum_Corner());
    VECTOR<T,3> shear(v.x+v.y,v.x+v.z,v.y+v.z);
    int i_j=(int)shear.x,i_ij=(int)shear.y,j_ij=(int)shear.z,i=i_j+i_ij-j_ij,j=i_j+j_ij-i_ij,k=i_ij+j_ij-i_j;
    T xy=shear.x-i_j,xz=shear.y-i_ij,yz=shear.z-j_ij;int x_gt_y=xz>yz,x_gt_z=xy>yz,y_gt_z=xy>xz;
    int tet_i=(i<<1)+x_gt_y+x_gt_z+1,tet_j=(j<<1)+2-x_gt_y+y_gt_z,tet_k=(k<<1)+3-x_gt_z-y_gt_z;
    assert(!elements(tet_i,tet_j,tet_k)->Outside(location,uniform_grid.one_over_dX.x*(T)1e-4));
    return elements(tet_i,tet_j,tet_k);
}
//#####################################################################
// Function Create_Overlay_Grid
//#####################################################################
template<class T> void RED_GREEN_GRID_3D<T>::
Create_Overlay_Grid(OCTREE_GRID<T>& overlay_grid,const int number_of_ghost_cells,const bool use_nodes,const bool use_faces) const
{
    ARRAY<VECTOR<T,3> >& node_locations=Node_Locations();
    GRID<TV> overlay_uniform_grid((uniform_grid.counts-1)/2+1,uniform_grid.domain);
    overlay_grid.Initialize(overlay_uniform_grid,maximum_depth,number_of_ghost_cells,use_nodes,use_faces);
    T small_number=(T)1e-4*uniform_grid.domain.Edge_Lengths().Max();
    for(int i=1;i<=number_of_nodes;i++){
        OCTREE_CELL<T>* cell=overlay_grid.Leaf_Cell(node_locations(i));
        while(true){
            VECTOR<T,3> center=cell->Center(),dx_over_two=(T).5*cell->DX();
            if((abs(center.x-dx_over_two.x-node_locations(i).x)<small_number||abs(center.x+dx_over_two.x-node_locations(i).x)<small_number) &&
               (abs(center.y-dx_over_two.y-node_locations(i).y)<small_number||abs(center.y+dx_over_two.y-node_locations(i).y)<small_number) &&
               (abs(center.z-dx_over_two.z-node_locations(i).z)<small_number||abs(center.z+dx_over_two.z-node_locations(i).z)<small_number))
                break;
            else{
                cell->Create_Children(overlay_grid.number_of_cells,0,overlay_grid.number_of_nodes,0,overlay_grid.number_of_faces,0,&overlay_grid);
                int child_to_use=0;if(node_locations(i).x>center.x)child_to_use+=1;if(node_locations(i).y>center.y)child_to_use+=2;if(node_locations(i).z>center.z)child_to_use+=4;
                cell=cell->Child(child_to_use);}}}
}
//#####################################################################
template class RED_GREEN_GRID_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_GREEN_GRID_3D<double>;
#endif
#endif
