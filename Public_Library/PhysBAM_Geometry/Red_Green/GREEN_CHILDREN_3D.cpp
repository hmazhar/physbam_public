//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Cell_Compaction_Mapping
//#####################################################################
template<class T> void GREEN_CHILDREN_3D<T>::
Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count)
{
    for(int i=1;i<=cells.m;i++){mapping_array(cells(i))=++cell_count;cells(i)=cell_count;}
}
//#####################################################################
// Function Create_Node_Compaction_Mapping_Helper
//#####################################################################
template<class T> void GREEN_CHILDREN_3D<T>::
Create_Node_Compaction_Mapping_Helper(ARRAY<int>& mapping_array,int& node_count)
{
    for(int i=0;i<6;i++)if(midpoints[i]){
        if(!mapping_array(midpoints[i]))mapping_array(midpoints[i])=++node_count;
        midpoints[i]=mapping_array(midpoints[i]);}
    for(int t=1;t<=cells.m;t++)for(int i=1;i<=4;i++){
        assert(mapping_array(elements(t)(i)));elements(t)(i)=mapping_array(elements(t)(i));}
}
//#####################################################################
// Function Simulate_Green_Refine
//#####################################################################
template<class T> bool GREEN_CHILDREN_3D<T>::
Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid)
{
    bool new_midpoints=false;
    for(int i=0;i<6;i++)if(!midpoints[i]){
        midpoints[i]=parent->Get_Neighbor_Midpoint(i,grid);
        if(midpoints[i])new_midpoints=true;}
    if(new_midpoints){elements.Clean_Memory();cells.Clean_Memory();}else return true;
    int ij=midpoints[0],kl=midpoints[1],ik=midpoints[2],il=midpoints[3],jl=midpoints[4],jk=midpoints[5];
    int count=(ij!=0)+(ik!=0)+(il!=0)+(jk!=0)+(jl!=0)+(kl!=0);
    switch(count){
        case 0:PHYSBAM_FATAL_ERROR();
        case 1:return true;
        case 2:int new_midpoint;
               if     ((ik&&jk) || (il&&jl))new_midpoint=0;else if((ij&&jk) || (il&&kl)) new_midpoint=2;
               else if((ij&&jl) || (ik&&kl))new_midpoint=3;else if((ij&&ik) || (jl&&kl)) new_midpoint=5;
               else if((ij&&il) || (jk&&kl))new_midpoint=4;else if((ik&&il) || (jk&&jl)) new_midpoint=1;
               else return true;
               midpoints[new_midpoint]=++number_of_nodes;parent->Propogate_Refinement(new_midpoint,number_of_cells,number_of_nodes,grid);return true;
        case 3:if((ij&&ik&&jk) || (ij&&il&&jl) || (ik&&il&&kl) || (jk&&jl&&kl)) return true;
        default:return false;}
}
//#####################################################################
// Function Finalize_Green_Refinement
//#####################################################################
template<class T> void GREEN_CHILDREN_3D<T>::
Finalize_Green_Refinement(int& number_of_cells)
{
    if(elements.m)return; // old tetrahedrons would have been deleted if we need more refinement
    int i=parent->Node(0),j=parent->Node(1),k=parent->Node(2),l=parent->Node(3);
    int ij=midpoints[0],kl=midpoints[1],ik=midpoints[2],il=midpoints[3],jl=midpoints[4],jk=midpoints[5];
    int count=(ij!=0)+(ik!=0)+(il!=0)+(jk!=0)+(jl!=0)+(kl!=0);
    switch(count){
        case 1:{int a,b,c,d,ab;
               if     (ij){a=i;b=j;c=k;d=l;ab=ij;}else if(ik){a=i;b=k;c=l;d=j;ab=ik;}
               else if(il){a=i;b=l;c=j;d=k;ab=il;}else if(jk){a=j;b=k;c=i;d=l;ab=jk;}
               else if(jl){a=j;b=l;c=k;d=i;ab=jl;}else       {a=k;b=l;c=i;d=j;ab=kl;assert(kl);}
               elements.Resize(2);elements(1).Set(b,d,c,ab);elements(2).Set(a,c,d,ab);break;}
        case 2:{int a,b,c,d,ab,cd;
               if     (ij&&kl){ab=ij;cd=kl;a=i;b=j;c=k;d=l;}
               else if(ik&&jl){ab=ik;cd=jl;a=i;b=k;c=l;d=j;}
               else           {ab=il;cd=jk;a=i;b=l;c=j;d=k;assert(il&&jk);}
               elements.Resize(4);elements(1).Set(d,b,ab,cd);elements(2).Set(b,c,ab,cd);elements(3).Set(a,d,ab,cd);elements(4).Set(c,a,ab,cd);break;}
        case 3:{int a,b,c,d,ab,ac,bc;
               if     (ij&&ik&&jk){a=i;b=j;c=k;d=l;ab=ij;ac=ik;bc=jk;}
               else if(ij&&il&&jl){a=j;b=i;c=l;d=k;ab=ij;ac=jl;bc=il;}
               else if(ik&&il&&kl){a=i;b=k;c=l;d=j;ab=ik;ac=il;bc=kl;}
               else               {a=j;b=l;c=k;d=i;ab=jl;ac=jk;bc=kl;assert(jk&&jl&&kl);}
               elements.Resize(4);elements(1).Set(c,d,ac,bc);elements(2).Set(b,d,bc,ab);elements(3).Set(a,d,ab,ac);elements(4).Set(d,ab,ac,bc);break;}
        default:PHYSBAM_FATAL_ERROR();}
    cells.Resize(elements.m);for(int t=1;t<=cells.m;t++)cells(t)=++number_of_cells;
}
//#####################################################################
// Function Build_Tetrahedron_Mesh
//#####################################################################
template<class T> void GREEN_CHILDREN_3D<T>::
Build_Tetrahedron_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const
{
    for(int t=1;t<=cells.m;t++){
         int n[4];elements(t).Get(n[0],n[1],n[2],n[3]);
        if(!phi||(*phi)(n[0])<=0||(*phi)(n[1])<=0||(*phi)(n[2])<=0||(*phi)(n[3])<=0){
            int p[4];
            if(node_to_particle_mapping) for(int i=0;i<4;i++){if(!(*node_to_particle_mapping)(n[i])) (*node_to_particle_mapping)(n[i])=++tetrahedron_mesh.number_nodes;p[i]=(*node_to_particle_mapping)(n[i]);}
            else for(int i=0;i<4;i++) p[i]=n[i];
            tetrahedron_mesh.elements.Append(VECTOR<int,4>(p[0],p[1],p[2],p[3]));
            cell_to_tetrahedron_mapping(cells(t))=tetrahedron_mesh.elements.m;}}
}
//#####################################################################
// Function Green_Leaf_Tetrahedron
//#####################################################################
template<class T> int GREEN_CHILDREN_3D<T>::
Green_Leaf_Tetrahedron(const VECTOR<T,3>& location,const ARRAY<VECTOR<T,3> >& node_locations) const
{
    T max_inside=-FLT_MAX;int result=0;
    for(int t=1;t<=elements.m;t++){
        int i,j,k,l;elements(t).Get(i,j,k,l);
        VECTOR<T,3> w=TETRAHEDRON<T>::Barycentric_Coordinates(location,node_locations(i),node_locations(j),node_locations(k),node_locations(l));
        T inside=min(w.x,w.y,w.z,1-w.x-w.y-w.z);if(max_inside<inside){max_inside=inside;result=t;}}
    return result;
}
//#####################################################################
template class GREEN_CHILDREN_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GREEN_CHILDREN_3D<double>;
#endif
#endif
