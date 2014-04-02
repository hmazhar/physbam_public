//#####################################################################
// Copyright 2004, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_3D
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Generate_Topology
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Generate_Topology()
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    topology.Preallocate(5000);
    vertices.Resize(1,m-1,1,n-1,1,mn-1);
    int z=vertices.counts.z,yz=vertices.counts.y*z;
    TV_INT i;
    for(i.x=1;i.x<m;i.x++)for(i.y=2;i.y<n;i.y++)for(i.z=2;i.z<mn;i.z++){ // generate one triangle pair per x edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x+1,i.y,i.z);
        int index=vertices.Standard_Index(i),index_j=index-z,index_ij=index-1,index_j_ij=index_j-1;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_j_ij,index_ij,index,index_j));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_j,index,index_ij,index_j_ij));}
    for(i.x=2;i.x<m;i.x++)for(i.y=1;i.y<n;i.y++)for(i.z=2;i.z<mn;i.z++){ // generate one triangle pair per y edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x,i.y+1,i.z);
        int index=vertices.Standard_Index(i),index_i=index-yz,index_ij=index-1,index_i_ij=index_i-1;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_i_ij,index_i,index,index_ij));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_ij,index,index_i,index_i_ij));}
    for(i.x=2;i.x<m;i.x++)for(i.y=2;i.y<n;i.y++)for(i.z=1;i.z<mn;i.z++){ // generate one triangle pair per z edge
        T phi1=levelset.phi(i),phi2=levelset.phi(i.x,i.y,i.z+1);
        int index=vertices.Standard_Index(i),index_i=index-yz,index_j=index-z,index_i_j=index_i-z;
        if(phi1<=0&&phi2>0)topology.Append(VECTOR<int,4>(index_i_j,index_j,index,index_i));
        else if(phi1>0&&phi2<=0)topology.Append(VECTOR<int,4>(index_i,index,index_j,index_i_j));}
    topology.Compact();
    for(int t=1;t<=topology.m;t++){
        int i,j,k,l;topology(t).Get(i,j,k,l);
        vertices.array(i)=vertices.array(j)=vertices.array(k)=vertices.array(l)=1;}
}
//#####################################################################
// Function Generate_Vertices
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Generate_Vertices()
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    int count=vertices.Sum();
    geometry.Resize(count);
    normals.Resize(count);
    levelset.Compute_Normals();
    int vertex=0;
    for(int i=1;i<m;i++)for(int j=1;j<n;j++)for(int ij=1;ij<mn;ij++)if(vertices(i,j,ij)){ // generate vertices where needed
        vertices(i,j,ij)=++vertex;
        TV position=grid.Center(i,j,ij);
        TV normal=levelset.Normal(position);
        T phi=levelset.Phi(position);
        int iterations=0;
        while(abs(phi)>1e-5*grid.min_dX && (iterations++)<10){
            position-=phi*normal;
            phi=levelset.Phi(position);
            normal=levelset.Normal(position);}
        geometry(vertex)=position;normals(vertex)=normal;}
}
//#####################################################################
// Function Ensure_Vertices_In_Correct_Cells
//#####################################################################
template<class T> void DUALCONTOUR_3D<T>::
Ensure_Vertices_In_Correct_Cells()
{
    int vertex=0;
    TV_INT i;
    for(i.x=1;i.x<grid.counts.x;i.x++)for(i.y=1;i.y<grid.counts.y;i.y++)for(i.z=1;i.z<grid.counts.z;i.z++) if(vertices(i)){
        vertex++;TV_INT v=grid.Clamp_To_Cell(geometry(vertex));
        if(i!=v){
            TV cell_center=grid.Center(i);TV offset((T).5*grid.dX);
            geometry(vertex)=BOX<TV>(cell_center-offset,cell_center+offset).Surface(geometry(vertex));}}
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* DUALCONTOUR_3D<T>::
Get_Triangulated_Surface()
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    particles.array_collection->Add_Elements(geometry.m);
    particles.X=geometry;
    ARRAY<TV>* vertex_normals=new ARRAY<TV>(normals);
    TRIANGLE_MESH& mesh=surface->mesh;
    mesh.number_nodes=geometry.m;
    mesh.elements.Exact_Resize(2*topology.m);
    int current_triangle=1;
    for(int t=1;t<=topology.m;t++){
        int i,j,k,l;topology(t).Get(i,j,k,l);i=vertices.array(i);j=vertices.array(j);k=vertices.array(k);l=vertices.array(l);
        mesh.elements(current_triangle++).Set(i,j,l);mesh.elements(current_triangle++).Set(j,k,l);}
    surface->Update_Triangle_List();surface->Use_Vertex_Normals();surface->vertex_normals=vertex_normals;
    return surface;
}
//#####################################################################
template class DUALCONTOUR_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DUALCONTOUR_3D<double>;
#endif
