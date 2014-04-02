//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_REFRESH__
#define __TETRAHEDRALIZED_VOLUME_REFRESH__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;
template<class TV> class GRID;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Update_Tetrahedron_List
//#####################################################################
template<class T>
void Update_Tetrahedron_List(TETRAHEDRALIZED_VOLUME<T>& tv)
{
    if(!tv.tetrahedron_list) tv.tetrahedron_list=new ARRAY<TETRAHEDRON<T> >;
    tv.tetrahedron_list->Resize(tv.mesh.elements.m,false,false);
    for(int t=1;t<=tv.mesh.elements.m;t++){
        (*tv.tetrahedron_list)(t).X=tv.particles.X.Subset(tv.mesh.elements(t));
        (*tv.tetrahedron_list)(t).Create_Triangles();}
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T>
void Initialize_Hierarchy(TETRAHEDRALIZED_VOLUME<T>& tv,const bool update_boxes) // creates and updates the boxes as well
{
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    delete tv.hierarchy;
    if(tv.tetrahedron_list) tv.hierarchy=new TETRAHEDRON_HIERARCHY<T>(tv.mesh,tv.particles,*tv.tetrahedron_list,update_boxes);
    else tv.hierarchy=new TETRAHEDRON_HIERARCHY<T>(tv.mesh,tv.particles,update_boxes);
}
//#####################################################################
// Funcion Initialize_Octahedron_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Octahedron_Mesh_And_Particles(TETRAHEDRALIZED_VOLUME<T>& tv,const GRID<VECTOR<T,3> >& grid)
{
    int i,j,k,m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    tv.particles.array_collection->Delete_All_Elements();
    tv.mesh.Initialize_Octahedron_Mesh(m,n,mn);
    tv.particles.array_collection->Add_Elements(m*n*mn+(m+1)*(n+1)*(mn+1));
    for(k=1;k<=mn;k++) for(j=1;j<=n;j++) for(i=1;i<=m;i++) tv.particles.X(++particle)=grid.X(i,j,k);
    for(k=0;k<=mn;k++) for(j=0;j<=n;j++) for(i=0;i<=m;i++) tv.particles.X(++particle)=grid.X(i,j,k)+(T).5*grid.dX;
}
//#####################################################################
// Funcion Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Cube_Mesh_And_Particles(TETRAHEDRALIZED_VOLUME<T>& tv,const GRID<VECTOR<T,3> >& grid)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    tv.particles.array_collection->Delete_All_Elements();
    tv.mesh.Initialize_Cube_Mesh(m,n,mn);
    tv.particles.array_collection->Add_Elements(m*n*mn);
    for(int k=1;k<=mn;k++) for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) tv.particles.X(++particle)=grid.X(i,j,k);
}
//#####################################################################
// Function Initialize_Prismatic_Cube_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Prismatic_Cube_Mesh_And_Particles(TETRAHEDRALIZED_VOLUME<T>& tv,const GRID<VECTOR<T,3> >& grid)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z,particle=0;
    tv.particles.array_collection->Delete_All_Elements();
    tv.mesh.Initialize_Prismatic_Cube_Mesh(m,n,mn);
    tv.particles.array_collection->Add_Elements(m*n*mn);
    for(int k=1;k<=mn;k++) for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) tv.particles.X(++particle)=grid.X(i,j,k);
}
//#####################################################################
// Function Initialize_Triangulated_Surface
//#####################################################################
template<class T>
void Initialize_Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T>& tv)
{
    if(tv.mesh.number_nodes!=tv.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    delete tv.triangulated_surface;
    if(!tv.mesh.boundary_mesh) tv.mesh.Initialize_Boundary_Mesh();
    tv.triangulated_surface=new TRIANGULATED_SURFACE<T>(*tv.mesh.boundary_mesh,tv.particles);
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T>
void Rescale(TETRAHEDRALIZED_VOLUME<T>& tv,const T scaling_x,const T scaling_y,const T scaling_z)
{
    typedef VECTOR<T,3> TV;
    for(int k=1;k<=tv.particles.array_collection->Size();k++) tv.particles.X(k)*=TV(scaling_x,scaling_y,scaling_z);
}
//#####################################################################
// Function Compute_Tetrahedron_Volumes
//#####################################################################
template<class T>
void Compute_Tetrahedron_Volumes(TETRAHEDRALIZED_VOLUME<T>& tv)
{
    if(!tv.tetrahedron_volumes) tv.tetrahedron_volumes=new ARRAY<T>;
    tv.tetrahedron_volumes->Resize(tv.mesh.elements.m,false,false);
    for(int t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
        (*tv.tetrahedron_volumes)(t)=TETRAHEDRON<T>::Signed_Volume(tv.particles.X(i),tv.particles.X(j),tv.particles.X(k),tv.particles.X(l));}
}
//#####################################################################
// Function Compute_Nodal_Volumes
//#####################################################################
template<class T>
void Compute_Nodal_Volumes(TETRAHEDRALIZED_VOLUME<T>& tv,bool save_tetrahedron_volumes)
{
    if(!tv.nodal_volumes) tv.nodal_volumes=new ARRAY<T>;
    *tv.nodal_volumes=CONSTANT_ARRAY<T>(tv.particles.array_collection->Size(),0);
    if(save_tetrahedron_volumes){
        if(!tv.tetrahedron_volumes) tv.tetrahedron_volumes=new ARRAY<T>;
        tv.tetrahedron_volumes->Resize(tv.mesh.elements.m);
        for(int t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
            (*tv.tetrahedron_volumes)(t)=TETRAHEDRON<T>::Signed_Volume(tv.particles.X(i),tv.particles.X(j),tv.particles.X(k),tv.particles.X(l));
            T volume=(T).25*(*tv.tetrahedron_volumes)(t);
            (*tv.nodal_volumes)(i)+=volume;(*tv.nodal_volumes)(j)+=volume;(*tv.nodal_volumes)(k)+=volume;(*tv.nodal_volumes)(l)+=volume;}}
    else
        for(int t=1;t<=tv.mesh.elements.m;t++){int i,j,k,l;tv.mesh.elements(t).Get(i,j,k,l);
            T volume=(T).25*TETRAHEDRON<T>::Signed_Volume(tv.particles.X(i),tv.particles.X(j),tv.particles.X(k),tv.particles.X(l));
            (*tv.nodal_volumes)(i)+=volume;(*tv.nodal_volumes)(j)+=volume;(*tv.nodal_volumes)(k)+=volume;(*tv.nodal_volumes)(l)+=volume;}
}
}
}
#endif
