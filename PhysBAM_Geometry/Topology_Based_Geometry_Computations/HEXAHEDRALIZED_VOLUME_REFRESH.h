//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __HEXAHEDRALIZED_VOLUME_REFRESH__
#define __HEXAHEDRALIZED_VOLUME_REFRESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>

namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Update_Hexahedron_List
//#####################################################################
template<class T>
void Update_Hexahedron_List(HEXAHEDRALIZED_VOLUME<T>& hv)
{
    if(!hv.hexahedron_list) hv.hexahedron_list=new ARRAY<HEXAHEDRON<T> >(hv.mesh.elements.m);
    for(int t=1;t<=hv.mesh.elements.m;t++){
        int p1,p2,p3,p4,p5,p6,p7,p8;hv.mesh.elements(t).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        (*hv.hexahedron_list)(t).x1=hv.particles.X(p1);(*hv.hexahedron_list)(t).x2=hv.particles.X(p2);
        (*hv.hexahedron_list)(t).x3=hv.particles.X(p3);(*hv.hexahedron_list)(t).x4=hv.particles.X(p4);
        (*hv.hexahedron_list)(t).x5=hv.particles.X(p5);(*hv.hexahedron_list)(t).x6=hv.particles.X(p6);
        (*hv.hexahedron_list)(t).x7=hv.particles.X(p7);(*hv.hexahedron_list)(t).x8=hv.particles.X(p8);}
}
//#####################################################################
// Funcion Initialize_Tetrahedralized_Volume
//#####################################################################
template<class T>
void Initialize_Tetrahedralized_Volume(HEXAHEDRALIZED_VOLUME<T>& hv)
{
    typedef VECTOR<T,3> TV;
    hv.mesh.Initialize_Faces();hv.mesh.Initialize_Face_Hexahedrons();
    ARRAY<int> face_particle_indices(hv.mesh.faces->m);ARRAY<int> hex_particle_indices(hv.mesh.elements.m);ARRAY<VECTOR<int,4> > tetrahedron_list;
    //add node in the center of each hex
    for(int h=1;h<=hv.mesh.elements.m;h++){
        ARRAY<int> p(8);hv.mesh.elements(h).Get(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
        TV hex_center;for(int i=1;i<=8;i++) hex_center+=hv.particles.X(p(i));hex_center*=(T).125;
        hv.particles.X(hv.particles.array_collection->Add_Element())=hex_center;hex_particle_indices(h)=hv.particles.array_collection->Size();}
    //add node in the center of each boundary face
    for(int f=1;f<=hv.mesh.faces->m;f++){
        int h=(*hv.mesh.face_hexahedrons)(f)(2),node1,node2,node3,node4;(*hv.mesh.faces)(f).Get(node1,node2,node3,node4);
        if(!h){hv.particles.X(hv.particles.array_collection->Add_Element())=(T).25*(hv.particles.X(node1)+hv.particles.X(node2)+hv.particles.X(node3)+hv.particles.X(node4));face_particle_indices(f)=hv.particles.array_collection->Size();}}
    //for each face, add in four tets from the associated octahedron
    for(int f=1;f<=hv.mesh.faces->m;f++){
        int h_outward,h_inward,h1,h2,p1,p2,p3,p4,ph_outward,ph_inward;(*hv.mesh.faces)(f).Get(p1,p2,p3,p4);
        //find which hexahedron the face is outwardly oriented with
        h1=(*hv.mesh.face_hexahedrons)(f)(1);h2=(*hv.mesh.face_hexahedrons)(f)(2);h_outward=h1;h_inward=h2;
        if(h2){
            h_outward=h2;h_inward=h1;
            for(int k=0;k<6;k++){// loop over faces of h1
                int i1=hv.mesh.face_indices[k][0],i2=hv.mesh.face_indices[k][1],i3=hv.mesh.face_indices[k][2],i4=hv.mesh.face_indices[k][3];
                if(p1 == hv.mesh.elements(h1)(i1) && p2 == hv.mesh.elements(h1)(i2) &&
                    p3 == hv.mesh.elements(h1)(i3) && p4 == hv.mesh.elements(h1)(i4)){h_outward=h1;h_inward=h2;break;}}}
        ph_outward=hex_particle_indices(h_outward);if(h_inward == 0) ph_inward=face_particle_indices(f); else ph_inward=hex_particle_indices(h_inward);
        tetrahedron_list.Append(VECTOR<int,4>(p1,ph_inward,p2,ph_outward));tetrahedron_list.Append(VECTOR<int,4>(p2,ph_inward,p3,ph_outward));
        tetrahedron_list.Append(VECTOR<int,4>(p3,ph_inward,p4,ph_outward));tetrahedron_list.Append(VECTOR<int,4>(p4,ph_inward,p1,ph_outward));}
    if(hv.tetrahedralized_volume) delete &(hv.tetrahedralized_volume->mesh);delete hv.tetrahedralized_volume;
    hv.tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(hv.particles);
    hv.tetrahedralized_volume->mesh.Initialize_Mesh(hv.particles.array_collection->Size(),tetrahedron_list);
}
//#####################################################################
// Funcion Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T,class TV>
void Initialize_Cube_Mesh_And_Particles(HEXAHEDRALIZED_VOLUME<T>& hv,const GRID<TV>& grid)
{
    hv.particles.array_collection->Delete_All_Elements();
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    hv.mesh.Initialize_Cube_Mesh(m,n,mn);
    hv.particles.array_collection->Preallocate(m*n*mn);
    for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++)for(int i=1;i<=m;i++) hv.particles.X(hv.particles.array_collection->Add_Element())=grid.X(i,j,ij);
}
}
}
#endif
