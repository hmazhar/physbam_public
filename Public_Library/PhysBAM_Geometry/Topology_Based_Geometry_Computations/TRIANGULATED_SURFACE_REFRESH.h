//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_REFRESH__
#define __TRIANGULATED_SURFACE_REFRESH__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T>
void Initialize_Hierarchy(TRIANGULATED_SURFACE<T>& ts,const bool update_boxes,const int triangles_per_group) // creates and updates the boxes as well
{
    delete ts.hierarchy;
    if(ts.triangle_list) ts.hierarchy=new TRIANGLE_HIERARCHY<T>(ts.mesh,ts.particles,*ts.triangle_list,update_boxes,triangles_per_group);
    else ts.hierarchy=new TRIANGLE_HIERARCHY<T>(ts.mesh,ts.particles,update_boxes,triangles_per_group);
}
//#####################################################################
// Function Initialize_Particle_Hierarchy
//#####################################################################
template<class T>
void Initialize_Particle_Hierarchy(TRIANGULATED_SURFACE<T>& ts,const INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,3> > >& array_input,const bool update_boxes,const int particles_per_group) // creates and updates the boxes as well
{
    typedef VECTOR<T,3> TV;
    delete ts.particle_hierarchy;
    ts.particle_hierarchy=new PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >(array_input,update_boxes,particles_per_group);
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T>
void Update_Triangle_List(TRIANGULATED_SURFACE<T>& ts,ARRAY_VIEW<const VECTOR<T,3> > X)
{
    if(!ts.triangle_list) ts.triangle_list=new ARRAY<TRIANGLE_3D<T> >;
    ts.triangle_list->Resize(ts.mesh.elements.m);
    for(int t=1;t<=ts.mesh.elements.m;t++)
        (*ts.triangle_list)(t)=TRIANGLE_3D<T>(X.Subset(ts.mesh.elements(t)));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T>
void Rescale(TRIANGULATED_SURFACE<T>& ts,const T scaling_x,const T scaling_y,const T scaling_z)
{
    typedef VECTOR<T,3> TV;
    if(scaling_x*scaling_y*scaling_z<=0) PHYSBAM_FATAL_ERROR();
    for(int k=1;k<=ts.particles.array_collection->Size();k++) ts.particles.X(k)*=TV(scaling_x,scaling_y,scaling_z);
    if(ts.triangle_list) ts.Update_Triangle_List();if(ts.hierarchy) ts.hierarchy->Update_Boxes();if(ts.bounding_box) ts.Update_Bounding_Box();
}
//#####################################################################
// Function Initialize_Segment_Lengths
//#####################################################################
template<class T>
void Initialize_Segment_Lengths(TRIANGULATED_SURFACE<T>& ts)
{
    bool segment_mesh_defined=ts.mesh.segment_mesh!=0;if(!segment_mesh_defined) ts.mesh.Initialize_Segment_Mesh();
    delete ts.segment_lengths;ts.segment_lengths=new ARRAY<T>(ts.mesh.segment_mesh->elements.m);
    for(int t=1;t<=ts.mesh.segment_mesh->elements.m;t++) 
        (*ts.segment_lengths)(t)=(ts.particles.X(ts.mesh.segment_mesh->elements(t)(1))-ts.particles.X(ts.mesh.segment_mesh->elements(t)(2))).Magnitude();
    if(!segment_mesh_defined){delete ts.mesh.segment_mesh;ts.mesh.segment_mesh=0;}
}
//#####################################################################
// Function Initialize_Torus_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Torus_Mesh_And_Particles(TRIANGULATED_SURFACE<T>& ts,const int m,const int n,const T major_radius,const T minor_radius)
{
    typedef VECTOR<T,3> TV;
    T di=T(2*pi)/m,dj=T(2*pi)/n;
    for(int j=1;j<=n;j++){
        T phi=-dj*j,radius=major_radius+minor_radius*cos(phi),z=minor_radius*sin(phi);
        for(int i=1;i<=m;i++){
            int p=ts.particles.array_collection->Add_Element();T theta=di*(i-(T).5*(j&1));
            ts.particles.X(p)=TV(radius*cos(theta),radius*sin(theta),z);}}
    ts.mesh.Initialize_Torus_Mesh(m,n);
}
//#####################################################################
// Function Initialize_Cylinder_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Cylinder_Mesh_And_Particles(TRIANGULATED_SURFACE<T>& ts,const int m,const int n,const T length,const T radius,const bool create_caps)
{
    typedef VECTOR<T,3> TV;
    ts.particles.array_collection->Delete_All_Elements();T dtheta=(T)two_pi/n;T dlength=length/(m-1);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        int p=ts.particles.array_collection->Add_Element();T theta=(j-1)*dtheta;
        ts.particles.X(p)=TV(dlength*(i-1),radius*sin(theta),radius*cos(theta));}
    if(create_caps){int p_1=ts.particles.array_collection->Add_Element();int p_2=ts.particles.array_collection->Add_Element();ts.particles.X(p_1)=TV(0,0,0);ts.particles.X(p_2)=TV(length,0,0);}
    ts.mesh.Initialize_Cylinder_Mesh(m,n,create_caps);
}
}
}
#endif
