//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_AREA_REFRESH__
#define __TRIANGULATED_AREA_REFRESH__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Funcion Initialize_Square_Mesh_And_Particles
//#####################################################################
template<class T,class TV>
void Initialize_Square_Mesh_And_Particles(TRIANGULATED_AREA<T>& ta,const GRID<TV>& grid,const bool reverse_triangles)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    ta.particles.array_collection->Delete_All_Elements();ta.mesh.Initialize_Square_Mesh(m,n,reverse_triangles);ta.particles.array_collection->Add_Elements(m*n);
    for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) ta.particles.X(++particle)=grid.X(i,j);
}
//#####################################################################
// Funcion Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Circle_Mesh_And_Particles(TRIANGULATED_AREA<T>& ta,const T outer_radius,const T inner_radius,const int num_radial,const int num_tangential)
{
    int particle=0;ta.particles.array_collection->Delete_All_Elements();
    ta.mesh.Initialize_Circle_Mesh(num_radial,num_tangential);ta.particles.array_collection->Add_Elements(num_radial*num_tangential);
    for(int j=1;j<=num_radial;j++) for(int i=1;i<=num_tangential;i++){
        T r=T(j-1)/T(num_tangential-1)*(outer_radius-inner_radius)+inner_radius,theta=T(i)/T(num_tangential)*(T)2*(T)pi; 
        ta.particles.X(++particle)=VECTOR<T,2>(r*cos(theta),r*sin(theta));}
}
//#####################################################################
// Funcion Initialize_Herringbone_Mesh_And_Particles
//#####################################################################
template<class T,class TV>
void Initialize_Herring_Bone_Mesh_And_Particles(TRIANGULATED_AREA<T>& ta,const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    ta.particles.array_collection->Delete_All_Elements();ta.mesh.Initialize_Herring_Bone_Mesh(m,n);ta.particles.array_collection->Add_Elements(m*n);
    for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) ta.particles.X(++particle)=grid.X(i,j);
}
//#####################################################################
// Funcion Initialize_Equilateral_Mesh_And_Particles
//#####################################################################
template<class T,class TV>
void Initialize_Equilateral_Mesh_And_Particles(TRIANGULATED_AREA<T>& ta,const GRID<TV>& grid)
{
    int m=grid.counts.x,n=grid.counts.y,particle=0;
    ta.particles.array_collection->Delete_All_Elements();ta.mesh.Initialize_Equilateral_Mesh(m,n);ta.particles.array_collection->Add_Elements(m*n);
    for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) ta.particles.X(++particle)=VECTOR<T,2>((j%2)?grid.Axis_X_minus_half(i,1):grid.Axis_X(i,1),grid.Axis_X(j,2));
}
//#####################################################################
// Function Initialize_Segmented_Curve
//#####################################################################
template<class T>
void Initialize_Segmented_Curve(TRIANGULATED_AREA<T>& ta)
{
    delete ta.segmented_curve;if(!ta.mesh.boundary_mesh) ta.mesh.Initialize_Boundary_Mesh();
    ta.segmented_curve=new SEGMENTED_CURVE_2D<T>(*ta.mesh.boundary_mesh,ta.particles);
}
//#####################################################################
// Function Initialize_Triangle_Area_Fractions_From_Voronoi_Regions
//#####################################################################
template<class T>
void Initialize_Triangle_Area_Fractions_From_Voronoi_Regions(TRIANGULATED_AREA<T>& ta)
{
    if(!ta.triangle_area_fractions) ta.triangle_area_fractions=new ARRAY<VECTOR<T,2> >;
    ta.triangle_area_fractions->Resize(ta.mesh.elements.m);
    for(int t=1;t<=ta.mesh.elements.m;t++){
        int i,j,k;ta.mesh.elements(t).Get(i,j,k);
        VECTOR<T,3> fractions=VECTOR<T,3>(1,1,1)-TRIANGLE_2D<T>::Circumcenter_Barycentric_Coordinates(ta.particles.X(i),ta.particles.X(j),ta.particles.X(k));
        fractions.x=max(T(0),fractions.x);fractions.y=max(T(0),fractions.y);fractions.z=max(T(0),fractions.z);
        fractions/=(fractions.x+fractions.y+fractions.z);
        (*ta.triangle_area_fractions)(t).Set(fractions.x,fractions.y);}
}
//#####################################################################
// Function Compute_Triangle_Areas
//#####################################################################
template<class T>
void Compute_Triangle_Areas(TRIANGULATED_AREA<T>& ta)
{
    if(!ta.triangle_areas) ta.triangle_areas=new ARRAY<T>;
    ta.triangle_areas->Resize(ta.mesh.elements.m);
    for(int t=1;t<=ta.mesh.elements.m;t++){
        int i,j,k;ta.mesh.elements(t).Get(i,j,k);
        (*ta.triangle_areas)(t)=TRIANGLE_2D<T>::Signed_Area(ta.particles.X(i),ta.particles.X(j),ta.particles.X(k));}
}
//#####################################################################
// Function Compute_Nodal_Areas
//#####################################################################
template<class T>
void Compute_Nodal_Areas(TRIANGULATED_AREA<T>& ta,bool save_triangle_areas)
{
    if(!ta.nodal_areas) ta.nodal_areas=new ARRAY<T>;
    ta.nodal_areas->Resize(ta.particles.array_collection->Size());
    if(save_triangle_areas){
        if(!ta.triangle_areas) ta.triangle_areas=new ARRAY<T>;
        ta.triangle_areas->Resize(ta.mesh.elements.m);}
    ARRAYS_COMPUTATIONS::Fill(*ta.nodal_areas,(T)0);
    for(int t=1;t<=ta.mesh.elements.m;t++){
        int i,j,k;ta.mesh.elements(t).Get(i,j,k);
        T area=TRIANGLE_2D<T>::Signed_Area(ta.particles.X(i),ta.particles.X(j),ta.particles.X(k));
        if(save_triangle_areas) (*ta.triangle_areas)(t)=area;
        if(!ta.triangle_area_fractions){area*=(T)one_third;(*ta.nodal_areas)(i)+=area;(*ta.nodal_areas)(j)+=area;(*ta.nodal_areas)(k)+=area;}
        else{T wi,wj;(*ta.triangle_area_fractions)(t).Get(wi,wj);(*ta.nodal_areas)(i)+=wi*area;(*ta.nodal_areas)(j)+=wj*area;(*ta.nodal_areas)(k)+=(1-wi-wj)*area;}}
}
}
}
#endif
