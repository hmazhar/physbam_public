//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SEGMENTED_CURVE_REFRESH__
#define __SEGMENTED_CURVE_REFRESH__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Update_Segment_List
//#####################################################################
template<class TV>
void Update_Segment_List(SEGMENTED_CURVE<TV>& sc) // updates the segments assuming the particle positions are already updated
{
    if(!sc.segment_list) sc.segment_list=new ARRAY<typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT>(sc.mesh.elements.m);
    for(int s=1;s<=sc.mesh.elements.m;s++){
        int i,j;sc.mesh.elements(s).Get(i,j);
        (*sc.segment_list)(s)=typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT(sc.particles.X(i),sc.particles.X(j));}
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class TV>
void Initialize_Hierarchy(SEGMENTED_CURVE<TV>& sc,const bool update_boxes) // creates and updates the boxes as well
{
    delete sc.hierarchy;
    if(sc.segment_list) sc.hierarchy=new SEGMENT_HIERARCHY<TV>(sc.mesh,sc.particles,*sc.segment_list,update_boxes);
    else sc.hierarchy=new SEGMENT_HIERARCHY<TV>(sc.mesh,sc.particles,update_boxes);
}
//#####################################################################
// Function Initialize_Straight_Mesh_And_Particles
//#####################################################################
template<class TV,class T>
void Initialize_Straight_Mesh_And_Particles(SEGMENTED_CURVE<TV>& sc,const GRID<VECTOR<T,1> >& grid)
{
    sc.Clean_Memory();
    sc.mesh.Initialize_Straight_Mesh(grid.counts.x);
    sc.particles.array_collection->Preallocate(grid.counts.x);
    for(int i=1;i<=grid.counts.x;i++) sc.particles.X(sc.particles.array_collection->Add_Element()).x=grid.X(VECTOR<int,1>(i)).x;
}
//#####################################################################
// Function Initialize_Circle_Mesh_And_Particles
//#####################################################################
template<class TV,class T>
void Initialize_Circle_Mesh_And_Particles(SEGMENTED_CURVE<TV>& sc,const int m,const T radius)
{
    sc.Clean_Memory();
    sc.mesh.Initialize_Straight_Mesh(m,true);
    sc.particles.array_collection->Add_Elements(m);
    for(int p=1;p<=m;p++){
        COMPLEX<T> X=COMPLEX<T>::Polar(radius,(T)2*(T)pi/m*p);
        sc.particles.X(p)=TV(VECTOR<T,2>(X.re,X.im));}
}
}
}
#endif
