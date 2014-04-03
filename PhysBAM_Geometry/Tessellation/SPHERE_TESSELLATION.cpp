//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const SPHERE<VECTOR<T,3> >& sphere,int levels)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(6);
    particles.X(1)=TV(-1,0,0);particles.X(2)=TV(1,0,0);particles.X(3)=TV(0,-1,0);
    particles.X(4)=TV(0,1,0);particles.X(5)=TV(0,0,-1);particles.X(6)=TV(0,0,1);
    ARRAY<VECTOR<int,3> >& triangles=surface->mesh.elements;triangles.Exact_Resize(8);
    triangles(1).Set(1,6,4);triangles(2).Set(1,3,6);triangles(3).Set(6,2,4);triangles(4).Set(6,3,2);
    triangles(5).Set(5,1,4);triangles(6).Set(5,3,1);triangles(7).Set(2,3,5);triangles(8).Set(2,5,4);
    surface->mesh.number_nodes=6;
    surface->mesh.Initialize_Neighbor_Nodes();
    for(int i=1;i<=levels;i++) surface->Root_Three_Subdivide();
    for(int p=1;p<=particles.array_collection->Size();p++) particles.X(p)=sphere.center+sphere.radius*particles.X(p).Normalized();
    return surface;
}
template<class T> TRIANGULATED_AREA<T>* Generate_Triangles(const SPHERE<VECTOR<T,2> >& circle,int levels)
{
    assert(levels>=1);
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> E;
    TRIANGULATED_AREA<T>* area=TRIANGULATED_AREA<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=area->particles;particles.array_collection->Add_Elements(1+3*levels*(levels+1));
    particles.X(1)=TV(0,0);
    for(int i=1,k=1;i<=levels;i++)
        for(int j=0;j<i*6;j++){
            T a=(T)pi/3*(j-(T).5*(i-1))/i;
            particles.X(++k)=(T)i*TV(cos(a),sin(a));}
    for(int i=2;i<=7;i++) area->mesh.elements.Append(E(1,i,((i-1)%6+2)));
    for(int i=2,p=2;i<=levels;i++){
        int n=i*6;
        int c=p+n-6;
        area->mesh.elements.Append(E(c-1,c+n-1,c));
        area->mesh.elements.Append(E(c,p,c-1));
        for(int j=0;j<n-1;j++){
            area->mesh.elements.Append(E(p,c+j,c+j+1));
            if(j%i!=(i-1)/2){
                area->mesh.elements.Append(E(c+j+1,p+1,p));
                p++;}}
        p=c;}
    particles.X=particles.X*(circle.radius/levels)+circle.center;
    area->Update_Number_Nodes();
    return area;
}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const SPHERE<VECTOR<T,2> >& sphere,int levels)
{
    assert(levels>=1);
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    int n=1<<levels;
    curve->particles.array_collection->Add_Elements(n);
    for(int i=1;i<=n;i++) curve->particles.X(i)=VECTOR<T,2>(cos(i*two_pi/n),sin(i*two_pi/n))*sphere.radius+sphere.center;
    for(int i=1;i<n;i++) curve->mesh.elements.Append(VECTOR<int,2>(i,i+1));
    curve->mesh.elements.Append(VECTOR<int,2>(n,1));
    curve->Update_Number_Nodes();
    return curve;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const SPHERE<VECTOR<float,3> >&,int);
template TRIANGULATED_AREA<float>* Generate_Triangles(const SPHERE<VECTOR<float,2> >&,int);
template SEGMENTED_CURVE_2D<float>* Tessellate_Boundary(const SPHERE<VECTOR<float,2> >& sphere,int levels);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const SPHERE<VECTOR<double,3> >&,int);
template TRIANGULATED_AREA<double>* Generate_Triangles(const SPHERE<VECTOR<double,2> >&,int);
template SEGMENTED_CURVE_2D<double>* Tessellate_Boundary(const SPHERE<VECTOR<double,2> >& sphere,int levels);
#endif
}
}
