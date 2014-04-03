//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_AREA_PRUNE__
#define __TRIANGULATED_AREA_PRUNE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Discard_Triangles_Outside_Implicit_Curve
//#####################################################################
// uses Whitney-like criterion to discard only those triangles that are for sure outside levelset (assuming accurate signed distance)
template<class T,class TV> void
Discard_Triangles_Outside_Implicit_Curve(const TRIANGULATED_AREA<T>& ta,IMPLICIT_OBJECT<TV>& implicit_curve)
{
    VECTOR<T,2> xi,xj,xk;int t=1;
    while(t <= ta.mesh.elements.m){
        int i,j,k;ta.mesh.elements(t).Get(i,j,k);
        xi=ta.particles.X(i);xj=ta.particles.X(j);xk=ta.particles.X(k);
        T max_length=sqrt(max((xi-xj).Magnitude_Squared(),(xj-xk).Magnitude_Squared(),(xk-xi).Magnitude_Squared()));
        T min_phi=min(implicit_curve(xi),implicit_curve(xj),implicit_curve(xk));
        if(min_phi > max_length) ta.mesh.elements.Remove_Index_Lazy(t);else t++;}
}
}
}
#endif
