//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOPOLOGY_BASED_GEOMETRY_POLICY 
//#####################################################################
#ifndef __TOPOLOGY_BASED_GEOMETRY_POLICY__
#define __TOPOLOGY_BASED_GEOMETRY_POLICY__

#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> struct TOPOLOGY_BASED_GEOMETRY_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,0> >
{
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,1> >
{
    struct UNUSABLE{};
    typedef POINT_SIMPLICES_1D<T> POINT_SIMPLICES;
    typedef PhysBAM::SEGMENTED_CURVE<VECTOR<T,1> > SEGMENTED_CURVE;
    typedef UNUSABLE TRIANGULATED_OBJECT;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,2> >
{
    typedef SEGMENTED_CURVE_2D<T> SEGMENTED_CURVE;
    typedef TRIANGULATED_AREA<T> TRIANGULATED_OBJECT;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T,3> >
{
    typedef PhysBAM::SEGMENTED_CURVE<VECTOR<T,3> > SEGMENTED_CURVE;
    typedef TRIANGULATED_SURFACE<T> TRIANGULATED_OBJECT;
};
//#####################################################################
}
#endif
