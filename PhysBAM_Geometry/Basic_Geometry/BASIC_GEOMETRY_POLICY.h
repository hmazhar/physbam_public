//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_GEOMETRY_POLICY 
//#####################################################################
#ifndef __BASIC_GEOMETRY_POLICY__
#define __BASIC_GEOMETRY_POLICY__

#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> struct BASIC_GEOMETRY_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct BASIC_GEOMETRY_POLICY<VECTOR<T,0> >
{
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct BASIC_GEOMETRY_POLICY<VECTOR<T,1> >
{
    typedef VECTOR<T,1> POINT;
    typedef BOX<VECTOR<T,1> > ORIENTED_BOX;
    typedef POINT_SIMPLEX_1D<T> HYPERPLANE;
    typedef POINT_SIMPLEX_1D<T> POINT_SIMPLEX;
    typedef SEGMENT_1D<T> SEGMENT;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct BASIC_GEOMETRY_POLICY<VECTOR<T,2> >
{
    typedef POINT_2D<T> POINT;
    typedef PhysBAM::ORIENTED_BOX<VECTOR<T,2> > ORIENTED_BOX;
    typedef LINE_2D<T> HYPERPLANE;
    typedef SEGMENT_2D<T> SEGMENT;
    typedef TRIANGLE_2D<T> TRIANGLE;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct BASIC_GEOMETRY_POLICY<VECTOR<T,3> >
{
    typedef VECTOR<T,3> POINT;
    typedef PhysBAM::ORIENTED_BOX<VECTOR<T,3> > ORIENTED_BOX;
    typedef PLANE<T> HYPERPLANE;
    typedef SEGMENT_3D<T> SEGMENT;
    typedef TRIANGLE_3D<T> TRIANGLE;
};
//#####################################################################
}
#endif
