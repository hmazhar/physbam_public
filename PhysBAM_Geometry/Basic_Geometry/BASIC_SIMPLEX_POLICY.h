//#####################################################################
// Copyright 2006-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_SIMPLEX_POLICY 
//#####################################################################
#ifndef __BASIC_SIMPLEX_POLICY__
#define __BASIC_SIMPLEX_POLICY__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
namespace PhysBAM{

template<class T> class TETRAHEDRON;

template<class TV,int d> struct BASIC_SIMPLEX_POLICY;

//#####################################################################
// Points
//#####################################################################
template<class T>
struct BASIC_SIMPLEX_POLICY<VECTOR<T,1>,0>
{
    typedef VECTOR<T,1> TV;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::POINT_SIMPLEX SIMPLEX_FACE;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::POINT_SIMPLEX SIMPLEX;
};
template<class TV>
struct BASIC_SIMPLEX_POLICY<TV,0>
{
};
//#####################################################################
// Segments
//#####################################################################
template<class TV>
struct BASIC_SIMPLEX_POLICY<TV,1>
{
    typedef typename BASIC_GEOMETRY_POLICY<TV>::POINT SIMPLEX_FACE;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT SIMPLEX;
};
//#####################################################################
// Triangles
//#####################################################################
template<class TV>
struct BASIC_SIMPLEX_POLICY<TV,2>
{
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT SIMPLEX_FACE;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::TRIANGLE SIMPLEX;
};
//#####################################################################
// Tetrahedra
//#####################################################################
template<class T>
struct BASIC_SIMPLEX_POLICY<VECTOR<T,3>,3>
{
    typedef typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::TRIANGLE SIMPLEX_FACE;
    typedef TETRAHEDRON<T> SIMPLEX;
};
//#####################################################################
}
#endif
