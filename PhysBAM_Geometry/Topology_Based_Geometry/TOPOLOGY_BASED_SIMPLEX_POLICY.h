//#####################################################################
// Copyright 2006-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOPOLOGY_BASED_SIMPLEX_POLICY 
//#####################################################################
#ifndef __TOPOLOGY_BASED_SIMPLEX_POLICY__
#define __TOPOLOGY_BASED_SIMPLEX_POLICY__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;

template<class TV,int d> struct TOPOLOGY_BASED_SIMPLEX_POLICY;

//#####################################################################
// Points
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<T,1>,0>
{
    typedef VECTOR<T,1> TV;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::POINT_SIMPLICES OBJECT;
    typedef PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > HIERARCHY;
};
template<class TV>
struct TOPOLOGY_BASED_SIMPLEX_POLICY<TV,0>
{
    typedef PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > HIERARCHY;
};
//#####################################################################
// Segments
//#####################################################################
template<class TV>
struct TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>
{
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE OBJECT;
    typedef SEGMENT_HIERARCHY<TV> HIERARCHY;
};
//#####################################################################
// Triangles
//#####################################################################
template<class TV>
struct TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>
{
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT OBJECT;
    typedef TRIANGLE_HIERARCHY<typename TV::SCALAR> HIERARCHY;
};
//#####################################################################
// Tetrahedra
//#####################################################################
template<class T>
struct TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<T,3>,3>
{
    typedef TETRAHEDRALIZED_VOLUME<T> OBJECT;
    typedef TETRAHEDRON_HIERARCHY<T> HIERARCHY;
};
//#####################################################################
}
#endif
