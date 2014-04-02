//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MESH_POLICY 
//#####################################################################
#ifndef __MESH_POLICY__
#define __MESH_POLICY__

namespace PhysBAM{

class SEGMENT_MESH;
class TRIANGLE_MESH;
class TETRAHEDRON_MESH;
class HEXAHEDRON_MESH;
class POINT_SIMPLEX_MESH;

template<int d> struct MESH_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<>
struct MESH_POLICY<0>
{
    typedef POINT_SIMPLEX_MESH MESH;
};
//#####################################################################
// 1D
//#####################################################################
template<>
struct MESH_POLICY<1>
{
    typedef SEGMENT_MESH MESH;
};
//#####################################################################
// 2D
//#####################################################################
template<>
struct MESH_POLICY<2>
{
    typedef TRIANGLE_MESH MESH;
};
//#####################################################################
// 3D
//#####################################################################
template<>
struct MESH_POLICY<3>
{
    typedef TETRAHEDRON_MESH MESH;
};
//#####################################################################
}
#endif
