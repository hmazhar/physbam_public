//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_GEOMETRY_FORWARD 
//#####################################################################
#ifndef __BASIC_GEOMETRY_FORWARD__
#define __BASIC_GEOMETRY_FORWARD__

namespace PhysBAM{

template<class TV> class BOX;
template<class TV> class RANGE;
template<class TV> class RAY;
template<class T,int d> class VECTOR;
template<class TV> class ORIENTED_BOX;
template<class T> class POINT_SIMPLEX_1D;
template<class T> class LINE_2D;
template<class T> class PLANE;
template<class T> class SEGMENT_1D;
template<class T> class SEGMENT_2D;
template<class T> class SEGMENT_3D;
template<class T> class TRIANGLE_2D;
template<class T> class TRIANGLE_3D;
template<class T> class TETRAHEDRON;
template<class T> class CYLINDER;
template<class T> class RING;
template<class T> class ELLIPSOID;
template<class TV> class SPHERE;
template<class T> class TORUS;
template<class T> class POINT_2D;
template<class TV> class BOUNDED_HORIZONTAL_PLANE;

enum POINT_SIMPLEX_COLLISION_TYPE {POINT_SIMPLEX_NO_COLLISION,POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE,POINT_SIMPLEX_COLLISION_ENDS_INSIDE,POINT_SIMPLEX_UNKNOWN_COLLISION};

void Initialize_Read_Write_Structures();

}
#endif
