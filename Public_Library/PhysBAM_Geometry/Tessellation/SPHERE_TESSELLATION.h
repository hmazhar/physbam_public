//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __SPHERE_TESSELLATION__
#define __SPHERE_TESSELLATION__
 
namespace PhysBAM{
template<class TV> class SPHERE;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGULATED_AREA;
template<class T> class SEGMENTED_CURVE_2D;
template<class T,int d> class VECTOR;

namespace TESSELLATION{
//#####################################################################
template<class T> TRIANGULATED_AREA<T>* Generate_Triangles(const SPHERE<VECTOR<T,2> >& circle,int levels=4);
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const SPHERE<VECTOR<T,3> >& sphere,int levels=4);
template<class T> TRIANGULATED_SURFACE<T>* Tessellate_Boundary(const SPHERE<VECTOR<T,3> >& sphere,int levels=4)
{return Generate_Triangles(sphere,levels);}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const SPHERE<VECTOR<T,2> >& sphere,int levels=4);
//#####################################################################
}
}
#endif
