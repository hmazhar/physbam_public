//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __IMPLICIT_OBJECT_TESSELLATION__
#define __IMPLICIT_OBJECT_TESSELLATION__
 
namespace PhysBAM{
template<class TV> class IMPLICIT_OBJECT;
template<class T> class POINT_SIMPLICES_1D;
template<class T> class SEGMENTED_CURVE_2D;
template<class T> class TRIANGULATED_SURFACE;
template<class T,int d> class VECTOR;

namespace TESSELLATION{
//#####################################################################
    template<class T> POINT_SIMPLICES_1D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,1> >& implicit);
    template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,2> >& implicit);
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,3> >& implicit);
//#####################################################################
}
}
#endif
