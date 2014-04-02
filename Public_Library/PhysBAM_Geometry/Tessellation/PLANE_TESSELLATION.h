//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __PLANE_TESSELLATION__
#define __PLANE_TESSELLATION__
 
namespace PhysBAM{
template<class T> class PLANE;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const PLANE<T>& plane);
//#####################################################################
}
}
#endif
