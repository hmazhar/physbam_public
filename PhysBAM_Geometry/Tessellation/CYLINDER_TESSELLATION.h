//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __CYLINDER_TESSELLATION__
#define __CYLINDER_TESSELLATION__
 
namespace PhysBAM{
template<class T> class CYLINDER;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const CYLINDER<T>& cylinder,const int radius_height=4,const int resolution_radius=16);
//#####################################################################
}
}
#endif
