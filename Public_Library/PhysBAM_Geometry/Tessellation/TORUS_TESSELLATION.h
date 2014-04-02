//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __TORUS_TESSELLATION__
#define __TORUS_TESSELLATION__
 
namespace PhysBAM{
template<class T> class TORUS;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const TORUS<T>& torus,const int inner_resolution=8,const int outer_resolution=16);
//#####################################################################
}
}
#endif
