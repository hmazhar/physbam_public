//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __ORIENTED_BOX_TESSELLATION__
#define __ORIENTED_BOX_TESSELLATION__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{
template<class TV> class ORIENTED_BOX;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const ORIENTED_BOX<VECTOR<T,3> >& box);
//#####################################################################
}
}
#endif
