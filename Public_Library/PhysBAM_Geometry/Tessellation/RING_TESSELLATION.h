//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __RING_TESSELLATION__
#define __RING_TESSELLATION__
 
namespace PhysBAM{
template<class T> class RING;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RING<T>& ring,const int n=40);
//#####################################################################
}
}
#endif
